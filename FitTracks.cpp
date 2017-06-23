#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>

#include "TH1.h"
#include "TTree.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"

#include "makeNoisyPixelMask.h"
#include "Hit.h"
#include "PositionHit.h"
#include "HoughTransformer.h"
#include "transformHits.h"

#include "DetectorConfiguration.h"

#if 0
#include "ResidualHistogrammer.h"
#include "linearRegressionFit.h"
#else
#include "ResidualHistogrammer.cpp"
#include "linearRegressionFit.cpp"
#endif

using namespace std;

const DetectorConfiguration mimosa= {
	6, //planes
	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other!
};

//return true if the event is passed
bool passEvent( std::vector<std::vector<PositionHit> > spaceHit ) {
	int nPlanesHit=0;
	int nTotalHits=0;
	for(auto& v:spaceHit) {
		int size=v.size();
		if(size) {
			nPlanesHit++;
			nTotalHits+=size;
		}
	}
	const int nMinPlanesHit=3;//require at least 4 planes hit
	return nPlanesHit>nMinPlanesHit;
}


void FitTracks (std::string inputfile) {

	//get tree from file
	TFile* file=openFile(inputfile);
	TTree* hitTable=getObjectFromFile<TTree>("Hits", file);

	//first mask noisy pixels
	std::vector<pixelMask> mask;
	for(int iplane=0; iplane<mimosa.nPlanes; iplane++) {
		const double ntimesThreshold=1e4;
		mask.emplace_back( makeNoisyPixelMask(hitTable,iplane,ntimesThreshold,{mimosa.pixelColumns, mimosa.pixelRows} ) );
	}

	//setup tree for reading
	const std::vector<std::vector<Hit>>* mimosaHit=nullptr;
	hitTable->SetBranchAddress("mimosa", &mimosaHit);
	//    unsigned short triggerNumberBegin, triggerNumberEnd;
//	hitTable->SetBranchAddress("triggerNumberBegin", &triggerNumberBegin);
//	hitTable->SetBranchAddress("triggerNumberEnd", &triggerNumberEnd);

	//setup houghtransformer
	int binsx=50,binsy=25;
	HoughTransformer houghTransform(mimosa.planexmax(), mimosa.planeymax(), binsx, binsy);

	ResidualHistogrammer residualHistogramsUnaligned("residualHistogramsUnaligned.root");

	//loop over all entries
	const long long nEvents=hitTable->GetEntriesFast(); //std::min( (long long) 2,);
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%1000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		hitTable->GetEntry(iEvent);

		std::vector<std::vector<PositionHit> > spaceHit;

		//apply mask
		auto itmask=mask.begin();
		for(unsigned plane=0; plane<mimosaHit->size(); plane++ ) {
			auto maskedHit = applyPixelMask( *itmask++ , mimosaHit->at(plane) ) ;
			//convert hits to positions
			spaceHit.push_back( convertHits(maskedHit, mimosa.planePosition[plane], mimosa.pixelsize, mimosa.pixelsize, plane) );
		}

		if( !passEvent(spaceHit) ) continue;
		++nPassed;

		//hough transform
		auto houghClusters = houghTransform(spaceHit);

		//require at least 5/6 planes to be hit
		const int nMinPlanesHit=5;
		houghClusters.remove_if([](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
		if(!houghClusters.size()) {
			continue;
		}

//		std::cout<<"found "<<houghClusters.size()<<" clusters of size ";
//		for(auto& cl : houghClusters) {
//			std::cout<<cl.clusterSize<<", ";
//		}
//		std::cout<<std::endl;
//		HoughTransformer::drawClusters(houghClusters, mimosa);


		//first fit on outer planes and translate inner hits
		std::vector<SimpleFitResult> fits;
		for(auto& hitCluster : houghClusters) {

			//remove all hits except first and last one
//			auto hitsOnOuterPlanes=hitCluster;
//			const int firstPlane=0, lastPlane=5;
//			hitsOnOuterPlanes.remove_if([](PositionHit& hit){return hit.plane != firstPlane && hit.plane != lastPlane;});

			//fit track
			auto fit=linearRegressionFit(hitCluster);
//			fit.draw(0, mimosa.planePosition.back() );

			//remove outliers
			auto residuals=calculateResiduals(hitCluster, fit);
			const double maxResidual=0.2;
			hitCluster=cutOnResiduals(hitCluster, residuals, maxResidual);

			++nClusters;
			//refit on outer planes again
			fit=linearRegressionFit(hitCluster);
			residuals=calculateResiduals(hitCluster, fit);
			residualHistogramsUnaligned.fill(residuals);
			fits.push_back(fit);
		}


//		HoughTransformer::drawClusters(spaceHit, mimosa);
////		HoughTransformer::drawClusters(houghClusters, mimosa);
//		for(auto& f : fits) f.draw(0, mimosa.planePosition.back());
//		gPad->Update();
//		auto signal=std::cin.get();
//		if(signal=='q') break;
//		else if(signal=='l') {
//			while(!gSystem->ProcessEvents()) {
//				gSystem->Sleep(50);
//			}
//			break;
//		}

	}
	std::cout<<"passed: "<<nPassed<<"/"<<nEvents<<" with "<<nClusters<<"\n";


	auto means=residualHistogramsUnaligned.getMeansOfPlanes();

	for(auto& m : means ) cout<<m.first<<" "<<m.second<<endl;


	//loop over all entries again now aligning in x,y
	ResidualHistogrammer residualHistogramsTranslated("residualHistogramsTranslated.root");

	nPassed=0;
	nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%1000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		hitTable->GetEntry(iEvent);

		std::vector<std::vector<PositionHit> > spaceHit;

		//apply mask
		auto itmask=mask.begin();
		for(unsigned plane=0; plane<mimosaHit->size(); plane++ ) {
			auto maskedHit = applyPixelMask( *itmask++ , mimosaHit->at(plane) ) ;
			//convert hits to positions
			spaceHit.push_back( convertHits(maskedHit, mimosa.planePosition[plane], mimosa.pixelsize, mimosa.pixelsize, plane) );
		}

		//apply translation
		spaceHit=translateHits(spaceHit,means);

		if( !passEvent(spaceHit) ) continue;
		++nPassed;

		//hough transform
		auto houghClusters = houghTransform(spaceHit);

		//require at least 5/6 planes to be hit
		const int nMinPlanesHit=5;
		houghClusters.remove_if([](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
		if(!houghClusters.size()) continue;

		//now fit all tracks with translated planes
		std::vector<SimpleFitResult> fits;
		for(auto& hitCluster : houghClusters) {
			++nClusters;

			//fit track
			auto fit=linearRegressionFit(hitCluster);
			auto residuals=calculateResiduals(hitCluster, fit);

			//remove outliers
			const double maxResidual=0.5;
			hitCluster=cutOnResiduals(hitCluster, residuals, maxResidual);

			//refit track
			fit=linearRegressionFit(hitCluster);
			residuals=calculateResiduals(hitCluster, fit);

			residualHistogramsTranslated.fill(residuals);
//			fit.draw(0, mimosa.planePosition.back() );
			fits.push_back(fit);
		}
	}
	std::cout<<"passed: "<<nPassed<<"/"<<nEvents<<" with "<<nClusters<<"\n";

	auto rotations=residualHistogramsTranslated.getRotationOfPlanes();
	for(auto& r : rotations ) cout<<r<<" = "<<r/M_PI*180.<<endl;


	ResidualHistogrammer residualHistograms("residualHistograms.root");
	nPassed=0;
	nClusters=0;
	//loop over all entries again now also aligning rotations
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%1000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		hitTable->GetEntry(iEvent);

		std::vector<std::vector<PositionHit> > spaceHit;

		//apply mask
		auto itmask=mask.begin();
		for(unsigned plane=0; plane<mimosaHit->size(); plane++ ) {
			auto maskedHit = applyPixelMask( *itmask++ , mimosaHit->at(plane) ) ;
			//convert hits to positions
			spaceHit.push_back( convertHits(maskedHit, mimosa.planePosition[plane], mimosa.pixelsize, mimosa.pixelsize, plane) );
		}

		//apply translation
		spaceHit=translateHits(spaceHit,means);
		spaceHit=rotateHits(spaceHit, rotations);

		if( !passEvent(spaceHit) ) continue;
		++nPassed;

		//hough transform
		auto houghClusters = houghTransform(spaceHit);

		//require at least 5/6 planes to be hit
		const int nMinPlanesHit=4;
		houghClusters.remove_if([](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
		if(!houghClusters.size()) continue;

		//now fit all tracks with translated planes
		std::vector<SimpleFitResult> fits;
		for(auto& hitCluster : houghClusters) {
			++nClusters;

			//fit track
			auto fit=linearRegressionFit(hitCluster);
			auto residuals=calculateResiduals(hitCluster, fit);

			//remove outliers
			const double maxResidual=0.2;
			hitCluster=cutOnResiduals(hitCluster, residuals, maxResidual);

			//cut on planes again
			if(hitCluster.recalculateNPlanesHit()<nMinPlanesHit) {continue;}

			//refit track
			fit=linearRegressionFit(hitCluster);
			residuals=calculateResiduals(hitCluster, fit);

			residualHistograms.fill(residuals);
//			fit.draw(0, mimosa.planePosition.back() );
			fits.push_back(fit);
		}

	}
	std::cout<<"passed: "<<nPassed<<"/"<<nEvents<<" with "<<nClusters<<"\n";

	rotations=residualHistograms.getRotationOfPlanes();
	for(auto& r : rotations ) cout<<r<<" = "<<r/M_PI*180.<<endl;

}



//redirect to root main function
int main(int argc, const char* argv[]) {
	gROOT->ProcessLine(".L Hit.h+"); //compile hit
	if(argc>1) FitTracks(argv[1]);
}
