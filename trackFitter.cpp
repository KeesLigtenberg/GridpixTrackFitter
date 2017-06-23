/*
 * trackFitter.cpp
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#include "trackFitter.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"

using namespace std;

trackFitter::trackFitter(std::string inputfile, const DetectorConfiguration& detector) :
	detector(detector),
	houghTransform(detector.planexmax(), detector.planeymax(), 50/*xbins*/, 25 /*ybins*/ ),
	residualHistograms(nullptr)
{

	//get tree from file
	file=openFile(inputfile);
	hitTable=getObjectFromFile<TTree>("Hits", file);

	//setup tree for reading
	hitTable->SetBranchAddress("mimosa", &mimosaHit);
	//    unsigned short triggerNumberBegin, triggerNumberEnd;
//	hitTable->SetBranchAddress("triggerNumberBegin", &triggerNumberBegin);
//	hitTable->SetBranchAddress("triggerNumberEnd", &triggerNumberEnd);

}

trackFitter::~trackFitter() {
	file->Close();
}

int trackFitter::makeMask(double ntimesThreshold) {
	//first mask noisy pixels
	for(int iplane=0; iplane<detector.nPlanes; iplane++) {
		mask.emplace_back( makeNoisyPixelMask(hitTable,iplane,ntimesThreshold,{detector.pixelColumns, detector.pixelRows} ) );
	}
	return 0; //should return number of pixels masked
}

bool trackFitter::passEvent(std::vector<std::vector<PositionHit> > spaceHit) {
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

void trackFitter::fitTracks(std::string outputfilename) {
	residualHistograms=unique_ptr(new ResidualHistogrammer(outputfilename));

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
			residualHistograms->fill(residuals);
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


	auto means=residualHistograms->getMeansOfPlanes();
	auto rotation=residualHistograms->getRotationOfPlanes();

	for(auto& m : means ) cout<<m.first<<" "<<m.second<<endl;
}
