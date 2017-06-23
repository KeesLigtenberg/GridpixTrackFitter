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


void FitTracks (std::string inputfile) {



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
