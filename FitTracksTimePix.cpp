/*
 * FitTracksTimePix.cpp
 *
 *  Created on: 11 jul. 2017
 *      Author: KeesL
 */

#include <iostream>

#include "TTree.h"

#include "PositionHit.h"
#include "HoughTransformer.h"
#include "linearRegressionFit.h"
#include "Hit.h"
#include "getObjectFromFile.h"
#include "ResidualHistogrammer.h"
#include "makeNoisyPixelMask.h"

#if 1 //root?
#include "linearRegressionFit.cpp"
#include "ResidualHistogrammer.cpp"
#include "makeNoisyPixelMask.cpp"
#endif

using namespace std;

struct TimePixDetectorConfiguration : DetectorConfiguration {
	TimePixDetectorConfiguration() : DetectorConfiguration{
		1, {0}, //nplanes, planeposition
		1, 256, 256 //pixelsize, xpixels, ypixels
	} {};
	virtual double xmin() const {return -200; }
	virtual double xmax() const {return 200; }
	virtual double zmin() const {return 0; }
	virtual double zmax() const {return 256; };
} timePixChip;

void FitTracksTimePix(std::string inputfile) {

	//get tree from file
	TFile* file=openFile(inputfile);
	TTree* hitTable=getObjectFromFile<TTree>("Hits", file);

	auto mask=makeNoisyPixelMask(hitTable, 2, {256,256} );

	//set adresses
	ULong64_t trigger=0;
	const std::vector<TimePixHit>* timepixHits=nullptr;
	hitTable->SetBranchAddress("timepix", &timepixHits );
//	hitTable->SetBranchAddress("trigger",&trigger);

	HoughTransformer houghTransform(timePixChip.xmin(),timePixChip.xmax(), timePixChip.ymin(), timePixChip.ymax(), 25/*xbins*/ ,24 /*ybins*/ );
	houghTransform.angleOfTracksX=-0.005;
	houghTransform.angleOfTracksY=-0.162;
	houghTransform.minCandidateSize=10;
	houghTransform.minClusterSize=30;

	ResidualHistogrammer residualsHistograms("residualHistogramsTimePix.root", timePixChip);

	//slopes
	double slope1Sum=0, slope2Sum=0;

	//loop over all entries
	const long long nEvents=hitTable->GetEntriesFast();
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%1000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		hitTable->GetEntry(iEvent);

		//mask
		auto maskedHits=applyPixelMask(mask, *timepixHits);

//		cout<<"got entry"<<endl;
		std::vector<PositionHit> spaceHit;
		double driftScale=25./4096;

		//convert hits changes x <-> z!!
		spaceHit=convertHits( *timepixHits, timePixChip.pixelsize, timePixChip.pixelsize, driftScale );

//		cout<<"drawing cluster"<<endl;

//		HoughTransformer::drawCluster(spaceHit, timePixChip);

		auto houghClusters = houghTransform(spaceHit);

		if(!houghClusters.size()) {
			continue;
		}


		//fit clusters
		std::vector<SimpleFitResult> fits;
		for(auto& hitCluster : houghClusters) {

			if(hitCluster.size()<2) continue;

			//fit track
			auto fit=linearRegressionFit(hitCluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}

			auto residuals = calculateResiduals(hitCluster, fit);

			residualsHistograms.fill(residuals);

			slope1Sum+=fit.slope1;
			slope2Sum+=fit.slope2;

//			std::cout<<fit<<endl;
			fits.push_back(fit);
			++nClusters;
		}

		++nPassed;

		static TCanvas* canv=nullptr;
		if(!canv) canv=new TCanvas("eventCanv", "canvas for event with fits", 600,400);
		canv->cd();
		HoughTransformer::drawCluster(spaceHit, timePixChip);
//		HoughTransformer::drawClusters(houghClusters, detector);
		for(auto& f : fits) f.draw( timePixChip.zmin(), timePixChip.zmax() );
		gPad->Update();
		if(std::cin.get()=='q') break;



	}

	cout<<"Passed "<<nPassed<<" events with "<<nClusters<<" clusters"<<endl;

	cout<<" slopes ("<<slope1Sum/nClusters<<", "<<slope2Sum/nClusters<<")"<<endl;

}


