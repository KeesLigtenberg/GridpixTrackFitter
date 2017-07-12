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

#if 1 //root?
#include "linearRegressionFit.cpp"
#endif

using namespace std;

DetectorConfiguration timePixChip {
	1, {0}, //nplanes, planeposition
	1, 256, 256 //pixelsize, xpixels, ypixels
};

void FitTracksTimePix(std::string inputfile) {

	//get tree from file
	TFile* file=openFile(inputfile);
	TTree* hitTable=getObjectFromFile<TTree>("Hits", file);

	//set adresses
	ULong64_t trigger=0;
	const std::vector<TimePixHit>* timepixHits=nullptr;
	hitTable->SetBranchAddress("timepix", &timepixHits );
//	hitTable->SetBranchAddress("trigger",&trigger);

	TimePixHoughTransformer houghTransform(-150E3, 150E3, timePixChip.planeymax(), 100,8);

	//loop over all entries
	const long long nEvents=hitTable->GetEntriesFast();
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%1))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		hitTable->GetEntry(iEvent);

		cout<<"got entry"<<endl;
		std::vector<PositionHit> spaceHit;
		double driftScale=1E-3;
		spaceHit=convertHits( *timepixHits, timePixChip.pixelsize, timePixChip.pixelsize, driftScale );

		cout<<"drawing cluster"<<endl;

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

			fits.push_back(fit);
		}

		HoughTransformer::drawClusters(spaceHit, timePixChip);
//		HoughTransformer::drawClusters(houghClusters, detector);
		for(auto& f : fits) f.draw(0, timePixChip.planePosition.back());
		if(std::cin.get()=='q') break;

		}
	}

}


