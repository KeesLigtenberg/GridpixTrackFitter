/*
 * CombineTracks.cpp
 *
 *  Created on: Jul 19, 2017
 *      Author: cligtenb
 */


#include <string>

#include "TTree.h"

#include "TrackFitter.h"
#include "TimePixFitter.h"
#include "linearRegressionFit.h"

#if 1 //root?
#include "linearRegressionFit.cpp"
#include "TrackFitter.cpp"
#include "TimePixFitter.cpp"
#include "makeNoisyPixelMask.cpp"
#endif

#include "mimosaAlignment.h"

using namespace std;

struct TimePixDetectorConfiguration : DetectorConfiguration {
	static constexpr double driftSpeed=0.075; //mm/ns
	TimePixDetectorConfiguration() : DetectorConfiguration{
		1, {0}, //nplanes, planeposition
		0.055, 256, 256 //pixelsize, xpixels, ypixels
	} {};
	virtual double xmin() const {return -200*driftSpeed; }
	virtual double xmax() const {return 200*driftSpeed; }
	virtual double zmin() const {return 0; }
	virtual double zmax() const {return 256*pixelsize; };
} timePixChip;

const DetectorConfiguration mimosa= {
	6, //planes
//	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	{0, 15.8, 31.8, 143.1, 161.55, 179.91 }, //plane positions as measured
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other!
};


void CombineTracks(std::string mimosaInput, std::string timepixInput, int offset=0) {

	//mimosa
	trackFitter telescopeFitter(mimosaInput, mimosa);

	telescopeFitter.makeMask(5e3);

	telescopeFitter.setShifts( savedShifts );
	telescopeFitter.setAngles( savedAngles );
	telescopeFitter.setSlopes( savedSlopes );


	//timepix
	TimePixFitter tpcFitter(timepixInput,timePixChip);

	tpcFitter.makeMask(1e3);

	tpcFitter.setSlopes( {0.006,-0.154} );
	tpcFitter.houghTransform.minCandidateSize=6;
	tpcFitter.houghTransform.minClusterSize=10;


	vector<SimpleFitResult> telescopeFits;
	vector<SimpleFitResult> tpcFits;

	TFile outputFile("fitResults.root", "RECREATE");
	TTree fitResultTree( "fitResultTree", "Tree with telescope and timepix fit results") ;
	fitResultTree.Branch("telescopeFits", &telescopeFits);
	fitResultTree.Branch("timepixFits", &tpcFits);

	int nTelescopeTriggers=0,previousTriggerNumberBegin=0;
	for(int i=0;
//			i<1000 &&
			telescopeFitter.getEntry(i) &&
					tpcFitter.getEntry(telescopeFitter.triggerNumberBegin+offset);
			i++) {

		cout<<"entry "<<i<<"/"<<telescopeFitter.nEvents<<" ";
		cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd<<endl;

		if(previousTriggerNumberBegin!=telescopeFitter.triggerNumberBegin) {
			nTelescopeTriggers++;
			previousTriggerNumberBegin=telescopeFitter.triggerNumberBegin;
		}

//		continue;

		auto telescopeHits=telescopeFitter.getSpaceHits();
		if( !telescopeFitter.passEvent(telescopeHits) ) continue;
		telescopeHits=telescopeFitter.rotateAndShift(telescopeHits);
		auto telescopeClusters = telescopeFitter.houghTransform(telescopeHits);
		for( auto& cluster : telescopeClusters) {
			if(cluster.size()<2 || cluster.getNPlanesHit()<=1) continue;
			auto fit= linearRegressionFit(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			telescopeFits.push_back(fit);
		}

		auto tpcHits=tpcFitter.getSpaceHits();
		if( !tpcFitter.passEvent(tpcHits) ) continue;
		tpcHits=tpcFitter.rotateAndShift(tpcHits);
		auto tpcClusters = tpcFitter.houghTransform(tpcHits);
		for( auto& cluster : tpcClusters) {
			if(cluster.size()<2) continue;
			auto fit= linearRegressionFit(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			tpcFits.push_back(fit);
		}

		//display event
		bool displayEvent=true;
		if( displayEvent ) {
			static TCanvas* timepixCanv=new TCanvas("timepix","Display of timepix event", 600,400);
			timepixCanv->cd();
			tpcFitter.drawEvent(tpcHits, tpcFits);

			static TCanvas* mimosaCanv=new TCanvas("mimosa","Display of mimosa event", 600,400);
			mimosaCanv->cd();
//			telescopeFitter.drawEvent( telescopeHits, telescopeFits );

			DetectorConfiguration combinedSetup{
					2, {0,200 }, //nplanes, planeposition
					0.001, int(1000*mimosa.xmax()), int(1000*mimosa.ymax()) //pixelsize, xpixels, ypixels
				};

			std::vector<std::vector<PositionHit> > combinedHits(1);
			for(auto& v: telescopeHits)
				for(auto& h : v)
					combinedHits.at(0).push_back(h);
			rotateHits()
			combinedHits.push_back(tpcHits);
			HoughTransformer::drawClusters(combinedHits, combinedSetup);
			for (auto& f : tpcFits)
				f.draw( combinedSetup.zmin(), combinedSetup.zmax() );
			for (auto& f : telescopeFits)
				f.draw( combinedSetup.zmin(), combinedSetup.zmax() );

			gPad->Update();

			if( telescopeFitter.processDrawSignals()  ) break;
		}

		fitResultTree.Fill();

		telescopeFits.clear();
		tpcFits.clear();
	}

	cout<<"number of telescope unique telescope triggers "<<nTelescopeTriggers<<endl;
	cout<<"number of entries in timepix "<<tpcFitter.nEvents<<endl;

	cout<<"entries in tree "<<fitResultTree.GetEntriesFast()<<endl;

	fitResultTree.DrawClone("timepixFits.slope1:telescopeFits.slope2");

	fitResultTree.Write();
	outputFile.Close();

}


