/*
 * CombineTracks.cpp
 *
 *  Created on: Jul 19, 2017
 *      Author: cligtenb
 */


#include <string>

#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"

#include "TrackFitter.h"
#include "TimePixFitter.h"
#include "linearRegressionFit.h"
#include "ResidualHistogrammer.h"

#if 1 //root?
#include "linearRegressionFit.cpp"
#include "TrackFitter.cpp"
#include "TimePixFitter.cpp"
#include "makeNoisyPixelMask.cpp"
#include "ResidualHistogrammer.cpp"
#endif

#include "mimosaAlignment.h"

using namespace std;

struct TimePixDetectorConfiguration : DetectorConfiguration {
	static constexpr double driftSpeed=0.075; //mm/ns
	TimePixDetectorConfiguration() : DetectorConfiguration{
		1, {0}, //nplanes, planeposition
		0.055, 256, 256 //pixelsize, xpixels, ypixels
	} {};
	virtual double xmin() const {return driftSpeed; }
	virtual double xmax() const {return 400*driftSpeed; }
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

struct ResidualTreeEntry {
	ResidualTreeEntry() : x(), y(), z(), ToT(), row(), col() {};
	ResidualTreeEntry(Residual r) :
		x(r.x),
		y(r.y),
		z(r.z),
		ToT(r.h.ToT),
		row(r.h.row),
		col(r.h.column)
	{};
	double x, y, z;
	int ToT;
	int row, col;
};
#pragma link C++ class std::vector<ResidualTreeEntry>+;

//returns correlation factor
double CombineTracks(std::string mimosaInput, std::string timepixInput, int triggerOffset=0,  bool displayEvent=false) {

	//mimosa
	trackFitter telescopeFitter(mimosaInput, mimosa);

	telescopeFitter.makeMask(5e3);

	telescopeFitter.setShifts( savedShifts ); //xy shift of planes
	telescopeFitter.setAngles( savedAngles ); //planes rotation along z-axis in xy plane
//	telescopeFitter.setSlopes( savedSlopes ); //slope along z, undone later by rotation of all hits


	//timepix
	TimePixFitter tpcFitter(timepixInput,timePixChip);

	tpcFitter.makeMask(1e3);

//	tpcFitter.setSlopes( {0.006,-0.154} );
	tpcFitter.houghTransform.minCandidateSize=6;
	tpcFitter.houghTransform.minClusterSize=10;

	vector<SimpleFitResult> telescopeFits;
	vector<SimpleFitResult> tpcFits;
	vector<ResidualTreeEntry> tpcResiduals;
	int ntpcHits, ntelescopeHits;

	std::unique_ptr<TFile> outputFile( displayEvent ? nullptr : new TFile("fitResults.root", "RECREATE") );
	TTree fitResultTree( "fitResultTree", "Tree with telescope and timepix fit results") ;
	fitResultTree.Branch("telescopeFits", &telescopeFits);
	fitResultTree.Branch("timepixFits", &tpcFits);
	fitResultTree.Branch("ntimepixHits", &ntpcHits);
	fitResultTree.Branch("ntelescopeHits", &ntelescopeHits);
	fitResultTree.Branch("timepixResiduals", &tpcResiduals);

	//histograms
	auto residualHistograms=unique_ptr<ResidualHistogrammer>(displayEvent ? nullptr : new ResidualHistogrammer("timepixTelescopeTrackResiduals.root", timePixChip));
	TH1D timeDifference("timeDifference", "telescope time - timepix time;time [s];entries", 100, -1E-3,1E-3 );

	int nTelescopeTriggers=0,previousTriggerNumberBegin=0;
	for(int i=0,j=0;
//			i<20
			;i++) {

		// Get Entry and match trigger Numbers
		previousTriggerNumberBegin=telescopeFitter.triggerNumberBegin;
		if( !telescopeFitter.getEntry(i) ) break;

		//if triggerNumberBegin decreased, we must get entries until the tpc triggernumber also decreases
		if(telescopeFitter.triggerNumberBegin<previousTriggerNumberBegin)
			while( (tpcFitter.triggerNumber+triggerOffset % 32768) > previousTriggerNumberBegin && tpcFitter.getEntry(j++) ) {};

		//get next entry until tpc trigger number is larger than begin
		while( (tpcFitter.triggerNumber+triggerOffset % 32768) <telescopeFitter.triggerNumberBegin && tpcFitter.getEntry(j++) ) {};

		//if also larger than end, continue;
		if( (tpcFitter.triggerNumber+triggerOffset % 32768) > telescopeFitter.triggerNumberEnd) continue;

//		cout<<"entry: "<<i<<"/"<<telescopeFitter.nEvents<<" ";
//		cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd;
//		cout<<" timepix triggerNumber: "<<tpcFitter.triggerNumber<<"="<<(tpcFitter.triggerNumber % 32768)<<" in entry "<<j<<endl;

		if(previousTriggerNumberBegin!=telescopeFitter.triggerNumberBegin) {
			nTelescopeTriggers++;
		}

		// Fit tracks

		//telescope
		telescopeFits.clear();
		auto telescopeHits=telescopeFitter.getSpaceHits();
		if( !telescopeFitter.passEvent(telescopeHits) ) continue;
		telescopeHits=telescopeFitter.rotateAndShift(telescopeHits);
		for(auto&v:telescopeHits) for(auto&h:v) {
			h.RotatePosition(-savedSlopes.first, {mimosa.getCentre().first,mimosa.getCentre().second,0}, {0,1,0} );
			h.RotatePosition(savedSlopes.second, {mimosa.getCentre().first,mimosa.getCentre().second,0}, {1,0,0} );
		}
		auto telescopeClusters = telescopeFitter.houghTransform(telescopeHits);
		for( auto& cluster : telescopeClusters) {
			if(cluster.size()<5 || cluster.getNPlanesHit()<=4) continue;
			auto fit= linearRegressionFit(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			telescopeFits.push_back(fit);
		}
		if(telescopeFits.empty()) continue;

//		cout<<"telescope passed!"<<endl;

		//timepix
		tpcFits.clear();
		auto tpcHits=tpcFitter.getSpaceHits();
		if( !tpcFitter.passEvent(tpcHits) ) continue;
		tpcHits=tpcFitter.rotateAndShift(tpcHits);
		for(auto& h: tpcHits) {
			h.y=-h.y;
//			h.RotatePosition(-0.0, {0,0,10}, {0,1,0});
			h.RotatePosition(0.29, {0,-7,6}, {1,0,0});
			h.SetPosition(h.getPosition() + TVector3(6,14,-380) );
		}
		auto tpcClusters = tpcFitter.houghTransform(tpcHits);
		for( auto& cluster : tpcClusters ) {
			if(cluster.size()<2) continue;
			auto fit= linearRegressionFit(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			tpcFits.push_back(fit);
//			cout<<"tpc fit: "<<fit<<endl;
		}
		if(tpcFits.empty()) continue;

//		cout<<"timepix passed!"<<endl;

		//display event
		if( displayEvent ) {

			static TCanvas* timepixCanv=new TCanvas("timepix","Display of timepix event", 600,400);
			timepixCanv->cd();
			tpcFitter.drawEvent(tpcHits, tpcFits);
			for (auto& f : telescopeFits)
				f.draw(timePixChip.zmin()-380, timePixChip.zmax()-380);
			for (auto& f : tpcFits)
				f.draw( timePixChip.zmin()-380, timePixChip.zmin()-380 );
			gPad->Update();

//			if( telescopeFitter.processDrawSignals()  ) break;
/*
			static TCanvas* mimosaCanv=new TCanvas("mimosa","Display of mimosa event", 600,400);
			mimosaCanv->cd();
			telescopeFitter.drawEvent( telescopeHits, telescopeFits );

			DetectorConfiguration combinedSetup{
					2, {-350,200 }, //nplanes, planeposition
					0.001, int(1000*mimosa.xmax()), int(1000*mimosa.ymax()) //pixelsize, xpixels, ypixels
				};



			std::vector<PositionHit> combinedHits;
			for(auto& v: telescopeHits)
				for(auto& h : v)
					combinedHits.push_back(h);
			for(auto& h: tpcHits) combinedHits.push_back(h);


			HoughTransformer::drawCluster(combinedHits, combinedSetup);
			for (auto& f : tpcFits)
				f.draw( combinedSetup.zmin(), combinedSetup.zmax() );
			for (auto& f : telescopeFits)
				f.draw( combinedSetup.zmin(), combinedSetup.zmax() );
//*/
			gPad->Update();

			if( telescopeFitter.processDrawSignals()  ) break;
		}

		auto residuals=calculateResiduals(tpcClusters.front(), telescopeFits[0]);
		//rotate back to frame of timepix
		for(auto& r:residuals) {
			auto v=r.getVector();
			v.Rotate(-0.29, {1,0,0});
			r.setVector(v);
		}
		tpcResiduals.clear();
		tpcResiduals.insert(tpcResiduals.begin(), residuals.begin(), residuals.end());
		if(residualHistograms) {
			if( averageResidual(residuals).Perp()<3 )
				residualHistograms->fill(residuals);
		}

		ntpcHits=tpcHits.size();
		ntelescopeHits=0;
		for(auto& v : telescopeHits ) ntelescopeHits+=v.size();
		fitResultTree.Fill();

		double telescopeSeconds=telescopeFitter.timestamp/40.E6, tpcSeconds=tpcFitter.timestamp/4096.*25E-9;
		static double firstTimeDifference=telescopeSeconds-tpcSeconds;
		timeDifference.Fill( telescopeSeconds-tpcSeconds-firstTimeDifference );

	}

	cout<<"highest telescope trigger (number of unique) "<< telescopeFitter.triggerNumberEnd<<" ("<<nTelescopeTriggers<<")"<<endl;
	cout<<"number of entries in timepix "<<tpcFitter.nEvents<<endl;

	cout<<"entries in tree "<<fitResultTree.GetEntriesFast()<<endl;

	fitResultTree.Draw("timepixFits[0].intersept1:telescopeFits[0].intersept1");
//	fitResultTree.Draw("timepixFits[0].slope1:telescopeFits[].slope2");//, "fabs(timepixFits.slope1)<1");
//	fitResultTree.Draw("ntelescopeHits:ntimepixHits","ntimepixHits<1000", "prof");//, "fabs(timepixFits.slope1)<1");
//	timeDifference.DrawClone();
//	residualHistograms->DrawClone();

	TGraph gr(fitResultTree.GetSelectedRows(),
    		fitResultTree.GetV2(), fitResultTree.GetV1());

	if(outputFile) {
		outputFile->cd();
		fitResultTree.Write();
		timeDifference.Write();
		outputFile->Close();
	}

	return gr.GetCorrelationFactor();
}


