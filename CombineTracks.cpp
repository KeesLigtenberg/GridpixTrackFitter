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
const double timepixZCenter=-374;

const DetectorConfiguration mimosa= {
	6, //planes
//	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	{0, 15.8, 31.8, 143.1, 161.55, 179.91 }, //plane positions as measured
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other!
};

struct HitEntry {
	HitEntry() : rx(), ry(), rz(),x(),y(),z(), ToT(), row(), col() {};
	HitEntry(Residual r) :
		rx(r.x),
		ry(r.y),
		rz(r.z),
		x(r.h.x),
		y(r.h.y),
		z(r.h.z),
		ToT(r.h.ToT),
		row(r.h.row),
		col(r.h.column)
	{};
	double rx, ry, rz;
	double x, y, z;
	int ToT;
	int row, col;
};
#pragma link C++ class std::vector<HitEntry>+;
#pragma link C++ class std::vector< std::vector<HitEntry> >+;

//returns correlation factor
double CombineTracks(std::string mimosaInput, std::string timepixInput, int triggerOffset=0,  bool displayEvent=false) {

	//mimosa
	trackFitter telescopeFitter(mimosaInput, mimosa);

	telescopeFitter.makeMask(5e3);

	telescopeFitter.setShifts( savedShifts ); //xy shift of planes
	telescopeFitter.setAngles( savedAngles ); //planes rotation along rz-axis in xy plane
//	telescopeFitter.setSlopes( savedSlopes ); //slope along rz, undone later by rotation of all hits


	//timepix
	TimePixFitter tpcFitter(timepixInput,timePixChip);

	tpcFitter.makeMask(1e3);

//	tpcFitter.setSlopes( {0.006,-0.154} );
	tpcFitter.houghTransform.minCandidateSize=6;
	tpcFitter.houghTransform.minClusterSize=10;

	vector<SimpleFitResult> telescopeFits;
	vector<SimpleFitResult> tpcFits;
	vector< vector<HitEntry> > tpcResiduals;
	int ntpcHits, ntelescopeHits;
	vector<int> tpcClusterSize;

	std::unique_ptr<TFile> outputFile( displayEvent ? nullptr : new TFile("fitResults.root", "RECREATE") );
	TTree fitResultTree( "fitResults", "Tree with telescope and timepix fit results") ;
	fitResultTree.Branch("telescopeFits", &telescopeFits);
	fitResultTree.Branch("timepixFits", &tpcFits);
	fitResultTree.Branch("ntimepixHits", &ntpcHits);
	fitResultTree.Branch("ntelescopeHits", &ntelescopeHits);
	fitResultTree.Branch("timepixClusterSize", &tpcClusterSize);
	fitResultTree.Branch("timepixHits", &tpcResiduals);

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
		vector<HoughTransformer::HitCluster*> tpcFittedClusters;
		auto tpcHits=tpcFitter.getSpaceHits();
		if( !tpcFitter.passEvent(tpcHits) ) continue;
		tpcHits=tpcFitter.rotateAndShift(tpcHits);
		tpcHits=tpcFitter.correctTimeWalk(tpcHits, 0.1209 /*mm/ns correction*/, 0.05 /*min ToT*/);
		for(auto& h: tpcHits) {
			h.y=-h.y;
//			h.RotatePosition(-0.0, {0,0,10}, {0,1,0});
			h.RotatePosition(0.29, {0,-7,6}, {1,0,0});
			h.SetPosition(h.getPosition() + TVector3(6+0.62,14-0.0543605,timepixZCenter) );
		}
		auto tpcClusters = tpcFitter.houghTransform(tpcHits);
		for( auto& cluster : tpcClusters ) {
			if(cluster.size()<2) continue;
			auto fit= linearRegressionFit(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			tpcFits.push_back(fit);
			tpcFittedClusters.push_back(&cluster);
//			cout<<"tpc fit: "<<fit<<endl;
		}
		if(tpcFits.empty()) continue;

//		cout<<"timepix passed!"<<endl;

		//match fits and clusters
		tpcResiduals.clear();
		int nmatched=0;
		std::vector<bool> tpcFitIsMatched(tpcFits.size()), telescopeFitIsMatched(telescopeFits.size());
		for(unsigned iFit=0; iFit<tpcFits.size();++iFit) {
			const auto& tpcFit=tpcFits[iFit];
			for(unsigned jFit=0; jFit<telescopeFits.size(); jFit++) {
				const auto& telescopeFit = telescopeFits[jFit];
				//should we use the actual errors of the fit here?
				if( fabs( tpcFit.xAt(timepixZCenter)-telescopeFit.xAt(timepixZCenter) ) < 1.5
				    && fabs( tpcFit.yAt(timepixZCenter)-telescopeFit.yAt(timepixZCenter) ) < 1.5  ) {

					if(telescopeFitIsMatched[jFit]) { std::cerr<<"telescope fit is matched to 2 timepix clusters!"<<std::endl;}

					auto residuals=calculateResiduals(*tpcFittedClusters.at(iFit), telescopeFit);
					//rotate back to frame of timepix
					for(auto& r:residuals) {
						auto v=r.getVector();
						v.Rotate(-0.29, {1,0,0});
						r.setVector(v);
					}
//					tpcResiduals.insert(tpcResiduals.end(), residuals.begin(), residuals.end() );
					tpcResiduals.emplace_back( residuals.begin(), residuals.end() );//construct new vector with entries just for this cluster in tpcResiduals

					if(residualHistograms) {
						residualHistograms->fill(residuals);
					}
					++nmatched;
					telescopeFitIsMatched[jFit]=true;
					tpcFitIsMatched[iFit]=true;
					break;
				}
			}
		}
		cout<<"matched "<<nmatched<<" clusters"<<endl;
		//remove unmatched clusters
		int iFit=0;
		for(auto it=tpcFits.begin(); it!=tpcFits.end();) {
			if(!tpcFitIsMatched[iFit++]) {
				it=tpcFits.erase(it);
			} else {
				++it;
			}
		}
		int jFit=0;
		for(auto it=telescopeFits.begin(); it!=telescopeFits.end();) {
			if(!telescopeFitIsMatched[jFit++]) {
				it=telescopeFits.erase(it);
			} else {
				++it;
			}
		}

		//display event
		if( displayEvent ) {

			static TCanvas* timepixCanv=new TCanvas("timepix","Display of timepix event", 600,400);
			timepixCanv->cd();
			tpcFitter.drawEvent(tpcHits, tpcFits);
			for (auto& f : telescopeFits)
				f.draw(timePixChip.zmin()+timepixZCenter, timePixChip.zmax()+timepixZCenter);
			for (auto& f : tpcFits)
				f.draw( timePixChip.zmin()+timepixZCenter, timePixChip.zmin()+timepixZCenter );
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

		tpcClusterSize.clear();
		for(unsigned iClust=0; iClust<tpcFittedClusters.size(); ++iClust) {
			if(tpcFitIsMatched[iClust])	tpcClusterSize.push_back( tpcFittedClusters[iClust]->size());//count cluster size of only fitted clusters
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

//	fitResultTree.Draw("timepixFits[0].slope1:telescopeFits[].slope2");//, "fabs(timepixFits.slope1)<1");
//	fitResultTree.Draw("ntelescopeHits:ntimepixHits","ntimepixHits<1000", "prof");//, "fabs(timepixFits.slope1)<1");
//	timeDifference.DrawClone();
//	residualHistograms->DrawClone();

//	fitResultTree.Draw("timepixFits[0].intersept2-374*timepixFits[0].slope2:telescopeFits[0].intersept2-374*telescopeFits[0].slope2");
//	TGraph gr(fitResultTree.GetSelectedRows(),
//    		fitResultTree.GetV2(), fitResultTree.GetV1());
//	TFitResultPtr fit=gr.Fit("x+[0]", "S");
//	cout<<"x shift is "<<fit->Parameter(0);

	if(outputFile) {
		outputFile->cd();
		fitResultTree.Write();
		timeDifference.Write();
		outputFile->Close();
	}

//	return gr.GetCorrelationFactor();
	return 0;
}


