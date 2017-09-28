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
#include "HitEntry.h"

#if 1 //root?
#include "linearRegressionFit.cpp"
#include "TrackFitter.cpp"
#include "TimePixFitter.cpp"
#include "makeNoisyPixelMask.cpp"
#include "ResidualHistogrammer.cpp"
#endif

#include "testBeamSetup.h"
#include "mimosaAlignment.h"
#include "relativeAlignment.h"

using namespace std;

//returns correlation factor
double CombineTracks(std::string mimosaInput, std::string timepixInput, int triggerOffset=0,  bool displayEvent=false) {

	//mimosa
	trackFitter telescopeFitter(mimosaInput, mimosa);

	telescopeFitter.makeMask(5e3);

	telescopeFitter.setShifts( savedShifts ); //xy shift of planes
	telescopeFitter.setAngles( savedAngles ); //planes rotation along rz-axis in xy plane
	telescopeFitter.setCentres( savedCOM );
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
	double dxz, dyz;

	std::unique_ptr<TFile> outputFile( displayEvent ? nullptr : new TFile("fitResults.root", "RECREATE") );
	TTree fitResultTree( "fitResults", "Tree with telescope and timepix fit results") ;
	fitResultTree.Branch("telescopeFits", &telescopeFits);
	fitResultTree.Branch("timepixFits", &tpcFits);
	fitResultTree.Branch("ntimepixHits", &ntpcHits);
	fitResultTree.Branch("ntelescopeHits", &ntelescopeHits);
	fitResultTree.Branch("timepixClusterSize", &tpcClusterSize);
	fitResultTree.Branch("timepixHits", &tpcResiduals);
	fitResultTree.Branch("dxz", &dxz);
	fitResultTree.Branch("dyz", &dyz);

	//histograms
//	auto residualHistograms=unique_ptr<ResidualHistogrammer>(displayEvent ? nullptr : new ResidualHistogrammer("timepixTelescopeTrackResiduals.root", timePixChip));
	TH1D timeDifference("timeDifference", "telescope time - timepix time;time [s];entries", 100, -1E-3,1E-3 );
	TH1D triggerStatus("triggerStatus", "status of trigger", 1, 0, 1);

	int nTelescopeTriggers=0,previousTriggerNumberBegin=0, previous2TriggerNumberBegin=0;
	for(int i=0,j=0;
//			i<100000
			;i++) {

		// Get Entry and match trigger Numbers
		previous2TriggerNumberBegin=previousTriggerNumberBegin;
		previousTriggerNumberBegin=telescopeFitter.triggerNumberBegin;
		if( !telescopeFitter.getEntry(i) ) break;

		if(i && !(i%10000) ) {
			cout<<"entry: "<<i<<"/"<<telescopeFitter.nEvents<<" ";
			cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd;
			cout<<" timepix triggerNumber: "<<tpcFitter.triggerNumber<<"="<<(tpcFitter.triggerNumber+triggerOffset) % 32768<<" in entry "<<j<<endl;
//			if(cin.get()=='q') break;
		}

		//if previous frame did not increase and this frame did not increase, there is no related trigger to this frame


		//if triggerNumberBegin decreased, we must get entries until the tpc triggernumber also decreases
		if(telescopeFitter.triggerNumberBegin<previousTriggerNumberBegin)
			while( tpcFitter.getEntry(j++) && (tpcFitter.triggerNumber+triggerOffset) % 32768 > previousTriggerNumberBegin ) {};

		//get next entry until tpc trigger number is larger than or equal to begin
		while( (tpcFitter.triggerNumber+triggerOffset) % 32768 <telescopeFitter.triggerNumberBegin && tpcFitter.getEntry(j++) ) {};

		//if also larger than end, continue;
		if( (tpcFitter.triggerNumber+triggerOffset) % 32768 > telescopeFitter.triggerNumberEnd) { triggerStatus.Fill("Trigger numbers do not match", 1); continue;}



		if(previousTriggerNumberBegin!=telescopeFitter.triggerNumberBegin) {
			nTelescopeTriggers++;
		}

		// Fit tracks

		//telescope
		telescopeFits.clear();
		auto telescopeHits=telescopeFitter.getSpaceHits();
		if( !telescopeFitter.passEvent(telescopeHits) ) { triggerStatus.Fill("Less than 4 planes hit in telescope", 1); continue;} //minimal 4 planes!
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
		if(telescopeFits.empty()) { triggerStatus.Fill("All telescope clusters failed fit", 1); continue; }

//		cout<<"telescope passed!"<<endl;

		//timepix
		tpcFits.clear();
		vector<HoughTransformer::HitCluster*> tpcFittedClusters;
		auto tpcHits=tpcFitter.getSpaceHits();
		if( !tpcFitter.passEvent(tpcHits) ) { triggerStatus.Fill("Less than 20 hits in tpc", 1); continue; }
		tpcHits=tpcFitter.rotateAndShift(tpcHits);
		tpcHits=tpcFitter.correctTimeWalk(tpcHits, 0.1209 /*mm/ns correction*/, 0.05 /*min ToT*/);
		auto tpcHistInTimePixFrame=tpcHits;//save hits before rotation
		for(auto& h: tpcHits) {
			h.y=-h.y;
			h.RotatePosition(timepixYAngle, {11,0,6}, {0,1,0});
			h.RotatePosition(timepixXAngle, {0,-7,6}, {1,0,0});
//			h.RotatePosition(0.29, {0,-7,6}, {1,0,0});
			h.SetPosition(h.getPosition() + timepixShift);
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
		if(tpcFits.empty()) { triggerStatus.Fill("All tpc clusters failed fit", 1); continue; }

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
				if( fabs( tpcFit.xAt(timepixShift.z())-telescopeFit.xAt(timepixShift.z()) ) < 1.5
				    && fabs( tpcFit.yAt(timepixShift.z())-telescopeFit.yAt(timepixShift.z()) ) < 1.5  ) {

//					if(telescopeFitIsMatched[jFit]) { std::cerr<<"telescope fit is matched to 2 timepix clusters!"<<std::endl; }

					auto residuals=calculateResiduals(*tpcFittedClusters.at(iFit), telescopeFit);
					//rotate back to frame of timepix
					for(auto& r:residuals) {
						auto v=r.getVector();
						v.Rotate(-timepixXAngle, {1,0,0});
						v.Rotate(-timepixYAngle, {0,1,0});
						r.setVector(v);
					}
//					tpcResiduals.insert(tpcResiduals.end(), residuals.begin(), residuals.end() );
					tpcResiduals.emplace_back( residuals.begin(), residuals.end() );//construct new vector with entries just for this cluster in tpcResiduals

					TVector3 average=tpcFittedClusters.at(iFit)->getAveragePosition();
					dyz=(average.z()*telescopeFit.slope2+telescopeFit.intersept2-average.y())/sqrt(1+telescopeFit.slope2*telescopeFit.slope2);
					dxz=(average.z()*telescopeFit.slope1+telescopeFit.intersept1-average.x())/sqrt(1+telescopeFit.slope1*telescopeFit.slope1);

					++nmatched;
					telescopeFitIsMatched[jFit]=true;
					tpcFitIsMatched[iFit]=true;
					break;
				}
			}
		}
//		cout<<"matched "<<nmatched<<" clusters"<<endl;
		//remove unmatched clusters
		bool removeUnmatched=true;
		if(removeUnmatched) {
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

			if( telescopeFits.empty() || tpcFits.empty() ) {
				triggerStatus.Fill("Telescope and tpc fits do not match", 1);
				continue;
			} else {
				triggerStatus.Fill("Successful", 1);
			}
		}

		//display event
		if( displayEvent ) {

			static TCanvas* timepixCanv=new TCanvas("timepix","Display of timepix event", 600,400);
			timepixCanv->cd();
			vector<SimpleFitResult> telescopeFitsInTimePixFrame;
			for(auto& f : telescopeFits)
				telescopeFitsInTimePixFrame.push_back( f.makeShifted(-timepixShift).makeRotated(-timepixXAngle, {0,-7,6}, {1,0,0}).makeMirrorY() ); //
			tpcFitter.drawEvent(tpcHistInTimePixFrame, telescopeFitsInTimePixFrame);
//			for (auto& f : telescopeFits)
//				f.draw(timePixChip.zmin()+timepixShift.z(), timePixChip.zmax()+timepixShift.z());
//			for (auto& f : tpcFits)
//				f.draw( timePixChip.zmin()+timepixShift.z(), timePixChip.zmin()+timepixShift.z() );
			gPad->Update();

//			if( telescopeFitter.processDrawSignals()  ) break;
//			static TCanvas* mimosaCanv=new TCanvas("mimosa","Display of mimosa event", 600,400);
//			mimosaCanv->cd();
//			telescopeFitter.drawEvent( telescopeHits, telescopeFits );

			/*
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
//		displayEvent=false;


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
		triggerStatus.Write();
		outputFile->Close();
	}

//	return gr.GetCorrelationFactor();
	return 0;
}


