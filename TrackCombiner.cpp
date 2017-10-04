/*
 * TrackCombiner.cpp
 *
 *  Created on: Sep 21, 2017
 *      Author: cligtenb
 */

#include "TrackCombiner.h"

using namespace std;


TrackCombiner::TrackCombiner(std::string mimosaInput, std::string timepixInput, const DetectorConfiguration& telescope, const DetectorConfiguration& tpc) :
	telescope(telescope),
	tpcDetector(tpc),
	telescopeFitter(mimosaInput, telescope),
	tpcFitter(timepixInput,tpc)
{
	setTreeBranches();

	telescopeFitter.makeMask(5e3);
	telescopeFitter.setShifts( savedShifts ); //xy shift of planes
	telescopeFitter.setAngles( savedAngles ); //planes rotation along rz-axis in xy plane
	telescopeFitter.setCentres( savedCOM );
//	telescopeFitter.setSlopes( savedSlopes ); //slope along rz, undone later by rotation of all hits

	tpcFitter.makeMask(1e3);
//	tpcFitter.setSlopes( {0.006,-0.154} );
	tpcFitter.houghTransform.minCandidateSize=6;
	tpcFitter.houghTransform.minClusterSize=10;
}

TrackCombiner::~TrackCombiner() {
	if(outputFile) {
		outputFile->cd();
		fitResultTree.Write();
		triggerStatus.Write();
		frameStatusHistogram.Write();
//		outputFile->Close(); taken care of in custom uniqueptr destructor
	}
}

void TrackCombiner::setTreeBranches() {
	fitResultTree.Branch("telescopeFits", &telescopeFits);
	fitResultTree.Branch("timepixFits", &tpcFits);
	fitResultTree.Branch("ntimepixHits", &ntpcHits);
	fitResultTree.Branch("ntelescopeHits", &ntelescopeHits);
	fitResultTree.Branch("timepixClusterSize", &tpcClusterSize);
	fitResultTree.Branch("timepixHits", &tpcResiduals);
	fitResultTree.Branch("dxz", &dxz);
	fitResultTree.Branch("dyz", &dyz);
}

void TrackCombiner::openFile(std::string filename) {
	outputFile.reset(new TFile(filename.c_str(), "RECREATE"));
	fitResultTree.SetDirectory(outputFile.get());
}

bool isInsideDetector(const TVector3& p, const DetectorConfiguration& dc) {
	return not(
		p.x()<dc.xmin()
		or p.x() > dc.xmax()
		or p.y() < dc.ymin()
		or p.y() > dc.ymax() );
}

std::ostream& operator<<(std::ostream& os, const TVector3& v) {
	return cout<<"("<<v.x()<<", "<<v.y()<<", "<<v.z()<<")";
}

void TrackCombiner::processTracks() {
	nTelescopeTriggers=0;
	telescopeFitter.getEntry(0);
	for(int telescopeEntryNumber=0,tpcEntryNumber=0;
//			i<100000
			;) {

		// Get Entry and match trigger Numbers
		auto matchStatus=getAndMatchEntries(telescopeEntryNumber,tpcEntryNumber);
//		printTriggers(i,j);
//		if( cin.get()=='q') break;
		if( matchStatus == MatchResult::end ) break;
		else if( matchStatus == MatchResult::noMatch) continue;

		// Fit tracks

		//telescope
		telescopeFits.clear();
		auto telescopeHits=telescopeFitter.getSpaceHits();
		if( !telescopeFitter.passEvent(telescopeHits) ) { triggerStatus.Fill("Less than 4 planes hit in telescope", 1); frameStatus.replace(1, "less than 4 hits in telescope"); continue;} //minimal 4 planes!
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
		if(telescopeFits.empty()) { triggerStatus.Fill("All telescope clusters failed fit", 1);frameStatus.replace(2, "All telescope clusters failed fit"); continue; }

//		cout<<"telescope passed!"<<endl;

		//timepix
		tpcFits.clear();
		vector<const HoughTransformer::HitCluster*> tpcFittedClusters; //todo: replace with actual cluster for removing hits from cluster
		auto tpcHits=tpcFitter.getSpaceHits();
		if( !tpcFitter.passEvent(tpcHits) ) { triggerStatus.Fill("Less than 20 hits in tpc", 1);frameStatus.replace(3, "Less than 20 hits in tpc"); continue; }
		tpcHits=tpcFitter.rotateAndShift(tpcHits);
		tpcHits=tpcFitter.correctTimeWalk(tpcHits, 0.1209 /*mm/ns correction*/, 0.05 /*min ToT*/);
		auto tpcHistInTimePixFrame=tpcHits;//copy hits before rotation
		for(auto& h: tpcHits) {
			h.y=-h.y;
			h.RotatePosition(timepixYAngle, {11,0,6}, {0,1,0});
			h.RotatePosition(timepixXAngle, {0,-7,6}, {1,0,0});
//			h.RotatePosition(0.29, {0,-7,6}, {1,0,0});
			h.SetPosition(h.getPosition() + timepixShift);
		}
		const auto tpcClusters = tpcFitter.houghTransform(tpcHits);
		if(tpcClusters.size()>1) { auto mes="More than one cluster in tpc"; triggerStatus.Fill(mes, 1); frameStatus.replace(5, mes); continue; };
		for( auto& cluster : tpcClusters ) {
			if(cluster.size()<2) continue;
			auto fit= linearRegressionFit(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			tpcFits.push_back(fit);
			tpcFittedClusters.push_back(&cluster);
//			cout<<"tpc fit: "<<fit<<endl;
		}
		if(tpcFits.empty()) { triggerStatus.Fill("All tpc clusters failed fit", 1); frameStatus.replace(4, "all tpc clusters failed fit"); continue; }

//		cout<<"timepix passed!"<<endl;

		//match fits and clusters + calculate residuals
		tpcResiduals.clear();
		int nmatched=0;
		bool atLeastOneThroughDetector=false;
		vector<SimpleFitResult> telescopeTPCLines;
		std::vector<bool> tpcFitIsMatched(tpcFits.size()), telescopeFitIsMatched(telescopeFits.size());
		for(unsigned iFit=0; iFit<tpcFits.size();++iFit) {
			const auto& tpcFit=tpcFits[iFit];
			for(unsigned jFit=0; jFit<telescopeFits.size(); jFit++) {
				const auto& telescopeFit = telescopeFits[jFit];

				//check if going through detector
				TVector3 telescopePoint{ telescopeFit.xAt(timepixShift.z()), telescopeFit.yAt(timepixShift.z()), timepixShift.z() };
				telescopePoint-=timepixShift;
				telescopePoint=RotateAroundPoint(telescopePoint,-timepixXAngle, {0,-7,6}, {1,0,0});
				telescopePoint=RotateAroundPoint(telescopePoint,-timepixYAngle, {11,0,6}, {0,1,0});
				telescopePoint.SetY( -telescopePoint.y() );
				if( isInsideDetector( telescopePoint, timePixChip ) ) {
					atLeastOneThroughDetector=true;
				}

				//should we use the actual errors of the fit here?
				if( fabs( tpcFit.xAt(timepixShift.z())-telescopeFit.xAt(timepixShift.z()) ) < 1.5
					&& fabs( tpcFit.yAt(timepixShift.z())-telescopeFit.yAt(timepixShift.z()) ) < 1.5  ) {

					if(telescopeFitIsMatched[jFit]) { std::cerr<<"telescope fit is matched to 2 timepix clusters!"<<std::endl; }

					TVector3 average=tpcFittedClusters.at(iFit)->getAveragePosition();

					//draw line from telescope at z=0 to tpc centre
					SimpleFitResult telescopeTPCLine{
							(average.x()-telescopeFit.xAt(0))/average.z(),
							telescopeFit.xAt(0),
							(average.y()-telescopeFit.yAt(0))/average.z(),
							telescopeFit.yAt(0),
							0,0,0,0
					};
					telescopeTPCLines.push_back(telescopeTPCLine);

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

					dyz=(average.z()*telescopeFit.slope2+telescopeFit.intersept2-average.y())/sqrt(1+telescopeFit.slope2*telescopeFit.slope2);
					dxz=(average.z()*telescopeFit.slope1+telescopeFit.intersept1-average.x())/sqrt(1+telescopeFit.slope1*telescopeFit.slope1);

					++nmatched;
					telescopeFitIsMatched[jFit]=true;
					tpcFitIsMatched[iFit]=true;
					break;
				}
			}
		}


		//display event
		if( displayEvent ) {

			static TCanvas* timepixCanv=new TCanvas("timepix","Display of timepix event", 600,400);
			timepixCanv->cd();
			vector<SimpleFitResult> telescopeFitsInTimePixFrame;
			for(auto& f :telescopeFits ) //telescopeTPCLines)
				telescopeFitsInTimePixFrame.push_back(
						f.makeShifted(-timepixShift)
						 .makeRotated(-timepixXAngle, {0,-7,6}, {1,0,0})
						 .makeRotated(-timepixYAngle, {11,0,6}, {0,1,0})
						 .makeMirrorY() ); //
			tpcFitter.drawEvent(tpcHistInTimePixFrame, telescopeFitsInTimePixFrame);
			gPad->Update();

			if( telescopeFitter.processDrawSignals()  ) break;
//			displayEvent=false;
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
				auto message= atLeastOneThroughDetector ? "Telescope and tpc fits do not match" : "Telescope fit missed tpc";

				triggerStatus.Fill(message, 1);
				frameStatus.replace(10, message);
				if(displayEvent) cout<<"telescope and tpc fits do not match!"<<endl;
				continue;
			} else {
							}
		}

		//check if tpcEntry already has a matching fit
		if( find(tpcEntryHasMatchingFit.begin(), tpcEntryHasMatchingFit.end(), tpcEntryNumber )!=tpcEntryHasMatchingFit.end() ) {
			auto mes="Tpc entry already has a matching cluster"; triggerStatus.Fill(mes, 1); frameStatus.replace(10, mes); continue;
		}

//		cout<<"Success!";
		triggerStatus.Fill("Successful", 1);
		frameStatus.replace(20, "Successful");

		tpcEntryHasMatchingFit.push_back(tpcEntryNumber);
		//cleanup old entry info
		while(not tpcEntryHasMatchingFit.empty() and tpcEntryHasMatchingFit.front()<tpcEntryNumber-(telescopeFitter.triggerNumberEnd-telescopeFitter.triggerNumberBegin) ) tpcEntryHasMatchingFit.pop_front();

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

	}

	cout<<"highest telescope trigger (number of unique) "<< telescopeFitter.triggerNumberEnd<<" ("<<nTelescopeTriggers<<")"<<endl;
	cout<<"number of entries in timepix "<<tpcFitter.nEvents<<endl;

	cout<<"entries in tree "<<fitResultTree.GetEntriesFast()<<endl;

}

TrackCombiner::MatchResult TrackCombiner::getAndMatchEntries(
		int& telescopeEntry,
		int& tpcStartEntry) {

	//if previous frame did not increase and this frame did not increase, there is no related trigger to this frame
//	if(previous2TriggerNumberBegin==telescopeFitter.triggerNumberEnd) { cout<<"unrelated!"<<endl; }
//		triggerStatus.Fill("Unrelated telescope frame",1); return MatchResult::noMatch; }

	//if triggerNumberBegin decreased, we must get entries until the tpc triggernumber also decreases
	if(telescopeFitter.triggerNumberBegin<previousTriggerNumberBegin)
		do {
			if( !tpcFitter.getEntry(tpcStartEntry++) ) return MatchResult::end;
		} while ( (tpcFitter.triggerNumber+triggerOffset) % 32768 > previousTriggerNumberBegin );

	//get next entry until tpc trigger number is larger than or equal to begin
	do {
		if( !tpcFitter.getEntry(tpcStartEntry++) ) return MatchResult::end;
	} while( (tpcFitter.triggerNumber+triggerOffset) % 32768 < telescopeFitter.triggerNumberBegin );

	//if also larger than end: reached the end of this telescope frame, continue with next telescope frame;
	if( (tpcFitter.triggerNumber+triggerOffset) % 32768 > telescopeFitter.triggerNumberEnd) {
//		triggerStatus.Fill("Trigger numbers do not match", 1);

		if(frameStatus.priority) {
			frameStatusHistogram.Fill(frameStatus.message.c_str(), 1);
		}
		frameStatus.reset();
		nTelescopeTriggers+=telescopeFitter.triggerNumberEnd-previousTriggerNumberEnd;

		previousTriggerNumberBegin=telescopeFitter.triggerNumberBegin;
		previousTriggerNumberEnd=telescopeFitter.triggerNumberEnd;
		if( !telescopeFitter.getEntry(++telescopeEntry) ) return MatchResult::end;

//		cout<<"increased telescopeEntry, first time pix match was: "<<timepixEntryFirstMatch<<endl;

		tpcStartEntry=timepixEntryFirstMatch;
		hadFirstMatch=false;

		return MatchResult::noMatch;
	}

	if(not hadFirstMatch) {
		timepixEntryFirstMatch=tpcStartEntry-1;
		hadFirstMatch=true;
	}

	return MatchResult::match;

}

void TrackCombiner::printTriggers(int telescopeEntry, int tpcEntry) const {
	const int printEveryN=1;
	if( !(telescopeEntry%printEveryN) ) {
		cout<<"entry: "<<telescopeEntry<<"/"<<telescopeFitter.nEvents<<" ";
		cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd;
		cout<<" timepix triggerNumber: "<<tpcFitter.triggerNumber<<"="<<(tpcFitter.triggerNumber+triggerOffset) % 32768<<" in entry "<<tpcEntry<<endl;
	}
}
