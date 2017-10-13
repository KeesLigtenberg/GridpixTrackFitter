/*
 * TrackCombiner.cpp
 *
 *  Created on: Sep 21, 2017
 *      Author: cligtenb
 */

#include "TrackCombiner.h"

using namespace std;


//note the template causes multiple canvases, because each one has its own static variable
template <class T> //std::vector<PositionHit>
void TrackCombiner::drawEvent(const T& hits,
		const std::vector<FitResult3D>& fits) {
	static TCanvas* timepixCanv=new TCanvas(typeid(T).name(),"Display of timepix event", 600,400);
	timepixCanv->cd();
	auto& ra=alignment.relativeAlignment;
	vector<FitResult3D> fitsInTimePixFrame;
	for(auto& f :tpcFits ) //telescopeTPCLines)
		fitsInTimePixFrame.push_back(
				f.makeShifted(-ra.shift)
				 .makeRotated(-ra.angle[2], ra.getCOM(), {0,0,1})
				 .makeRotated(-ra.angle[0], ra.getCOM(), {1,0,0})
				 .makeRotated(-ra.angle[1], ra.getCOM(), {0,1,0})
				 .makeMirrorY() ); //
	auto hitsInTimepixFrame=hits;
	for(auto& h: hitsInTimepixFrame) {
		h.SetPosition(h.getPosition() - ra.shift);
		h.RotatePosition(-ra.angle[2], ra.getCOM(), {0,0,1});
		h.RotatePosition(-ra.angle[0], ra.getCOM(), {1,0,0});
		h.RotatePosition(-ra.angle[1], ra.getCOM(), {0,1,0});
		h.y=-h.y;
	}
	HoughTransformer::drawCluster(hitsInTimepixFrame, timePixChip);
	for (auto& f : fitsInTimePixFrame)
		f.draw(timePixChip.zmin(), timePixChip.zmax());
	gPad->Update();
}



TrackCombiner::TrackCombiner(std::string mimosaInput, std::string timepixInput, const DetectorConfiguration& telescope, const DetectorConfiguration& tpc) :
	telescope(telescope),
	tpcDetector(tpc),
	telescopeFitter(mimosaInput, telescope),
	tpcFitter(timepixInput,tpc)
{
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
		treeBuffer.Write();
		triggerStatusHistogram.Write();
		frameStatusHistogram.Write();
		timepixStatusHistogram.Write();
	}
}


void BufferedTreeFiller::Fill(const TreeEntry& entry) {
	currentEntry=entry;
	fitResultTree.Fill();
}

void BufferedTreeFiller::Write() {
	buffer.writeBuffer( [this](TreeEntry&t){this->Fill(t);} );
	fitResultTree.Write();
}

void BufferedTreeFiller::setTreeBranches() {
	fitResultTree.Branch("telescopeFits", &currentEntry.telescopeFits);
	fitResultTree.Branch("timepixFits", &currentEntry.tpcFits);
	fitResultTree.Branch("ntimepixHits", &currentEntry.ntpcHits);
	fitResultTree.Branch("ntelescopeHits", &currentEntry.ntelescopeHits);
	fitResultTree.Branch("timepixClusterSize", &currentEntry.tpcClusterSize);
	fitResultTree.Branch("timepixHits", &currentEntry.tpcResiduals);
	fitResultTree.Branch("dxz", &currentEntry.dxz);
	fitResultTree.Branch("dyz", &currentEntry.dyz);
	fitResultTree.Branch("nresiduals", &currentEntry.nresiduals);
	fitResultTree.Branch("nfitted", &currentEntry.nfitted);
}

void TrackCombiner::openFile(std::string filename) {
	outputFile.reset(new TFile(filename.c_str(), "RECREATE"));
	treeBuffer.SetTreeDirectory(outputFile.get());
}

bool isInsideDetector(const TVector3& p, const DetectorConfiguration& dc, double threshold=0) {
	return not( p.x()<dc.xmin()-threshold
		or p.x() > dc.xmax()+threshold
		or p.y() < dc.ymin()-threshold
		or p.y() > dc.ymax()+threshold );
}

//returns 0 if not through
//returns 1 if through either front or back
//returns 2 if through both front and back
int goesThroughTimepix(const FitResult3D& fit, const Alignment& alignment, const DetectorConfiguration& dc=timePixChip, double threshold=0) {
	int nPlanes=0;
	auto& ra=alignment.relativeAlignment;
	for(double z : {dc.zmin(), dc.zmax() }) {
		TVector3 telescopePoint{ fit.xAt(ra.shift.z()+z), fit.yAt(ra.shift.z()+z), ra.shift.z()+z };
		telescopePoint-=ra.shift;
		telescopePoint=RotateAroundPoint(telescopePoint,-ra.angle[2], ra.getCOM(), {0,0,1});
		telescopePoint=RotateAroundPoint(telescopePoint,-ra.angle[0], ra.getCOM(), {1,0,0});
		telescopePoint=RotateAroundPoint(telescopePoint,-ra.angle[1], ra.getCOM(), {0,1,0});
		telescopePoint.SetY( -telescopePoint.y() );
		if( isInsideDetector( telescopePoint, timePixChip ) ) {
			++nPlanes;
		}
	}
	return nPlanes;
}

std::ostream& operator<<(std::ostream& os, const TVector3& v) {
	return cout<<"("<<v.x()<<", "<<v.y()<<", "<<v.z()<<")";
}

std::vector<FitResult3D>& eraseUnmatched(std::vector<FitResult3D>& fits, const std::vector<bool>& fitIsMatched) {
	int iFit=0;
	for(auto it=fits.begin(); it!=fits.end();) {
		if(!fitIsMatched[iFit++]) {
			it=fits.erase(it);
		} else {
			++it;
		}
	}
	return fits;
}

std::pair<HoughTransformer::HitCluster, HoughTransformer::HitCluster> splitCluster(const HoughTransformer::HitCluster& cluster, std::function<bool(const PositionHit&)> predicate){
	HoughTransformer::HitCluster isTrue, isFalse;
	for(auto& h : cluster) {
		if(predicate(h)){
			isTrue.push_back(h);
		} else {
			isFalse.push_back(h);
		}
	}
	return {isTrue,isFalse};
}

void TrackCombiner::processTracks() {
	const TVector3& timepixShift=alignment.relativeAlignment.shift;
	const TVector3& rotationCOM=alignment.relativeAlignment.getCOM();

	nTelescopeTriggers=0;
	telescopeFitter.getEntry(0);
	for(int telescopeEntryNumber=0,tpcEntryNumber=0;
//			i<100000
			;) {

		triggerStatusHistogram.reset();

		// Get Entry and match trigger Numbers
		auto matchStatus=getAndMatchEntries(telescopeEntryNumber,tpcEntryNumber);
//		printTriggers(telescopeEntryNumber,tpcEntryNumber);
//		if( cin.get()=='q') break;
		if( matchStatus == MatchResult::end ) break;
		else if( matchStatus == MatchResult::noMatch) continue;

		// Fit tracks

		//telescope
		telescopeFits.clear();
		auto telescopeHits=telescopeFitter.getSpaceHits();
		if( !telescopeFitter.passEvent(telescopeHits) ) { replaceStatus(1, "less than 4 planes hit in telescope", tpcEntryNumber); continue;} //minimal 4 planes!
		telescopeHits=telescopeFitter.rotateAndShift(telescopeHits);
		for(auto&v:telescopeHits) for(auto&h:v) {
			h.RotatePosition(-savedSlopes.first, {mimosa.getCentre().first,mimosa.getCentre().second,0}, {0,1,0} );
			h.RotatePosition(savedSlopes.second, {mimosa.getCentre().first,mimosa.getCentre().second,0}, {1,0,0} );
		}
		auto telescopeClusters = telescopeFitter.houghTransform(telescopeHits);
		for( auto& cluster : telescopeClusters) {
			if(cluster.size()<4 or cluster.getNPlanesHit()<=3) { replaceStatus(1, "less than 4 planes hit in telescope", tpcEntryNumber); continue;}
			auto fit=regressionFit3d(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			auto residuals=calculateResiduals(cluster, fit);
			cluster=cutOnResiduals(cluster, residuals, 0.5 /*mm*/);
			if(cluster.size()<4  or cluster.getNPlanesHit()<=3) { replaceStatus(1, "less than 4 planes hit in telescope", tpcEntryNumber); continue;}
			fit=regressionFit3d(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			telescopeFits.push_back(fit);
		}
		if(telescopeFits.empty()) { replaceStatus(2, "All telescope clusters failed fit", tpcEntryNumber); continue; }

//		cout<<"telescope passed!"<<endl;

		//timepix
		tpcFits.clear();
		auto tpcHits=tpcFitter.getSpaceHits();
		if( !tpcFitter.passEvent(tpcHits) ) { replaceStatus(3, "Less than 20 hits in tpc", tpcEntryNumber); continue; }
		tpcHits=tpcFitter.rotateAndShift(tpcHits);
		tpcHits=alignment.timeWalkCorrection.correct(tpcHits);
		auto tpcHitsInTimePixFrame=tpcHits;//copy hits before rotation
		for(auto& h: tpcHits) {
			h.y=-h.y;
			h.RotatePosition(alignment.relativeAlignment.angle[1], rotationCOM, {0,1,0});
			h.RotatePosition(alignment.relativeAlignment.angle[0], rotationCOM, {1,0,0});
			h.RotatePosition(alignment.relativeAlignment.angle[2], rotationCOM, {0,0,1});
			h.SetPosition(h.getPosition() + timepixShift);
		}
		auto tpcClusters = tpcFitter.houghTransform(tpcHits);
		if(tpcClusters.size()>1) { auto mes="More than one cluster in tpc"; replaceStatus(5, mes, tpcEntryNumber); continue; };
		vector<HoughTransformer::HitCluster> tpcFittedClusters; //todo: replace with actual cluster for removing hits from cluster
		for( auto& cluster : tpcClusters ) {
			if(cluster.size()<2) continue;
			auto fit=regressionFit3d(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			auto residuals=calculateResiduals(cluster, fit);
//			cout<<cluster.size();
			cluster=cutOnResiduals(cluster, residuals, 1 /*mm*/);
//			cout<<" - "<<cluster.size()<<"\n";
			fit=regressionFit3d(cluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}
			tpcFits.push_back(fit);
			tpcFittedClusters.push_back(cluster);
		}
		if(tpcFits.empty()) { replaceStatus(4, "all tpc clusters failed fit", tpcEntryNumber); continue; }
		if(none_of(tpcFits.begin(), tpcFits.end(), [this](const FitResult3D& f) { return goesThroughTimepix(f, alignment) == 2; } )) {
			replaceStatus(5, "tpc fit leaves or enters on the side", tpcEntryNumber); continue;
		}

//		cout<<"timepix passed!"<<endl;

		//make new entry
		BufferedTreeFiller::TreeEntry treeEntry;

		//match fits and clusters + calculate residuals
		int nmatched=0;
		bool atLeastOneThroughDetector=false;
		vector<FitResult3D> telescopeTPCLines, combinedFits;
		std::vector<bool> tpcFitIsMatched(tpcFits.size()), telescopeFitIsMatched(telescopeFits.size());
		for(unsigned iFit=0; iFit<tpcFits.size();++iFit) {
			const auto& tpcFit=tpcFits[iFit];
			for(unsigned jFit=0; jFit<telescopeFits.size(); jFit++) {
				const auto& telescopeFit = telescopeFits[jFit];

				//check if going through detector
				if( goesThroughTimepix(telescopeFit, alignment) ) {
					atLeastOneThroughDetector=true;
				}

				if( fabs( tpcFit.xAt(timepixShift.z()+timePixChip.zmean())-telescopeFit.xAt(timepixShift.z()+timePixChip.zmean()) ) < 1
					&& fabs( tpcFit.yAt(timepixShift.z()+timePixChip.zmean())-telescopeFit.yAt(timepixShift.z()+timePixChip.zmean()) ) < 1	  ) {

					if(telescopeFitIsMatched[jFit]) { std::cerr<<"telescope fit is matched to 2 timepix clusters!"<<std::endl; }

					TVector3 average=tpcFittedClusters.at(iFit).getAveragePosition();

					//draw line from telescope at z=0 to tpc centre
					FitResult3D telescopeTPCLine{
							FitResult2D{
								(average.x()-telescopeFit.xAt(0))/average.z(),
								telescopeFit.xAt(0),
								{0.,0.,0.}
							},
							FitResult2D{
								(average.y()-telescopeFit.yAt(0))/average.z(),
								telescopeFit.yAt(0),
								{0.,0.,0.}
							}
					};
					telescopeTPCLines.push_back(telescopeTPCLine);

					//make alternative fit using an extra point in telescope
//					auto splitTpcCluster=splitCluster(tpcFittedClusters.at(iFit), [](const PositionHit& h) {
//						return h.column<128;
////						return (h.column+h.row)%2;
//					});
////					cout<<splittedTpcCluster.first.size()<<" - "<<splittedTpcCluster.second.size()<<"\n";
//					treeEntry.nfitted=splitTpcCluster.first.size();
//					treeEntry.nresiduals=splitTpcCluster.second.size();
//					if(splitTpcCluster.first.size()<1 or splitTpcCluster.second.size()<1) {
//						cout<<"split cluster has no hits!\n";
//						continue;
//					}
//
//					PositionHit lastPlaneCrossing( telescopeFit.xAt(0), telescopeFit.yAt(0), 0 );
//					lastPlaneCrossing.error2x=lastPlaneCrossing.error2y=1E-4;//=0.01mm
//					splitTpcCluster.first.add( lastPlaneCrossing );
//					auto combinedFit=regressionFit3d(splitTpcCluster.first);
//					combinedFits.push_back(combinedFit);
//					auto residuals=calculateResiduals( splitTpcCluster.second, combinedFit);

					auto residuals=calculateResiduals(tpcFittedClusters.at(iFit), telescopeFit);
					//rotate back to frame of timepix
					for(auto& r:residuals) {
						auto v=r.getVector();
						v.Rotate(-alignment.relativeAlignment.angle[2], {0,0,1});
						v.Rotate(-alignment.relativeAlignment.angle[0], {1,0,0});
						v.Rotate(-alignment.relativeAlignment.angle[1], {0,1,0});
						r.setVector(v);
					}
					treeEntry.tpcResiduals.emplace_back( residuals.begin(), residuals.end() );//construct new vector with entries just for this cluster in tpcResiduals

//					average=splitTpcCluster.second.getAveragePosition();
					auto& fit=telescopeFit; //telescopeFit,combinedFit
					treeEntry.dyz=(average.z()*fit.YZ.slope+fit.YZ.intercept-average.y())/sqrt(1+fit.YZ.slope*fit.YZ.slope);
					treeEntry.dxz=(average.z()*fit.XZ.slope+telescopeFit.XZ.intercept-average.x())/sqrt(1+fit.XZ.slope*fit.XZ.slope);

					++nmatched;
					telescopeFitIsMatched[jFit]=true;
					tpcFitIsMatched[iFit]=true;
					break;
				}
			}
		}

		//remove unmatched clusters
		bool removeUnmatched=true;
		if(removeUnmatched) {
			tpcFits=eraseUnmatched(tpcFits, tpcFitIsMatched);
			telescopeFits=eraseUnmatched(telescopeFits, telescopeFitIsMatched);

			if( telescopeFits.empty() || tpcFits.empty() ) {
				auto message= atLeastOneThroughDetector ? "Telescope and tpc fits do not match" : "Telescope fit missed tpc";
				replaceStatus(10, message, tpcEntryNumber);
				if(displayEvent) cout<<"telescope and tpc fits do not match!"<<endl;
				continue;
			}
		}

		if(not atLeastOneThroughDetector) cerr<<"fit did not go through detector but did match!?"<<endl;

		//display event
		if( displayEvent and tpcHits.size()<300 and tpcHits.size()*1./tpcFittedClusters.front().size()>0.95) {
			drawEvent(tpcHits, telescopeFits);
//			drawEvent(tpcFittedClusters.front(), telescopeTPCLines);
			if( telescopeFitter.processDrawSignals()  ) break;
		}

		//check if tpcEntry already has a matching fit
		if( find(tpcEntryHasMatchingFit.begin(), tpcEntryHasMatchingFit.end(), tpcEntryNumber )!=tpcEntryHasMatchingFit.end() ) {
			if( treeBuffer.isInBuffer(tpcEntryNumber)){
				treeBuffer.removeFromBuffer(tpcEntryNumber); //and remove other entry from buffer if so
			} else {
				cout<<"Entry that should be removed from buffer, was already removed.\n"
				"(The entry could have matched to 3 clusters, but this message can also indicate errors)\n";
			}
			auto mes="Tpc entry already has a matching cluster"; replaceStatus(21, mes, tpcEntryNumber);
			continue;
		}

//		cout<<"Success!\n";
		replaceStatus(20, "Successful", tpcEntryNumber);

		tpcEntryHasMatchingFit.push_back(tpcEntryNumber);
		//cleanup old entry info
		while(not tpcEntryHasMatchingFit.empty() and tpcEntryHasMatchingFit.front()<tpcEntryNumber-(telescopeFitter.triggerNumberEnd-telescopeFitter.triggerNumberBegin) ) tpcEntryHasMatchingFit.pop_front();


		for(unsigned iClust=0; iClust<tpcFittedClusters.size(); ++iClust) {
			if(tpcFitIsMatched[iClust])	treeEntry.tpcClusterSize.push_back( tpcFittedClusters[iClust].size());//count cluster size of only fitted clusters
		}
		treeEntry.ntpcHits=tpcHits.size();
		treeEntry.ntelescopeHits=0;
		for(auto& v : telescopeHits ) treeEntry.ntelescopeHits+=v.size();
		treeEntry.telescopeFits=telescopeFits;
		treeEntry.tpcFits=tpcFits;
		treeBuffer.placeInBuffer(tpcEntryNumber, treeEntry);
		int oldestRelevantTpcEntryNumber=tpcEntryNumber-(telescopeFitter.triggerNumberEnd-telescopeFitter.triggerNumberBegin);
//		cout<<"oldestRelevantTpcEntryNumber: "<<oldestRelevantTpcEntryNumber<<"\n";
		treeBuffer.emptyBufferUpTo(oldestRelevantTpcEntryNumber);
		timepixStatusKeepers.writeBufferUpTo(oldestRelevantTpcEntryNumber, [](StatusKeeper&s){s.reset();});

		double telescopeSeconds=telescopeFitter.timestamp/40.E6, tpcSeconds=tpcFitter.timestamp/4096.*25E-9;
		static double firstTimeDifference=telescopeSeconds-tpcSeconds;

	}

	cout<<"highest telescope trigger (number of unique) "<< telescopeFitter.triggerNumberEnd<<" ("<<nTelescopeTriggers<<")"<<endl;
	cout<<"number of entries in timepix "<<tpcFitter.nEvents<<endl;

	treeBuffer.emptyBuffer();
	cout<<"entries in tree "<<treeBuffer.GetEntries()<<endl;


}

TrackCombiner::MatchResult TrackCombiner::getAndMatchEntries(
		int& telescopeEntry,
		int& tpcStartEntry) {

	//if previous frame did not increase and this frame did not increase, there is no related trigger to this frame
//	if(previous2TriggerNumberBegin==telescopeFitter.triggerNumberEnd) { cout<<"unrelated!"<<endl; }
//		triggerStatusHistogram.Fill("Unrelated telescope frame",1); return MatchResult::noMatch; }

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
//		triggerStatusHistogram.Fill("Trigger numbers do not match", 1);

		frameStatusHistogram.reset();
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

void TrackCombiner::loadAlignment(std::string alignmentFile) {
	alignment=Alignment(alignmentFile);
}

void TrackCombiner::saveAlignment(std::string alignmentFile) {
	//calculate alignment
	alignment.timeWalkCorrection.calculate(&treeBuffer.getTree());
	alignment.relativeAlignment.calculate(&treeBuffer.getTree());

	alignment.saveToFile(alignmentFile);
}

void TrackCombiner::printTriggers(int telescopeEntry, int tpcEntry) const {
	const int printEveryN=1;
	if( !(telescopeEntry%printEveryN) ) {
		cout<<"entry: "<<telescopeEntry<<"/"<<telescopeFitter.nEvents<<" ";
		cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd;
		cout<<" timepix triggerNumber: "<<tpcFitter.triggerNumber<<"="<<(tpcFitter.triggerNumber+triggerOffset) % 32768<<" in entry "<<tpcEntry<<endl;
	}
}
