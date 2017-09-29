/*
 * TrackCombiner.h
 *
 *  Created on: Sep 21, 2017
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_H_
#define TRACKCOMBINER_H_

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

class TrackCombiner {
public:
	TrackCombiner(std::string mimosaInput, std::string timepixInput, const DetectorConfiguration& telescope=mimosa, const DetectorConfiguration& tpc=timePixChip);
	virtual ~TrackCombiner();

	void openFile(std::string filename="fitResults.root");
	void setTriggerOffset(int offset) { triggerOffset=offset; };
	void setDisplayEvent(bool doDisplay=true) { displayEvent=doDisplay; };


	void processTracks();
	void drawEvent();

private:
	void setTreeBranches();
	enum class MatchResult { match, noMatch, end };
	MatchResult getAndMatchEntries(int& telescopeEntry, int& tpcStartEntry);
	void printTriggers(int telescopeEntry, int tpcEntry) const;

	const DetectorConfiguration& telescope;
	const DetectorConfiguration& tpcDetector;

	trackFitter telescopeFitter;
	TimePixFitter tpcFitter;

	int triggerOffset=-1;
	int previousTriggerNumberBegin=0;
	bool hadFirstMatch=false;
	int timepixEntryFirstMatch=0;
//	int previous2TriggerNumberBegin=0;
	int nTelescopeTriggers=0;
	bool displayEvent=false;

	std::unique_ptr<TFile, std::function<void(TFile*)> > outputFile{nullptr, [](TFile* f) { f->Close(); delete f;} };
	TTree fitResultTree{ "fitResults", "Tree with telescope and timepix fit results"};

	//for tree
	vector<SimpleFitResult> telescopeFits{};
	vector<SimpleFitResult> tpcFits{};
	vector< vector<HitEntry> > tpcResiduals{};
	int ntpcHits=0, ntelescopeHits=0;
	vector<int> tpcClusterSize{};
	double dxz=0., dyz=0.;
	TH1D triggerStatus{"triggerStatus", "status of trigger", 1, 0, 1};

};

#endif /* TRACKCOMBINER_H_ */
