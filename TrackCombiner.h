/*
 * TrackCombiner.h
 *
 *  Created on: Sep 21, 2017
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_H_
#define TRACKCOMBINER_H_

#include <string>
#include <deque>

#include "EntryBuffer.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"

#include "linearRegressionFit.h"
#include "Alignment.h"
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

//helper class to fill tree buffered.
class BufferedTreeFiller {
public:
	BufferedTreeFiller() { setTreeBranches(); };

	struct TreeEntry {
		vector<FitResult3D> telescopeFits{};
		vector<FitResult3D> tpcFits{};
		vector< vector<HitEntry> > tpcResiduals{};
		vector<int> tpcClusterSize{}; //additional information because tpcResidualsl[i].size() is not easily accessible
		int ntpcHits=0, ntelescopeHits=0;
		double dxz=0., dyz=0.;
		int nfitted=0, nresiduals=0;
	};

	void Write();
	int64_t GetEntries() { return fitResultTree.GetEntriesFast(); }
	void SetTreeDirectory(TFile* f) { fitResultTree.SetDirectory(f); }
	void placeInBuffer(int entryNumber, const TreeEntry& entry) {buffer.placeInBuffer(entryNumber, entry);}
	void emptyBufferUpTo(int entryNumber) { buffer.writeBufferUpTo(entryNumber, [this](TreeEntry&t){this->Fill(t);}); }
	void emptyBuffer() { buffer.writeBuffer([this](TreeEntry&t){this->Fill(t);}); }
	void removeFromBuffer(int entryNumber) { buffer.removeFromBuffer(entryNumber); }
	bool isInBuffer(int tpcEntryNumber) const { return buffer.isInBuffer(tpcEntryNumber); }
	TTree& getTree() { return fitResultTree; }

private:
	TreeEntry currentEntry;
	EntryBuffer<int, TreeEntry> buffer;

	TTree fitResultTree{ "fitResults", "Tree with telescope and timepix fit results"};

	void Fill(const TreeEntry&);
	void setTreeBranches();
};


class TrackCombiner {
public:
	TrackCombiner(std::string mimosaInput, std::string timepixInput, const DetectorConfiguration& telescope=mimosa, const DetectorConfiguration& tpc=timePixChip);
	virtual ~TrackCombiner();

	void openFile(std::string filename="fitResults.root");
	void setTriggerOffset(int offset) { triggerOffset=offset; };
	void setDisplayEvent(bool doDisplay=true) { displayEvent=doDisplay; };

	void loadAlignment(std::string alignmentFile);
	void saveAlignment(std::string alignmentFile);

	void processTracks();

	template <class T>
	void drawEvent(const T& hits,
			const std::vector<FitResult3D>& fits);

	Alignment alignment;

	bool doSplitForResiduals=false;

private:
	enum class MatchResult { match, noMatch, end };
	MatchResult getAndMatchEntries(int& telescopeEntry, int& tpcStartEntry);
	void printTriggers(int telescopeEntry, int tpcEntry) const;

	const DetectorConfiguration& telescope;
	const DetectorConfiguration& tpcDetector;

	trackFitter telescopeFitter;
	TimePixFitter tpcFitter;

	int triggerOffset=-1;
	int previousTriggerNumberBegin=0, previousTriggerNumberEnd=0;
	bool hadFirstMatch=false;
	int timepixEntryFirstMatch=0;
//	int previous2TriggerNumberBegin=0;
	int nTelescopeTriggers=0;
	bool displayEvent=false;
	bool correctToTByCol=true;
	bool correctTimewalk=true;
	bool onlyUseInArea=false;
	std::deque<int> tpcEntryHasMatchingFit{};

	std::unique_ptr<TFile, std::function<void(TFile*)> > outputFile{nullptr, [](TFile* f) { f->Close(); delete f;} };

	//for tree
	BufferedTreeFiller treeBuffer;

	vector<FitResult3D> telescopeFits{};
	vector<FitResult3D> tpcFits{};

	struct StatusKeeper{
		StatusKeeper(std::string name) : statusHistogram( new TH1D( (name+"Status").c_str(), ("status of "+name).c_str(), 1,0,1) ) {};
		StatusKeeper(std::shared_ptr<TH1D> h) : statusHistogram(h) {};
		StatusKeeper(const TrackCombiner::StatusKeeper&) = default;
		virtual ~StatusKeeper() {};
		int priority=0; std::string message="";
		void replace(int messagePriority, std::string newMessage) {
			if(messagePriority>priority) { message=newMessage; priority=messagePriority;}
		}
		void reset() {
			if(priority) statusHistogram->Fill(message.c_str(),1);
			priority=0;
		}
		void Write() { statusHistogram->LabelsDeflate(); statusHistogram->Write(); };
		std::shared_ptr<TH1D> statusHistogram;
	} frameStatusHistogram{"frame"}, triggerStatusHistogram{"trigger"}, timepixStatusHistogram{"timepixTrigger"};
	EntryBuffer<int, StatusKeeper> timepixStatusKeepers; //special buffered statusKeeper function

	//one function to fill all statushistograms at once.
	void replaceStatus(int priority, std::string message, int tpcEntryNumber) {
		for(auto* s : {&frameStatusHistogram, &triggerStatusHistogram} ) s->replace(priority, message);
		if(!timepixStatusKeepers.isInBuffer(tpcEntryNumber)) {
			timepixStatusKeepers.placeInBuffer(tpcEntryNumber, timepixStatusHistogram);
		}
		timepixStatusKeepers.getFromBuffer(tpcEntryNumber).replace(priority, message);
	}

};



#endif /* TRACKCOMBINER_H_ */
