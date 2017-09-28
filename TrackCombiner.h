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
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other! e.g -1 for our data!
};

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

	const DetectorConfiguration& telescope;
	const DetectorConfiguration& tpcDetector;

	trackFitter telescopeFitter;
	TimePixFitter tpcFitter;

	int triggerOffset=-1;
	bool displayEvent=false;

	std::unique_ptr<TFile> outputFile{};
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
