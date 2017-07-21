/*
 * FitTracksTimePix.cpp
 *
 *  Created on: 11 jul. 2017
 *      Author: KeesL
 */

#include <iostream>

#include "TTree.h"

#include "TimePixFitter.h"
#include "PositionHit.h"
#include "HoughTransformer.h"
#include "linearRegressionFit.h"
#include "Hit.h"
#include "getObjectFromFile.h"
#include "ResidualHistogrammer.h"
#include "makeNoisyPixelMask.h"

#if 1 //root?
#include "TimePixFitter.cpp"
#include "linearRegressionFit.cpp"
#include "ResidualHistogrammer.cpp"
#include "makeNoisyPixelMask.cpp"
#endif

using namespace std;

struct TimePixDetectorConfiguration : DetectorConfiguration {
	TimePixDetectorConfiguration() : DetectorConfiguration{
		1, {0}, //nplanes, planeposition
		1, 256, 256 //pixelsize, xpixels, ypixels
	} {};
	virtual double xmin() const {return -200; }
	virtual double xmax() const {return 200; }
	virtual double zmin() const {return 0; }
	virtual double zmax() const {return 256; };
} timePixChip;

void FitTracksTimePix(std::string inputfile) {

	//get tree from file
	TimePixFitter tpcFitter(inputfile,timePixChip);

	tpcFitter.makeMask(1e3);

	tpcFitter.setSlopes( {0.006,-0.154} );
	tpcFitter.houghTransform.minCandidateSize=6;
	tpcFitter.houghTransform.minClusterSize=10;

	tpcFitter.displayEvent=true;

	tpcFitter.fitTracks("timepixHistograms.root");

}


