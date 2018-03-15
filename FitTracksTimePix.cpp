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
#include "TrackFitter.cpp"
#endif

#include "testBeamSetup.h"

using namespace std;

void FitTracksTimePix(std::string inputfile) {

	//get tree from file
	TimePixFitter tpcFitter(inputfile,timePixChip);

	tpcFitter.makeMask(); //mask noisy pixels

	//configure bin algorithm. Slope is needed to
	tpcFitter.setSlopes( {0.006,-0.26} );
	tpcFitter.houghTransform.minCandidateSize=6;
	tpcFitter.houghTransform.minClusterSize=10;

	tpcFitter.displayEvent=true;

	tpcFitter.fitTracks("timepixHistograms.root");

}


