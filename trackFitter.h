/*
 * trackFitter.h
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#ifndef TRACKFITTER_H_
#define TRACKFITTER_H_

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>

#include "TH1.h"
#include "TTree.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "makeNoisyPixelMask.h"
#include "Hit.h"
#include "PositionHit.h"
#include "HoughTransformer.h"
#include "transformHits.h"

#include "DetectorConfiguration.h"

class trackFitter {
public:
	trackFitter(std::string inputfile, const DetectorConfiguration& detector);
	virtual ~trackFitter();

	int makeMask(double ntimesThreshold=1e4);
private:
	TFile* file;
	TTree* hitTable;
	const std::vector<std::vector<Hit>>* mimosaHit=nullptr;

	const DetectorConfiguration& detector;
	HoughTransformer houghTransform;
	std::unique_ptr<ResidualHistogrammer> residualHistograms;

	std::vector<pixelMask> mask;
	std::vector<std::pair<double,double>> shifts;
	std::vector<double> angles;

	bool passEvent( std::vector<std::vector<PositionHit> > spaceHit ) ;	//return true if the event is passed
	void fitTracks( std::string outputfilename );

};

#endif /* TRACKFITTER_H_ */
