/*
 * trackFitter.h
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#ifndef TRACKFITTER_H_
#define TRACKFITTER_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Hit.h"
#include "HoughTransformer.h"
#include "PositionHit.h"
#include "makeNoisyPixelMask.h"

class ResidualHistogrammer;
class TrackHistogrammer;

class trackFitter {
public:
	trackFitter(std::string inputfile, const DetectorConfiguration& detector);
	virtual ~trackFitter();

	int makeMask(double ntimesThreshold=1e4);
	void fitTracks( std::string outputfilename );

	std::vector<std::pair<double,double>> getMeans(); //not constant because adds fits!
	const std::vector<std::pair<double,double>>& getShifts() const;
	std::vector<double> getRotations(); //not constant because adds fits!
	const std::vector<double>& getAngles() const;
	std::pair<double, double> getSlopes() const;

	void setShifts( const std::vector<std::pair<double,double>>& shifts);
	void addToShifts( const std::vector<std::pair<double,double>>& shifts );
	void setAngles( const std::vector<double>& angles);
	void addToAngles( const std::vector<double>& angles);
	void setSlopes( std::pair<double, double> slopes);
	int getEntry(int iEvent);

	bool displayEvent=false;
	bool makeTrackHistograms=false;
	bool recalculateCOM=true; //centre of mass
	bool constructLineParallelToZ=false;

	double maxResidual=0.2;

	std::function<bool(const PositionHit&)> selectHitForRefit = [](const PositionHit&){return true;}; //select all hits by default

protected:
	std::vector<std::vector<PositionHit> > getSpaceHits();
	std::vector<std::vector<PositionHit> >&  rotateAndShift(
			std::vector<std::vector<PositionHit> >& spaceHit);

private:
	TFile* file;
	TTree* hitTable;
	const std::vector<std::vector<Hit>>* mimosaHit=nullptr;

	const DetectorConfiguration& detector;
	HoughTransformer houghTransform;
	std::unique_ptr<ResidualHistogrammer> residualHistograms;
	std::unique_ptr<TrackHistogrammer> trackHistograms;

	bool passEvent( std::vector<std::vector<PositionHit> > spaceHit ) ;	//return true if the event is passed

	std::vector<pixelMask> mask;
	std::vector<std::pair<double,double>> shifts;
	std::vector<double> angles;
	std::vector<std::pair<double,double>> hitsCentre, averageResidualFromSum;
	std::vector<double> rotationZFromSum;

	double slope1FromSum, slope2FromSum;


};

#endif /* TRACKFITTER_H_ */
