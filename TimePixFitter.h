/*
 * TimePixFitter.h
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#ifndef TimePixFitter_H_
#define TimePixFitter_H_

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

class TimePixFitter {
public:
	TimePixFitter(std::string inputfile, const DetectorConfiguration& detector);
	virtual ~TimePixFitter();

	int makeMask(double ntimesThreshold=1e4);
	void fitTracks( std::string outputfilename );

	std::pair<double,double> getMeans(); //not constant because adds fits!
	const std::pair<double,double>& getShift() const;
	double getRotation(); //not constant because adds fits!
	const double& getAngle() const;
	std::pair<double, double> getSlopes() const;

	void setShift( const std::pair<double,double>& shift);
	void addToShift( const std::pair<double,double>& shift );
	void setAngle( double angles);
	void addToAngle( double angles);
	void setSlopes( std::pair<double, double> slopes);
	int getEntry(int iEvent);

	bool displayEvent=false;
	bool makeTrackHistograms=false;
	bool recalculateCOM=true; //centre of mass
	bool constructLineParallelToZ=false;

	double maxResidual=0.2;

	std::function<bool(const PositionHit&)> selectHitForRefit = [](const PositionHit&){return true;}; //select all rawHits by default

protected:
	std::vector<PositionHit> getSpaceHits();
	std::vector<PositionHit>&  rotateAndShift(
			std::vector<PositionHit>& spaceHit);

private:
	TFile* file;
	TTree* hitTable;
	const std::vector<TimePixHit>* rawHits=nullptr;

	const DetectorConfiguration& detector;
	HoughTransformer houghTransform;
	std::unique_ptr<ResidualHistogrammer> residualHistograms;
	std::unique_ptr<TrackHistogrammer> trackHistograms;

	bool passEvent(const std::vector<PositionHit>& spaceHit ) const;	//return true if the event is passed

	pixelMask mask;
	std::pair<double,double> shift;
	double angle;
	std::pair<double,double> hitsCentre, averageResidualFromSum;
	double rotationZFromSum;

	double slope1FromSum, slope2FromSum;


};

#endif /* TimePixFitter_H_ */
