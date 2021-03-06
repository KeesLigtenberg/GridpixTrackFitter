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
#include "linearRegressionFit.h"
#include "Alignment.h"

//forward declaration
class ResidualHistogrammer;
class TrackHistogrammer;

class TimePixFitter {
public:
	TimePixFitter(std::string inputfile, const DetectorConfiguration& detector);
	virtual ~TimePixFitter();

	int makeMask(double ntimesThreshold=1e4);
	void fitTracks( std::string outputfilename );

	std::pair<double,double> getMeans(); //not a constant function because adds fits!
	const std::pair<double,double>& getShift() const;
	double getRotation(); //not constant because adds fits!
//	const double& getAngle() const;
	std::pair<double, double> getSlopes() const;

	void setShift( const std::pair<double,double>& shift);
	void addToShift( const std::pair<double,double>& shift );
	void setSlopes( std::pair<double, double> slopes);
	int getEntry(int iEvent);

	void drawEvent(const std::vector<PositionHit>& spaceHit,
			const std::vector<FitResult3D>& fits);

	//process controls:
	bool displayEvent=false;
	bool makeTrackHistograms=true;
	bool recalculateCOM=true; //centre of mass
	bool constructLineParallelToZ=false;

	double maxResidual=1;

	std::function<bool(const PositionHit&)> selectHitForRefit = [](const PositionHit&){return true;}; //select all rawHits by default

	HoughTransformer houghTransform;

protected:
	std::vector<PositionHit> getSpaceHits();
	std::vector<PositionHit> getSpaceHitsWithCrossTalk(const TimeWalkCorrector&);
	std::vector<PositionHit>&  rotateAndShift(
			std::vector<PositionHit>& spaceHit);

private:
	TFile* file;
	TTree* hitTable;
	long long nEvents=0;
	const std::vector<TimePixHit>* rawHits=nullptr;
	long long triggerNumber=0;
	unsigned long long timestamp=0;

	const DetectorConfiguration& detector;
	std::unique_ptr<ResidualHistogrammer> residualHistograms;
	std::unique_ptr<TrackHistogrammer> trackHistograms;
	TH1D* numberOfHitsInTrack;
	TH1D* ToTOfHitsInTrack;


	bool passEvent(const std::vector<PositionHit>& spaceHit ) const;	//return true if the event is passed

	pixelMask mask;
	std::pair<double,double> shift;
	std::pair<double,double> hitsCentre, averageResidualFromSum;
	double rotationZFromSum;

	double slope1FromSum, slope2FromSum;

//friend function!
	friend class TrackCombiner;
};

#endif /* TimePixFitter_H_ */
