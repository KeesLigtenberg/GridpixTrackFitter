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
#include "linearRegressionFit.h"

class ResidualHistogrammer;
class TrackHistogrammer;

class trackFitter {
public:
	trackFitter(std::string inputfile, const DetectorConfiguration& detector);
	virtual ~trackFitter();

	int makeMask(double ntimesThreshold=1e4);
	void fitTracks( std::string outputfilename );

	std::vector<std::pair<double,double>> getMeanResiduals(); //not constant because adds fits!
	const std::vector<std::pair<double,double>>& getShifts() const;
	std::vector<double> getRotations(); //not constant because adds fits!
	const std::vector<double>& getAngles() const;
	std::pair<double, double> getSlopes() const;

	void setShifts( const std::vector<std::pair<double,double>>& shifts);
	void addToShifts( const std::vector<std::pair<double,double>>& shifts );
	void setAngles( const std::vector<double>& angles);
	void addToAngles( const std::vector<double>& angles);
	void setSlopes( std::pair<double, double> slopes);
	void setCentres( const std::vector<std::pair<double,double>>& COMs);
	int getEntry(int iEvent);

	void saveAlignment(std::string outputfile);
	void drawEvent(const std::vector<std::vector<PositionHit> >& spaceHit,
			const std::vector<SimpleFitResult>& fits);
	static bool processDrawSignals();

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
	long long nEvents=0;
	const std::vector<std::vector<Hit>>* mimosaHit=nullptr;
	unsigned short triggerNumberBegin=0, triggerNumberEnd=0;
	unsigned int timestamp=0;

	const DetectorConfiguration& detector;
	HoughTransformer houghTransform;
	std::unique_ptr<ResidualHistogrammer> residualHistograms;
	std::unique_ptr<TrackHistogrammer> trackHistograms;

	bool passEvent( const std::vector<std::vector<PositionHit> >& spaceHit ) const;	//return true if the event is passed

	std::vector<pixelMask> mask;
	std::vector<std::pair<double,double>> shifts;
	std::vector<double> angles;
	std::vector<std::pair<double,double>> hitsCentre, averageResidualFromSum;
	std::vector<double> rotationZFromSum;

	double slope1FromSum, slope2FromSum;

//friend function!
	friend double CombineTracks(std::string mimosaInput, std::string timepixInput, int offset, bool displayEvent);
};

#endif /* TRACKFITTER_H_ */
