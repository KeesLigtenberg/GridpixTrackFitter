/*
 * linearRegressionFit.h
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#ifndef LINEARREGRESSIONFIT_H_
#define LINEARREGRESSIONFIT_H_

#include "HoughTransformer.h"
#include "TPolyLine3D.h"

struct SimpleFitResult {
	double slope1, intersept1, slope2, intersept2;
	void draw(double zmin, double zmax) {
		const int npoints=2;
		double x[npoints] = {slope1 * zmin + intersept1, slope1 * zmax + intersept1};
		double y[npoints] = {slope2 * zmin + intersept2, slope2 * zmax + intersept2};
		double z[npoints] = {zmin, zmax};
		TPolyLine3D l( npoints, x, y, z );
//		std::cout<<"drawing fit "<<std::endl;
//		l.SavePrimitive(std::cout);
		l.DrawClone();
	}
};
struct TrackFitResult {
	TrackFitResult(const SimpleFitResult& fr) : //FIXME: check! this was before y <-> z
		phi(atan(1.)*2 - atan(fr.slope1)),
		d0(- fr.intersept1 * sin(phi)),
		tanlambda( (tan(fr.slope2 / fabs(fr.slope2)*acos(sqrt(fr.slope1*fr.slope1 + 1) / sqrt(fr.slope1*fr.slope1 + 1 + fr.slope2*fr.slope2)))) ),
		z0(- (fr.slope1*fr.intersept1*fr.slope2) + (fr.slope1*fr.slope1 + 1) + fr.intersept2)
	{}
	double phi, d0, tanlambda, z0;
};

SimpleFitResult linearRegressionFit(const HoughTransformer::HitCluster& cluster ) ;

TrackFitResult linearRegressionTrackFit(const HoughTransformer::HitCluster& cluster);

struct Residual { double x,y,z; const PositionHit& h; };
Residual calculateResidual( const PositionHit& h, const SimpleFitResult& fit );
std::vector<Residual> calculateResiduals( const HoughTransformer::HitCluster& cluster, const SimpleFitResult& fit);

HoughTransformer::HitCluster& cutOnResiduals( HoughTransformer::HitCluster& cluster, const std::vector<Residual>& residuals, double maxResidual );

#endif /* LINEARREGRESSIONFIT_H_ */