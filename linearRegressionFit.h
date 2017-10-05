/*
 * linearRegressionFit.h
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#ifndef LINEARREGRESSIONFIT_H_
#define LINEARREGRESSIONFIT_H_

#include <array>

#include <TVector3.h>
#include "TPolyLine3D.h"

#include "HoughTransformer.h"

struct FitResult2D {
	FitResult2D( double slope, double intercept, std::vector<double> error, double interceptz=0 ):
		slope(slope),
		intercept(intercept),
		error(error),
		interceptz(interceptz)
	{}
	FitResult2D() : slope(0), intercept(0), error({0,0,0}), interceptz(0) {}

	double slope, intercept;
	std::vector<double> error; //dslope^2, dslopeintercept, dintercept^2 todo: replace with std::array<3>
	double interceptz;

	double at(double z) const { return intercept+slope*(z-interceptz); }
	double error2At(double z) const { return z*z*error[0]+2*z*error[1]+error[2]; }
	double errorAt(double z) const { return sqrt(error2At(z)); }
	bool isValid() const { return !(std::isnan(slope) || std::isnan(intercept)); }

	FitResult2D makeShifted(double shift, double shiftz) const {
		return { slope, intercept+shift-shiftz*slope, error, /*errors do not change*/ interceptz };
	}
	FitResult2D makeMirror() const {
		return {-slope, -intercept, error, interceptz};
	}

};

struct FitResult3D {
	FitResult3D() : XZ(), YZ() {};
	FitResult3D(FitResult2D XZ, FitResult2D YZ) : XZ(XZ), YZ(YZ) {};

	FitResult2D XZ, YZ;

	void draw(double zmin, double zmax) const;
	bool isValid() const { return XZ.isValid() and YZ.isValid(); }
	double xAt(double z) const { return XZ.at(z); };
	double yAt(double z) const { return YZ.at(z); };

	FitResult3D makeShifted( const TVector3& shift ) const {
		return { XZ.makeShifted(shift.x(), shift.z()), YZ.makeShifted(shift.y(), shift.z())	};
	}
	FitResult3D makeMirrorY() const {
		return {XZ, YZ.makeMirror()};
	}
	FitResult3D makeRotated(double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis ) const;
};


//root dictionary for use in TTree
//#pragma link C++ class std::array<double, 3>+;
#pragma link C++ class std::vector<double>+;
#pragma link C++ class FitResult2D+;
#pragma link C++ class FitResult3D+;
#pragma link C++ class std::vector<FitResult3D>+;

//struct TrackFitResult {
//	TrackFitResult(const SimpleFitResult& fr) : //TODO: check! this was before y <-> z
//		phi(atan(1.)*2 - atan(fr.slope1)),
//		d0(fr.intersept1 * sin(phi)), //had -
//		tanlambda( (tan(fr.slope2 / fabs(fr.slope2)*acos(sqrt(fr.slope1*fr.slope1 + 1) / sqrt(fr.slope1*fr.slope1 + 1 + fr.slope2*fr.slope2)))) ),
//		z0(- (fr.slope1*fr.intersept1*fr.slope2) + (fr.slope1*fr.slope1 + 1) + fr.intersept2)
//	{}
//	double phi, d0, tanlambda, z0;
//};

std::ostream& operator<<(std::ostream& os, FitResult3D& fit);

FitResult3D regressionFit3d(const HoughTransformer::HitCluster& cluster, double interceptz=0);

FitResult3D makeLinesParallelToZ(double x, double y) {
	return {
		FitResult2D{ 0, x, {0,0,0} },
		FitResult2D{ 0, y, {0,0,0} }
	};
}

//TrackFitResult linearRegressionTrackFit(const HoughTransformer::HitCluster& cluster);

struct Residual {
	double x,y,z;
	const PositionHit& h;
	TVector3 getVector() { return TVector3(x,y,z); };
	void setVector(const TVector3& v) { x=v.x(); y=v.y(); z=v.z(); };
};

Residual calculateResidual( const PositionHit& h, const FitResult3D& fit );
std::vector<Residual> calculateResiduals( const HoughTransformer::HitCluster& cluster, const FitResult3D& fit);
TVector3 averageResidual(const std::vector<Residual>& residuals);

HoughTransformer::HitCluster& cutOnResiduals( HoughTransformer::HitCluster& cluster, const std::vector<Residual>& residuals, double maxResidual );

#endif /* LINEARREGRESSIONFIT_H_ */
