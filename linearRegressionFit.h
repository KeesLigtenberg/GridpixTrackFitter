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
	FitResult2D( double slope, double intercept, std::array<double, 3> error, double interceptz=0 ):
		slope(slope),
		intercept(intercept),
		error(error),
		interceptz(interceptz)
	{
	}

	double slope, intercept;
	std::array<double, 3> error; //dslope^2, dslopeintercept, dintercept^2
	double interceptz;

	double at(double z) const { return intercept+slope*(z-interceptz); }
	bool isValid() const { return !(std::isnan(slope) || std::isnan(intercept)); }

	FitResult2D makeShifted(double shift, double shiftz) const {
		return { slope, intercept+shift, error, /*errors do not change*/ interceptz+shiftz };
	}
	FitResult2D makeMirror() const {
		return {-slope, -intercept, error, interceptz};
	}

};

struct FitResult3D {
	FitResult3D(FitResult2D XZ, FitResult2D YZ) : XZ(XZ), YZ(YZ) {};

	FitResult2D XZ, YZ;

	void draw(double zmin, double zmax) const;
	bool isValid() const { return XZ.isValid() and YZ.isValid(); }

	FitResult3D makeShifted( const TVector3& shift ) const {
		return { XZ.makeShifted(shift.x(), shift.z()), YZ.makeShifted(shift.y(), shift.z())	};
	}
	FitResult3D makeMirrorY() const {
		return {XZ, YZ.makeMirror()};
	}
	FitResult3D makeRotated(double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis ) const;
};

void FitResult3D::draw(double zmin, double zmax) const {
	const int npoints=2;
	double x[npoints] = { XZ.at(zmin), XZ.at(zmax)};
	double y[npoints] = { YZ.at(zmin), YZ.at(zmax)};
	double z[npoints] = { zmin, zmax};
	TPolyLine3D l( npoints, z, y, x );
	l.SetLineColor(kOrange+7);
	l.SetLineWidth(2);
	l.DrawClone();
}
FitResult3D FitResult3D::makeRotated(double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis ) const {
	if( fabs(XZ.interceptz-YZ.interceptz) > 1E-10 ) throw "intercepts should be described at the same point in space!";
	//rotate
	TVector3 slope(XZ.slope, YZ.slope, 1), intercept(XZ.intercept, YZ.intercept, XZ.interceptz);
	slope.Rotate(rotation, rotationAxis);
	intercept=RotateAroundPoint(intercept, rotation, rotationPoint, rotationAxis);
	//return result
	return FitResult3D{
		FitResult2D{
			slope.x()/slope.z(),
			intercept.z()*slope.x()/slope.z()+intercept.x(),
			XZ.error,//todo:propagate errors!
			intercept.z()
		},
		FitResult2D{
			slope.y()/slope.z(),
			intercept.z()*slope.y()/slope.z()+intercept.y(),
			YZ.error,//todo:propagate errors!
			intercept.z()
		}
	};
}

struct SimpleFitResult {
	double slope1, intersept1, slope2, intersept2;
	double dslope1, dintersept1, dslope2, dintersept2; //uncertainties
//	const HoughTransformer::HitCluster* hitCluster; solve issue with streamer first
	void draw(double zmin, double zmax) const {
		const int npoints=2;
		double x[npoints] = {slope1 * zmin + intersept1, slope1 * zmax + intersept1};
		double y[npoints] = {slope2 * zmin + intersept2, slope2 * zmax + intersept2};
		double z[npoints] = {zmin, zmax};
		TPolyLine3D l( npoints, z, y, x );
		l.SetLineColor(kOrange+7);
		l.SetLineWidth(2);
		l.DrawClone();
	}
	bool isValid() const {
		return !(std::isnan(slope1) || std::isnan(slope2) || std::isnan(intersept1) || std::isnan(intersept2) );
	}
	double xAt(double z) const { return slope1*z+intersept1; };
	double yAt(double z) const { return slope2*z+intersept2; };
	SimpleFitResult makeShifted( const TVector3& shift ) {
		return SimpleFitResult{
			slope1, intersept1-slope1*shift.z()+shift.x(),
			slope2, intersept2-slope2*shift.z()+shift.y(),
			dslope1, dintersept1, //todo: propagate errors!
			dslope2, dintersept2 //todo:propagate errors!
		};
	}
	SimpleFitResult makeRotated(double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis ) {
		TVector3 slope(slope1, slope2, 1), intercept(intersept1, intersept2, 0);
		slope.Rotate(rotation, rotationAxis);
		for(TVector3* v : {&intercept} ) {
			*v-=rotationPoint;
			v->Rotate(rotation, rotationAxis);
			*v+=rotationPoint;
		}
		return SimpleFitResult{
			slope.x()/slope.z(),
			intercept.z()*slope.x()/slope.z()+intercept.x(),
			slope.y()/slope.z(),
			intercept.z()*slope.y()/slope.z()+intercept.y(),
			dslope1, dintersept1, //todo: propagate errors!
			dslope2, dintersept2 //todo:propagate errors!
		};
	}
	SimpleFitResult makeMirrorY() {
		return SimpleFitResult {
			slope1, intersept1, -slope2, -intersept2,
			dslope1, dintersept1, dslope2, dintersept2
		};
	}
};

//root dictionary for use in TTree
#pragma link C++ class SimpleFitResult+;
#pragma link C++ class std::vector<SimpleFitResult>+;

struct TrackFitResult {
	TrackFitResult(const SimpleFitResult& fr) : //TODO: check! this was before y <-> z
		phi(atan(1.)*2 - atan(fr.slope1)),
		d0(fr.intersept1 * sin(phi)), //had -
		tanlambda( (tan(fr.slope2 / fabs(fr.slope2)*acos(sqrt(fr.slope1*fr.slope1 + 1) / sqrt(fr.slope1*fr.slope1 + 1 + fr.slope2*fr.slope2)))) ),
		z0(- (fr.slope1*fr.intersept1*fr.slope2) + (fr.slope1*fr.slope1 + 1) + fr.intersept2)
	{}
	double phi, d0, tanlambda, z0;
};

std::ostream& operator<<(std::ostream& os, SimpleFitResult& fit);

SimpleFitResult linearRegressionFit(const HoughTransformer::HitCluster& cluster ) ;

TrackFitResult linearRegressionTrackFit(const HoughTransformer::HitCluster& cluster);

struct Residual {
	double x,y,z;
	const PositionHit& h;
	TVector3 getVector() { return TVector3(x,y,z); };
	void setVector(const TVector3& v) { x=v.x(); y=v.y(); z=v.z(); };
};

Residual calculateResidual( const PositionHit& h, const SimpleFitResult& fit );
std::vector<Residual> calculateResiduals( const HoughTransformer::HitCluster& cluster, const SimpleFitResult& fit);
TVector3 averageResidual(const std::vector<Residual>& residuals);

HoughTransformer::HitCluster& cutOnResiduals( HoughTransformer::HitCluster& cluster, const std::vector<Residual>& residuals, double maxResidual );

#endif /* LINEARREGRESSIONFIT_H_ */
