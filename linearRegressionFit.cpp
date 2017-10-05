/*
 * linearRegressionFit.cpp
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#include <iostream>

#include "linearRegressionFit.h"


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

FitResult2D regressionXZ(const HoughTransformer::HitCluster& cluster, double interceptz=0) {
    double sumX = 0;
    double sumZ = 0;
    double sumXZ = 0;
    double sumZsquare = 0;  // = Sum (Z^2)
    double sumW = 0;

    for(const auto& h : cluster) {
    	double errorx2=1;
    	double hiz=h.z-interceptz;
		sumX += h.x/errorx2;
		sumZ += hiz/errorx2;
		sumXZ += h.x*hiz/errorx2;
		sumZsquare += hiz*hiz/errorx2;
		sumW+=1/errorx2;
    }

    // z <-> y from LinearTrackRegression::doRegressionX of marlinTPC
    /** gives back parameters for \n
     *  x = slope1 * z + intersept1\n
     *  y = slope2 * z + intersept2\n
     *
     */

    double denominator=(sumZ * sumZ - sumW * sumZsquare);
    if(std::fabs(denominator)<1E-20){
    	std::cerr<<"error: (sumZ * sumZ - ntot * sumZsquare)<1E-20"<<std::endl;
    	std::cerr<<"sumZ="<<sumZ<<" sumZsquare="<<sumZsquare<<" ntot="<<sumW<<std::endl;
    	std::cerr<<cluster.size()<<" hits on "<<cluster.getNPlanesHit()<<" planes"<<std::endl;
    	throw "(sumZ * sumZ - ntot * sumZsquare)<1E-20";
    }

    double slope1     = (sumX * sumZ - sumW * sumXZ) / denominator;
    double intersept1 = (sumZ * sumXZ - sumZsquare * sumX) / denominator;

    double sigmaSlope2=sumZsquare/denominator;
    double sigmaIntercept2=sumW/denominator;
    double sigmaSlopeIntercept=-sumZ/denominator;
    std::vector<double> error={ sigmaSlope2, sigmaSlopeIntercept, sigmaIntercept2 };

    return FitResult2D(slope1, intersept1, error , interceptz);

}

FitResult2D regressionYZ(const HoughTransformer::HitCluster& cluster, double interceptz=0) {
    double sumY = 0;
    double sumZ = 0;
    double sumYZ = 0;
    double sumZsquare = 0;  // = Sum (Z^2)
    double sumW = 0;

    for(auto& h : cluster) {
    	double errory2=1;//add error!
    	double hiz=h.z-interceptz;
		sumY += h.y/errory2;
		sumZ += hiz/errory2;
		sumYZ += h.y*hiz/errory2;
		sumZsquare += hiz*hiz/errory2;
		sumW+=1/errory2;
    }

    // z <-> y from LinearTrackRegression::doRegressionX of marlinTPC
    /** gives back parameters for \n
     *  x = slope1 * z + intersept1\n
     *  y = slope2 * z + intersept2\n
     *
     */

    double denominator=(sumZ * sumZ - sumW * sumZsquare);
    if(std::fabs(denominator)<1E-20){
    	std::cerr<<"error: (sumZ * sumZ - ntot * sumZsquare)<1E-20"<<std::endl;
    	std::cerr<<"sumZ="<<sumZ<<" sumZsquare="<<sumZsquare<<" ntot="<<sumW<<std::endl;
    	std::cerr<<cluster.size()<<" hits on "<<cluster.getNPlanesHit()<<" planes"<<std::endl;
    	throw "(sumZ * sumZ - ntot * sumZsquare)<1E-20";
    }

    double slope1     = (sumY * sumZ - sumW * sumYZ) / denominator;
    double intersept1 = (sumZ * sumYZ - sumZsquare * sumY) / denominator;

    double sigmaSlope2=sumZsquare/denominator;
    double sigmaIntercept2=sumW/denominator;
    double sigmaSlopeIntercept=-sumZ/denominator;
    std::vector<double> error={ sigmaSlope2, sigmaSlopeIntercept, sigmaIntercept2 };

    return FitResult2D(slope1, intersept1, error, interceptz );

}


FitResult3D regressionFit3d(const HoughTransformer::HitCluster& cluster, double interceptz) {
	return FitResult3D {
		regressionXZ(cluster, interceptz),
		regressionYZ(cluster, interceptz)
	};
}


Residual calculateResidual( const PositionHit& h, const FitResult3D& fit ) {
	return Residual{
		h.x - fit.xAt(h.z),
		h.y - fit.yAt(h.z),
		0,
		h};
}
std::vector<Residual> calculateResiduals( const HoughTransformer::HitCluster& cluster, const FitResult3D& fit)  {
	std::vector<Residual> residuals;
	for(auto& h : cluster) {
		residuals.emplace_back( calculateResidual(h, fit) );
	}
	return residuals;
}

//function invalidates residuals reference to hits!
HoughTransformer::HitCluster&  cutOnResiduals( HoughTransformer::HitCluster& cluster, const std::vector<Residual>& residuals, double maxResidual ) {
	auto res=residuals.begin();
	int nremoved=cluster.size();
	cluster.remove_if( [&res, &maxResidual](const PositionHit&){
		bool cut= res->x>maxResidual || res->y>maxResidual;
		++res;
		return cut;
	} );
	nremoved-=cluster.size();
//	std::cout<<"Removed "<<nremoved<<" hits"<<std::endl;
	return cluster;
}
TVector3 averageResidual(const std::vector<Residual>& residuals) {
	double x=0, y=0, z=0;
	for(auto& r: residuals) {
		x+=r.x;
		y+=r.y;
		z+=r.z;
	}
	return { x/=residuals.size(),
			 y/=residuals.size(),
			 z/=residuals.size() };
}
std::ostream& operator <<(std::ostream& os, FitResult3D& fit) {
	return os<<"SimpleFitResult: slopes("<<fit.XZ.slope<<", "<<fit.YZ.slope<<") intercepts("<<fit.XZ.intercept<<", "<<fit.YZ.intercept<<")";
}
