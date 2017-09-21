/*
 * linearRegressionFit.cpp
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#include <iostream>

#include "linearRegressionFit.h"


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
    std::array<double,3> error={ sigmaSlope2, sigmaSlopeIntercept, sigmaIntercept2 };

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
    std::array<double,3> error={ sigmaSlope2, sigmaSlopeIntercept, sigmaIntercept2 };

    return FitResult2D(slope1, intersept1, error, interceptz );

}

FitResult3D regression3d(const HoughTransformer::HitCluster& cluster, double interceptz=0) {
	return FitResult3D {
		regressionXZ(cluster, interceptz),
		regressionYZ(cluster, interceptz)
	};
}

SimpleFitResult linearRegressionFit(const HoughTransformer::HitCluster& cluster) {

    double sumX = 0;
    double sumZ = 0;
    double sumY = 0;
    double sumXZ = 0;
    double sumYZ = 0;
    double sumZsquare = 0;  // = Sum (Z^2)
    int ntot = 0;

    for(auto& h : cluster) {
		sumX += h.x;
		sumY += h.y;
		sumZ += h.z;
		sumXZ += h.x*h.z;
		sumYZ += h.y*h.z;
		sumZsquare += h.z*h.z;
		ntot++;
    }

    // z <-> y from LinearTrackRegression::doRegressionX of marlinTPC
    /** gives back parameters for \n
     *  x = slope1 * z + intersept1\n
     *  y = slope2 * z + intersept2\n
     *
     */

    double denominator=(sumZ * sumZ - ntot * sumZsquare);
    if(std::fabs(denominator)<1E-20){
    	std::cerr<<"error: (sumZ * sumZ - ntot * sumZsquare)<1E-20"<<std::endl;
    	std::cerr<<"sumZ="<<sumZ<<" sumZsquare="<<sumZsquare<<" ntot="<<ntot<<std::endl;
    	std::cerr<<cluster.size()<<" hits on "<<cluster.getNPlanesHit()<<" planes"<<std::endl;
    	throw "(sumZ * sumZ - ntot * sumZsquare)<1E-20";
    }

    double slope1     = (sumX * sumZ - ntot * sumXZ) / denominator;
    double intersept1 = (sumZ * sumXZ - sumZsquare * sumX) / denominator;
    double slope2     = (sumY * sumZ - ntot * sumYZ) / denominator;
    double intersept2 = (sumZ * sumYZ - sumZsquare * sumY) / denominator;


	double dslope1 = 1 / sqrt( sumZsquare - sumZ*sumZ/ntot );
	double dintersept1 = dslope1 * sqrt( sumZsquare/ntot );
	double dslope2 = 1 / sqrt( sumZsquare - sumZ*sumZ/ntot );
	double dintersept2 = dslope1 * sqrt( sumZsquare/ntot );


    return SimpleFitResult {
    	slope1, intersept1, slope2, intersept2,
    	dslope1, dintersept1, dslope2, dintersept2};
}



TrackFitResult linearRegressionTrackFit(const HoughTransformer::HitCluster& cluster) {
	return TrackFitResult( linearRegressionFit(cluster) );
}

Residual calculateResidual( const PositionHit& h, const SimpleFitResult& fit ) {
	return Residual{
		h.x - fit.slope1 * h.z - fit.intersept1,
		h.y - fit.slope2 * h.z - fit.intersept2,
		0,
		h};
}
std::vector<Residual> calculateResiduals( const HoughTransformer::HitCluster& cluster, const SimpleFitResult& fit)  {
	std::vector<Residual> residuals;
	for(auto& h : cluster) {
		residuals.emplace_back( calculateResidual(h, fit) );
	}
	return residuals;
}
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
std::ostream& operator <<(std::ostream& os, SimpleFitResult& fit) {
	return os<<"SimpleFitResult: slopes("<<fit.slope1<<", "<<fit.slope2<<") intersepts("<<fit.intersept1<<", "<<fit.intersept2<<")";
}
