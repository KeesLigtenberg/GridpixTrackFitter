/*
 * linearRegressionFit.cpp
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#include <iostream>

#include "linearRegressionFit.h"

SimpleFitResult linearRegressionFit(const HoughTransformer::HitCluster& cluster ) {

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
     *  what kind of regression is this? should we do perpendicular regression?
     */
    double slope1     = (sumX * sumZ - ntot * sumXZ) / (sumZ * sumZ - ntot * sumZsquare);
    double intersept1 = (sumZ * sumXZ - sumZsquare * sumX) / (sumZ * sumZ - ntot * sumZsquare);
    double slope2     = (sumY * sumZ - ntot * sumYZ) / (sumZ * sumZ - ntot * sumZsquare);
    double intersept2 = (sumZ * sumYZ - sumZsquare * sumY) / (sumZ * sumZ - ntot * sumZsquare);

    return SimpleFitResult {slope1, intersept1, slope2, intersept2};
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


