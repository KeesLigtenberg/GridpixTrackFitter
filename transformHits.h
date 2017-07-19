/*
 * transformHits.h
 *
 *  Created on: Jun 16, 2017
 *      Author: cligtenb
 */

#ifndef TRANSFORMHITS_H_
#define TRANSFORMHITS_H_

#include "PositionHit.h"

inline PositionHit& translateHit(PositionHit& h, double dx, double dy ) {
	h.x-=dx;
	h.y-=dy;
	return h;
}
inline PositionHit& translateHit(PositionHit& h, const std::pair<double,double>& translation) { return translateHit(h,translation.first, translation.second); }

inline PositionHit& rotateHit(PositionHit& h, double rotation, const std::pair<double, double>& rotationPoint) {
	double sinr=sin(rotation), cosr=cos(rotation);

	const double xc=rotationPoint.first, yc=rotationPoint.second; //x and y center

	h.x=cosr*(h.x-xc)-sinr*(h.y-yc)+xc;
	h.y=cosr*(h.y-yc)+sinr*(h.x-xc)+yc;
	return h;
}

template<class hitCollection>
hitCollection& translateHits( hitCollection& hits, const std::vector<std::pair<double,double> >& translation ) {
	for(PositionHit& h : hits) translateHit(h, translation[h.plane] );
	return hits;
};

template<class hitCollection>
hitCollection& translateHits( hitCollection& hits, const std::pair<double,double>& translation ) {
	for(PositionHit& h : hits) translateHit(h, translation );
	return hits;
};

template<>
std::vector<std::vector<PositionHit> >& translateHits( std::vector<std::vector<PositionHit> >& hits, const std::vector<std::pair<double,double> >& translation ) {
	for(auto& v: hits)
		for(PositionHit& h : v)
			h=translateHit(h, translation[h.plane] );
	return hits;
};

template<class hitCollection>
hitCollection& rotateHits(hitCollection& hits, const std::vector<double>& rotations, const std::vector<std::pair<double, double>>& rotationPoints ) {
	for(PositionHit& h : hits) {
		h=rotateHit(h,rotations[h.plane], rotationPoints[h.plane]);
	}
	return hits;
}

template<class hitCollection>
hitCollection& rotateHits(hitCollection& hits, double rotations, const std::pair<double, double>& rotationPoints ) {
	for(PositionHit& h : hits) {
		h=rotateHit(h,rotations, rotationPoints);
	}
	return hits;
}

template<>
std::vector<std::vector<PositionHit> >& rotateHits(std::vector<std::vector<PositionHit> >& hits, const std::vector<double>& rotations, const std::vector<std::pair<double, double>>& rotationPoints) {
	for(auto& v : hits) {
		v=rotateHits(v,rotations,rotationPoints);
	}
	return hits;
}

#endif /* TRANSFORMHITS_H_ */
