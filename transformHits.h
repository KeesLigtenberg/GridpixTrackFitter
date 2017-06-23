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

inline PositionHit& rotateHit(PositionHit& h, double rotation) {
	double sinr=sin(rotation), cosr=cos(rotation);

	const double xc=1153*0.0184/2, yc= 577*0.0184/2; //x and y center

	h.x=cosr*(h.x-xc)-sinr*(h.y-yc)+xc;
	h.y=cosr*(h.y-yc)+sinr*(h.x-xc)+yc;
	return h;
}

template<class hitCollection>
hitCollection& translateHits( hitCollection& hits, const std::vector<std::pair<double,double> >& translation ) {
	for(PositionHit& h : hits) translateHit(h, translation[h.plane] );
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
hitCollection& rotateHits(hitCollection& hits, const std::vector<double>& rotations) {
	for(PositionHit& h : hits) {
		h=rotateHit(h,rotations[h.plane]);
	}
	return hits;
}

template<>
std::vector<std::vector<PositionHit> >& rotateHits(std::vector<std::vector<PositionHit> >& hits, const std::vector<double>& rotations) {
	for(auto& v : hits) {
		v=rotateHits(v,rotations);
	}
	return hits;
}

#endif /* TRANSFORMHITS_H_ */
