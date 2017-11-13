/*
 * HitEntry.h
 *
 *  Created on: Sep 25, 2017
 *      Author: cligtenb
 */

#ifndef HITENTRY_H_
#define HITENTRY_H_

#include "PositionHit.h"

struct HitEntry {
	HitEntry() : rx(), ry(), rz(),x(),y(),z(), ToT(), row(), col(), flag() {};
	HitEntry(const Residual& r) :
		rx(r.x),
		ry(r.y),
		rz(r.z),
		x(r.h.x),
		y(r.h.y),
		z(r.h.z),
		ToT(r.h.ToT),
		row(r.h.row),
		col(r.h.column),
		flag(r.h.flag)
	{};
	double rx, ry, rz; //residuals in timepix frame
	double x, y, z; //position in telescope frame
	int ToT; //ToT in [0.025 ns]
	int row, col; //row col of hit
	int flag;

	PositionHit createPositionHit() {
		return PositionHit(x,y,z, 0,row,col,ToT,1,1);
	}
};
#pragma link C++ class std::vector<HitEntry>+;
#pragma link C++ class std::vector< std::vector<HitEntry> >+;


#endif /* HITENTRY_H_ */
