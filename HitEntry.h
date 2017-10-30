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
	HitEntry() : rx(), ry(), rz(),x(),y(),z(), ToT(), row(), col() {};
	HitEntry(const Residual& r) :
		rx(r.x),
		ry(r.y),
		rz(r.z),
		x(r.h.x),
		y(r.h.y),
		z(r.h.z),
		ToT(r.h.ToT),
		row(r.h.row),
		col(r.h.column)
	{};
	double rx, ry, rz;
	double x, y, z;
	int ToT;
	int row, col;
};
#pragma link C++ class std::vector<HitEntry>+;
#pragma link C++ class std::vector< std::vector<HitEntry> >+;


#endif /* HITENTRY_H_ */
