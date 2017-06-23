#ifndef POSITIONHIT_H
#define POSITIONHIT_H

#include "Hit.h"

struct PositionHit {
	PositionHit(double x, double y, double z, int plane=0) : x(x), y(y), z(z), plane(plane) {};
	double x,y,z;
	int plane;
};

inline std::vector<PositionHit> convertHits(const std::vector<Hit>& hv, double planePosition, double pixelwidth, double pixelheight, int plane=0) {
		std::vector<PositionHit> phv;
		for(auto& h : hv) {
			// .5 to center position in pixel
			phv.emplace_back( PositionHit{(h.column+.5)*pixelwidth, (h.row+.5)*pixelheight, planePosition, plane} );
		}
		return phv;
}


#endif
