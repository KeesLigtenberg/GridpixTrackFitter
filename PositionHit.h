#ifndef POSITIONHIT_H
#define POSITIONHIT_H

#include "Hit.h"
#include "TVector3.h"

struct PositionHit {
	PositionHit(double x, double y, double z, int plane=0, int ToT=1) : x(x), y(y), z(z), plane(plane), ToT(ToT) {};
	double x,y,z;
	int plane;
	int ToT=1;
	TVector3 getPosition() const { return {x,y,z}; };
	void SetPosition(const TVector3& v) { x=v.x(); y=v.y(); z=v.z(); };
};

inline std::vector<PositionHit> convertHits(const std::vector<Hit>& hv, double planePosition, double pixelwidth, double pixelheight, int plane=0) {
		std::vector<PositionHit> phv;
		for(auto& h : hv) {
			// .5 to center position in pixel
			phv.emplace_back( PositionHit{(h.column+.5)*pixelwidth, (h.row+.5)*pixelheight, planePosition, plane} );
		}
		return phv;
}

template< class container=std::vector<PositionHit> >
std::vector<PositionHit> convertHits(const std::vector<TimePixHit>& hv, double pixelwidth, double pixelheight, double driftSpeed, int plane=0 ) {
		container phv;
		for(auto& h : hv) {
			// .5 to center position in pixel
			int plane=0;
			phv.emplace_back( PositionHit{h.driftTime*driftSpeed /*z*/, (h.row+.5)*pixelheight /*y*/, (h.column+.5)*pixelwidth/*x*/, plane, h.charge /*charge=ToT*/} );
		}
		return phv;
}

#endif
