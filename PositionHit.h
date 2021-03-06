#ifndef POSITIONHIT_H
#define POSITIONHIT_H

#include "Hit.h"
#include "TVector3.h"

TVector3& RotateAroundPoint(TVector3& v, double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis) {
	v-=rotationPoint;
	v.Rotate(rotation, rotationAxis);
	v+=rotationPoint;
	return v;
};

struct PositionHit {
	PositionHit() : PositionHit(0,0,0) {};
	PositionHit(double x, double y, double z, int plane=0, int row=0, int column=0, int ToT=1, double error2x=1, double error2y=1) :
			x(x), y(y), z(z), row(row), column(column), plane(plane), ToT(ToT), error2x(error2x), error2y(error2y) {};
	double x,y,z;
	int row, column, plane;
	double ToT=1;
	double error2x=1, error2y=1;
	int flag=1;

	TVector3 getPosition() const { return {x,y,z}; };
	void SetPosition(const TVector3& v) { x=v.x(); y=v.y(); z=v.z(); };
	void RotatePosition(double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis) {
		TVector3 v=getPosition();
		v=RotateAroundPoint(v, rotation, rotationPoint, rotationAxis);
		SetPosition(v);
	};
};

struct Residual {
	double x,y,z;
	const PositionHit& h;
	TVector3 getVector() { return TVector3(x,y,z); };
	void setVector(const TVector3& v) { x=v.x(); y=v.y(); z=v.z(); };
};


inline std::vector<PositionHit> convertHits(const std::vector<Hit>& hv, double planePosition, double pixelwidth, double pixelheight, int plane=0) {
		std::vector<PositionHit> phv;
		for(auto& h : hv) {
			// .5 to center position in pixel
			int ToTTelescope=1;
			phv.emplace_back(
					PositionHit{
						(h.column+.5)*pixelwidth, (h.row+.5)*pixelheight, planePosition,
								plane, h.row, h.column,
								ToTTelescope, pixelwidth*pixelwidth/12, pixelheight*pixelheight/12
					} );
		}
		return phv;
}

template< class Container=std::vector<PositionHit> >
Container convertHits(const std::vector<TimePixHit>& hv, double pixelwidth, double pixelheight, double driftSpeed, int plane=0 ) {
		Container phv;
		for(auto& h : hv) {
			// .5 to center position in pixel
			int plane=0;
			phv.emplace_back( PositionHit{
				h.driftTime*driftSpeed /*z*/, (h.row+.5)*pixelheight /*y*/, (h.column+.5)*pixelwidth/*x*/,
						plane, h.row, h.column,
						h.charge /*charge=ToT*/} );
		}
		return phv;
}

#endif
