/*
 * goodArea.h
 *
 *  Created on: Nov 27, 2017
 *      Author: cligtenb
 */

#ifndef GOODAREA_H_
#define GOODAREA_H_

#include "TBox.h"

namespace goodArea {
	constexpr int xmin=20, xmax=236, ybl=80, ybr=36, yul=232, yur=180;
	//create box
	const TBox excluded1(100,116,116,140),
			excluded2(80,56,88,136);
	bool isInsideArea(int x, int y) {
		//area:|_--_
		//	     --_|
		if(x<=xmin or x>=xmax) {return false;}
		if( (y<= ybl-(ybl-ybr)/double(xmax-xmin)*(x-xmin))
		 or (y>= yul-(yul-yur)/double(xmax-xmin)*(x-xmin)) ) { return false; }
		if(excluded1.IsInside(x,y)) { return false; }
//		if(excluded2.IsInside(x,y)) { return false; }
		return true;
	}
}

#endif /* GOODAREA_H_ */
