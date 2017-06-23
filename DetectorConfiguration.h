/*
 * detectorConfiguration.h
 *
 *  Created on: Jun 23, 2017
 *      Author: cligtenb
 */

#ifndef DETECTORCONFIGURATION_H_
#define DETECTORCONFIGURATION_H_

struct DetectorConfiguration {

	int nPlanes;
	std::vector<double> planePosition;
	double pixelsize;//mm
	int pixelColumns, pixelRows;

	double planexmax() const { return pixelColumns*pixelsize; }
	double planeymax() const { return pixelRows*pixelsize; }

};



#endif /* DETECTORCONFIGURATION_H_ */
