/*
 * detectorConfiguration.h
 *
 *  Created on: Jun 23, 2017
 *      Author: cligtenb
 */

#ifndef DETECTORCONFIGURATION_H_
#define DETECTORCONFIGURATION_H_

#include <utility>

struct DetectorConfiguration {
	virtual ~DetectorConfiguration() {}
	DetectorConfiguration (int nPlanes, std::vector<double> planePosition, double pixelsize, int pixelColumns, int pixelRows) :
			nPlanes(nPlanes),
			planePosition(planePosition),
			pixelsize(pixelsize),
			pixelColumns(pixelColumns),
			pixelRows(pixelRows)
		{};

	int nPlanes;
	std::vector<double> planePosition; //mm
	double pixelsize;//mm
	int pixelColumns, pixelRows;

	double planexmax() const { return pixelColumns*pixelsize; }
	double planeymax() const { return pixelRows*pixelsize; }

	virtual double xmax() const { return planexmax(); }
	virtual double xmin() const { return 0.; }
	virtual double ymax() const { return planeymax(); }
	virtual double ymin() const { return 0.; }
	virtual double zmax() const { return planePosition.back(); }
	virtual double zmin() const { return planePosition.front(); }
	virtual double zmean() const { return (zmin()+zmax())/2; }

	virtual std::pair<double,double> getCentre() const { return { (xmax()-xmin())/2, (ymax()-ymin())/2 }; };

};



#endif /* DETECTORCONFIGURATION_H_ */
