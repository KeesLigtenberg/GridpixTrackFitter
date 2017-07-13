/*
 * detectorConfiguration.h
 *
 *  Created on: Jun 23, 2017
 *      Author: cligtenb
 */

#ifndef DETECTORCONFIGURATION_H_
#define DETECTORCONFIGURATION_H_

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

	virtual double xmax() const { return planexmax(); };
	virtual double xmin() const { return 0.; };
	virtual double ymax() const { return planeymax(); };
	virtual double ymin() const { return 0.; };

};



#endif /* DETECTORCONFIGURATION_H_ */
