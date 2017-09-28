/*
 * testBeamSetup.h
 *
 *  Created on: Sep 28, 2017
 *      Author: cligtenb
 */

#ifndef TESTBEAMSETUP_H_
#define TESTBEAMSETUP_H_


struct TimePixDetectorConfiguration : DetectorConfiguration {
	static constexpr double driftSpeed=0.075; //mm/ns
	TimePixDetectorConfiguration() : DetectorConfiguration{
		1, {0}, //nplanes, planeposition
		0.055, 256, 256 //pixelsize, xpixels, ypixels
	} {};
	virtual double xmin() const {return driftSpeed; }
	virtual double xmax() const {return 225*driftSpeed; }
	virtual double zmin() const {return 0; }
	virtual double zmax() const {return pixelColumns*pixelsize; };
} timePixChip;

const DetectorConfiguration mimosa= {
	6, //planes
//	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	{0, 15.8, 31.8, 143.1, 161.55, 179.91 }, //plane positions as measured
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other!
};


#endif /* TESTBEAMSETUP_H_ */