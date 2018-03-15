/*
 * testBeamSetup.h
 *
 *  Created on: Sep 28, 2017
 *      Author: cligtenb
 */

#ifndef TESTBEAMSETUP_H_
#define TESTBEAMSETUP_H_


struct TimePixDetectorConfiguration : DetectorConfiguration {
	constexpr static double driftSpeed=0.075; //mm/ns
	TimePixDetectorConfiguration() : DetectorConfiguration{
		1, {0}, //nplanes, planeposition
		0.055, 256, 256 //pixelsize, xpixels, ypixels
	} {};
	virtual double xmin() const {return 0.*driftSpeed; }
	virtual double xmax() const {return 300*driftSpeed; }
	virtual double zmin() const {return 0; }
	virtual double zmax() const {return pixelColumns*pixelsize; };
} timePixChip;

const DetectorConfiguration mimosa= {
	6, //planes
//	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	{0, 15.8, 31.8, 143.1, 161.55, 179.91 }, //plane positions as measured
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1, (our data does start at 0) //TODO: remove one extra when done with exampledata.
};

const DetectorConfiguration combinedSetupForDrawing {
		2, {-350,200 }, //nplanes, planeposition
		0.001, int(1000*mimosa.xmax()), int(1000*mimosa.ymax()) //pixelsize, xpixels, ypixels
	};


#endif /* TESTBEAMSETUP_H_ */
