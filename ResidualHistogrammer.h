/*
 * ResidualHistogrammer.h
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#ifndef RESIDUALHISTOGRAMMER_H_
#define RESIDUALHISTOGRAMMER_H_
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TProfile2D.h"

#include "linearRegressionFit.h"
#include "DetectorConfiguration.h"

class ResidualHistogrammer {
public:
	ResidualHistogrammer(std::string outputFileName, const DetectorConfiguration& detector);
	virtual ~ResidualHistogrammer();

	void fill(const Residual&);
	void fill(const std::vector<Residual>&);


	std::vector< std::pair<double, double> > getMeansOfPlanes();
	std::vector< double > getRotationOfPlanes();

	const DetectorConfiguration& detector;

	TFile outputFile;

	struct PlaneHistograms {
		static int n;
		PlaneHistograms() :
			xResidual( ("xResidual_"+std::to_string(n)).c_str(), ("xResidual_"+std::to_string(n)).c_str(), 20, -0.2, 0.2),
			yResidual( ("yResidual_"+std::to_string(n)).c_str(), ("yResidual_"+std::to_string(n)).c_str(), 20, -0.2, 0.2),
			xResidualByPixel( ("xResidualByPixel_"+std::to_string(n)).c_str(), ("xResidualByPixel_"+std::to_string(n)).c_str(), 20, 0, 1153*0.0184, 20, 0, 577*0.0184 ),
			yResidualByPixel( ("yResidualByPixel_"+std::to_string(n)).c_str(), ("yResidualByPixel_"+std::to_string(n)).c_str(), 20, 0, 1153*0.0184, 20, 0, 577*0.0184 ),
			rotation( ("rotation_"+std::to_string(n)).c_str(), ("rotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05)
		{ ++n; }
		std::pair<double,double> getMeansFromFit();
		double getRotationFromFit();
		TH1D xResidual, yResidual;
		TProfile2D xResidualByPixel, yResidualByPixel;
		TH1D rotation;
	} planeHist[6];
};


#endif /* RESIDUALHISTOGRAMMER_H_ */
