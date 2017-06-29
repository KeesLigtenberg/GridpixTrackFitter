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

	TFile outputFile;
	const DetectorConfiguration& detector;

	struct PlaneHistograms {
		static int n;
		PlaneHistograms() :
			xResidual( ("xResidual_"+std::to_string(n)).c_str(), ("xResidual_"+std::to_string(n)).c_str(), 40, -0.2, 0.2),
			yResidual( ("yResidual_"+std::to_string(n)).c_str(), ("yResidual_"+std::to_string(n)).c_str(), 40, -0.2, 0.2),
			xResidualByPixel( ("xResidualByPixel_"+std::to_string(n)).c_str(), ("xResidualByPixel_"+std::to_string(n)).c_str(), 20, 0, 1153*0.0184, 20, 0, 577*0.0184 ),
			yResidualByPixel( ("yResidualByPixel_"+std::to_string(n)).c_str(), ("yResidualByPixel_"+std::to_string(n)).c_str(), 20, 0, 1153*0.0184, 20, 0, 577*0.0184 ),
//			xRotation( ("xRotation_"+std::to_string(n)).c_str(), ("xRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
//			yRotation( ("yRotation_"+std::to_string(n)).c_str(), ("yRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
			zRotation( ("zRotation_"+std::to_string(n)).c_str(), ("zRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05)
		{ ++n; }
		std::pair<double,double> getMeansFromFit();
		double getRotationFromFit();
		TH1D xResidual, yResidual;
		TProfile2D xResidualByPixel, yResidualByPixel;
		TH1D zRotation;// xRotation, yRotation;
	} planeHist[6];
};


class TrackHistogrammer {
public:
	TrackHistogrammer(const DetectorConfiguration& detector);
	void fill(TrackFitResult entry);
private:
	TH1D phi, d0, tanLambda, z0;
	const DetectorConfiguration& detector;
};


#endif /* RESIDUALHISTOGRAMMER_H_ */
