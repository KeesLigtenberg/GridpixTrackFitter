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
#include "TProfile.h"
#include "TProfile2D.h"

#include "linearRegressionFit.h"
#include "DetectorConfiguration.h"

class ResidualHistogrammer {
public:
	ResidualHistogrammer(std::string outputFileName, const DetectorConfiguration& detector);
	virtual ~ResidualHistogrammer();

	void fill(const Residual&, const std::pair<double, double>& rotationPoint);
	void fill(const std::vector<Residual>&, const std::vector<std::pair<double, double> >& rotationPoints );
	void fill(const std::vector<Residual>& residuals, const std::pair<double, double>& rotationPoint);
	void fill(const std::vector<Residual>& residualVector) {
		std::vector< std::pair<double,double> > rotationPoints(detector.nPlanes, detector.getCentre());
		fill(residualVector, rotationPoints );
	}


	std::vector< std::pair<double, double> > getMeansOfPlanes();
	std::vector< double > getRotationOfPlanes();

	TFile outputFile;
	const DetectorConfiguration& detector;

	struct PlaneHistograms {
		static int n;
		PlaneHistograms(const DetectorConfiguration& det) :
			xResidual( ("xResidual_"+std::to_string(n)).c_str(), ("xResidual_"+std::to_string(n)).c_str(), 40, -0.1, 0.1),
			yResidual( ("yResidual_"+std::to_string(n)).c_str(), ("yResidual_"+std::to_string(n)).c_str(), 40, -0.1, 0.1),
			zRotation( ("zRotation_"+std::to_string(n)).c_str(), ("zRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
			xResidualByPixel( ("xResidualByPixel_"+std::to_string(n)).c_str(), ("xResidualByPixel_"+std::to_string(n)).c_str(), 20, det.xmin(), det.xmax(), 20, 0, det.ymax() ),
			yResidualByPixel( ("yResidualByPixel_"+std::to_string(n)).c_str(), ("yResidualByPixel_"+std::to_string(n)).c_str(), 20, det.xmin(), det.xmax(), 20, 0, det.ymax() ),
//			xRotation( ("xRotation_"+std::to_string(n)).c_str(), ("xRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
//			yRotation( ("yRotation_"+std::to_string(n)).c_str(), ("yRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
			zRotationByPixel( ("zRotationByPixel_"+std::to_string(n)).c_str(), ("zRotationByPixel_"+std::to_string(n)).c_str(), 20, det.xmin(), det.xmax(), 20, 0, det.ymax() )
		{ ++n; }
		std::pair<double,double> getMeansFromFit();
		double getRotationFromFit();
		TH1D xResidual, yResidual;
		TH1D zRotation;// xRotation, yRotation;
		TProfile2D xResidualByPixel, yResidualByPixel, zRotationByPixel;
	};
	std::vector<PlaneHistograms> planeHist;
	TProfile xResidualByToT;

};


class TrackHistogrammer {
public:
	TrackHistogrammer(const DetectorConfiguration& detector);
	void fill(SimpleFitResult entry);
private:
	TH1D slope1, slope2, intercept1, intercept2;
	TH1D phi, d0, tanLambda, z0;
	const DetectorConfiguration& detector;
};


#endif /* RESIDUALHISTOGRAMMER_H_ */
