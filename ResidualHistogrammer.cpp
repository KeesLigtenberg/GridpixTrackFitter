/*
 * ResidualHistogrammer.cpp
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#include "ResidualHistogrammer.h"

#include <iostream>

#include "TF1.h"
#include "TFitResult.h"

int ResidualHistogrammer::PlaneHistograms::n=0;

ResidualHistogrammer::ResidualHistogrammer(std::string outputFileName, const DetectorConfiguration& detector) :
		outputFile(outputFileName.c_str(), "RECREATE"),
		detector(detector),
		planeHist(),
		xResidualByToT( "xResidualByToT", "xResidualByToT", 40, 0,1)
{
	planeHist.reserve(detector.nPlanes);
	for(int i=0; i<detector.nPlanes; ++i) planeHist.emplace_back( detector );
	ResidualHistogrammer::PlaneHistograms::n=0;
}

ResidualHistogrammer::~ResidualHistogrammer() {
	outputFile.Write();
	outputFile.Close();
}

void ResidualHistogrammer::fill(const Residual& r, const std::pair<double, double>& rotationPoint) {
	int planeNumber=r.h.plane;
	auto& plane= planeHist[planeNumber];
	auto& h=r.h;
	plane.xResidual.Fill( r.x );
	plane.yResidual.Fill( r.y );
	plane.xResidualByPixel.Fill( h.x, h.y, r.x );
	plane.yResidualByPixel.Fill( h.x, h.y, r.y );

	xResidualByToT.Fill(h.ToT/4096.*25, r.x);

	double xc=rotationPoint.first, yc=rotationPoint.second; //x and y center
	double hx=h.x-xc, hy=h.y-yc;
	double phi=(hy*r.x-hx*r.y)/(hx*hx+hy*hy);
	double weight=hx*hx+hy*hy;
	plane.zRotation.Fill( phi , weight);
	plane.zRotationByPixel.Fill(h.x, h.y, phi);
}

void ResidualHistogrammer::fill(const std::vector<Residual>& residuals, const std::vector<std::pair<double, double> >& rotationPoint) {
	if(!residuals.size()) return;
	for(auto& r:residuals) fill(r, rotationPoint[r.h.plane]);
}
void ResidualHistogrammer::fill(const std::vector<Residual>& residuals, const std::pair<double, double>& rotationPoint) {
	if(!residuals.size()) return;
	for(auto& r:residuals) fill(r, rotationPoint);
}

//get x,y mean
std::vector<std::pair<double, double> > ResidualHistogrammer::getMeansOfPlanes() {
	std::vector<std::pair<double, double>> vec;
	for(auto& p : planeHist) {
		vec.push_back(p.getMeansFromFit());
//		vec.push_back( {p.xResidual.GetMean(), p.yResidual.GetMean()} );
	}
	return vec;
}

std::pair<double, double> ResidualHistogrammer::PlaneHistograms::getMeansFromFit() {
	return std::make_pair( getMeanFromGausFit(xResidual), getMeanFromGausFit(yResidual) );
}

std::vector<double> ResidualHistogrammer::getRotationOfPlanes() {
	std::vector<double> vec;
	for(auto& p: planeHist) {
		vec.push_back( p.getRotationFromFit() );
//		vec.push_back( p.zRotation.GetMean() ); //get from mean in histogram
	}
	return vec;
}

double ResidualHistogrammer::PlaneHistograms::getRotationFromFit() {
	return getMeanFromGausFit(zRotation);
}

TrackHistogrammer::TrackHistogrammer(const DetectorConfiguration& detector) :
//	phi("trackPhi", "track #Phi; #phi [rad.]; tracks", 20,M_PI/2.-0.01,M_PI/2.+0.01),
//	d0("trackd0", "track d_{0}; d_{0} [mm]; tracks", 20, 0, detector.planexmax() ),
//	tanLambda("trackTanLambda", "track tan(#lambda); tan(#lambda); tracks", 20,-0.01,0.01),
//	z0("trackz0", "track z_{0}; z_{0} [mm]; tracks", 20, 0, detector.planeymax() ),
	slope1("slope1", "slope 1 (X); slope; tracks", 20, -1, 1),
	slope2("slope2", "slope 2 (Y); slope; tracks", 20, -1, 1),
	intercept1("intercept1", "intercept 1 (X); intercept; tracks", 20, -0.01, 0.01),
	intercept2("intercept2", "intercept 2 (Y); intercept; tracks", 20, -0.01, 0.01),
	detector(detector)
		{};

void TrackHistogrammer::fill(const FitResult3D& fit) {
//	TrackFitResult entry(fit);
//	phi.Fill(entry.phi);
//	d0.Fill(entry.d0);
//	tanLambda.Fill(entry.tanlambda);
//	z0.Fill(entry.z0);

	slope1.Fill(fit.XZ.slope);
	slope2.Fill(fit.YZ.slope);
	intercept1.Fill(fit.XZ.intercept);
	intercept2.Fill(fit.YZ.intercept);
}
