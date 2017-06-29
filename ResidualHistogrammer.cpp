/*
 * ResidualHistogrammer.cpp
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#include "ResidualHistogrammer.h"

#include <iostream>

#include "TFitResult.h"

int ResidualHistogrammer::PlaneHistograms::n=0;

ResidualHistogrammer::ResidualHistogrammer(std::string outputFileName, const DetectorConfiguration& detector) :
	outputFile(outputFileName.c_str(), "RECREATE"),
	detector(detector)
{
	ResidualHistogrammer::PlaneHistograms::n=0;
}

ResidualHistogrammer::~ResidualHistogrammer() {
	outputFile.Write();
	outputFile.Close();
}

void ResidualHistogrammer::fill(const Residual& r) {
	int planeNumber=r.h.plane;
	auto& plane= planeHist[planeNumber];
	auto& h=r.h;
	plane.xResidual.Fill( r.x );
	plane.yResidual.Fill( r.y );
	plane.xResidualByPixel.Fill( h.x, h.y, r.x );
	plane.yResidualByPixel.Fill( h.x, h.y, r.y );

	double xc=detector.planexmax()/2, yc= detector.planeymax()/2; //x and y center
	double hx=h.x-xc, hy=h.y-yc;
	plane.zRotation.Fill( (hy*r.x-hx*r.y)/(hx*hx+hy*hy) );
}

void ResidualHistogrammer::fill(const std::vector<Residual>& residuals) {
	if(!residuals.size()) return;
	for(auto& r:residuals) fill(r);
}

//get x,y mean
std::vector<std::pair<double, double> > ResidualHistogrammer::getMeansOfPlanes() {
	std::vector<std::pair<double, double>> vec;
	for(auto& p : planeHist) {
		vec.push_back(p.getMeansFromFit());
	}
	return vec;
}

double getMeanFromGausFit( TH1& hist ) {
	TFitResultPtr fitresult=hist.Fit("gaus", "SQ"); //Store and Quiet
	if(!fitresult->IsValid()) {std::cerr<< "getMeanFromGausFit error: failed to fit histogram "<<hist.GetName()<<std::endl; return 0;}
	return fitresult->Parameter(1);
}

std::pair<double, double> ResidualHistogrammer::PlaneHistograms::getMeansFromFit() {
	return std::make_pair( getMeanFromGausFit(xResidual), getMeanFromGausFit(yResidual) );
}

std::vector<double> ResidualHistogrammer::getRotationOfPlanes() {
	std::vector<double> vec;
	for(auto& p: planeHist) {
		vec.push_back( p.getRotationFromFit() );
	}
	return vec;
}

double ResidualHistogrammer::PlaneHistograms::getRotationFromFit() {
	return getMeanFromGausFit(zRotation);
}

TrackHistogrammer::TrackHistogrammer(const DetectorConfiguration& detector) :
	phi("trackPhi", "track #Phi; #phi [rad.]; tracks", 20,M_PI/2.-0.01,M_PI/2.+0.01),
	d0("trackd0", "track d_{0}; d_{0} [mm]; tracks", 20, 0, detector.planexmax() ),
	tanLambda("trackTanLambda", "track tan(#lambda); tan(#lambda); tracks", 20,-0.01,0.01),
	z0("trackz0", "track z_{0}; z_{0} [mm]; tracks", 20, 0, detector.planeymax() ),
	detector(detector)
		{};

void TrackHistogrammer::fill(TrackFitResult entry) {
	phi.Fill(entry.phi);
	d0.Fill(entry.d0);
	tanLambda.Fill(entry.tanlambda);
	z0.Fill(entry.z0);
}
