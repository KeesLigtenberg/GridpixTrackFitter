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
		detector(detector)
{
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

//get x,y mean
std::vector<std::pair<double, double> > ResidualHistogrammer::getMeansOfPlanes() {
	std::vector<std::pair<double, double>> vec;
	for(auto& p : planeHist) {
		vec.push_back(p.getMeansFromFit());
//		vec.push_back( {p.xResidual.GetMean(), p.yResidual.GetMean()} );
	}
	return vec;
}

double getMeanFromSimpleGausFit( TH1& hist ) {
	TFitResultPtr fitresult=hist.Fit("gaus", "QS");
	if(!fitresult->IsValid()) {std::cerr<< "error: failed to fit histogram "<<hist.GetName()<<std::endl; return 0;}
	return fitresult->Parameter(1);
}
double getParameterFromFit( TH1& hist, TF1& fit, int param) {
	fit.SetParameters(hist.GetEntries()/hist.GetNbinsX()*2, hist.GetMean(),2*hist.GetBinWidth(1) ); //estimates with the correct order of magnitude
	TFitResultPtr fitresult=hist.Fit(&fit, "MSQ"); //More(try to find more than one minimum), Store and Quiet
	if(!fitresult->IsValid()) {
		std::cerr<< "getParameterFromFit error: failed to fit histogram "<<hist.GetName()<<std::endl;
		throw int(1);
	}
	return fitresult->Parameter(param);
}
double getMeanFromGausFit( TH1& hist ) {
	TF1 gaus( "myGaus", "[0]*exp(-0.5*((x-[1])/[2])^2)", hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
	int meanParameterNumber=1;
	double mean=0;
	try{
		mean=getParameterFromFit(hist, gaus, meanParameterNumber);
	} catch(const int& error ) {
		if(error==1) {
			std::cerr<<"retry with gaus default instead"<<std::endl;
			mean=getMeanFromSimpleGausFit(hist);
		}
	}
	return mean;
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
