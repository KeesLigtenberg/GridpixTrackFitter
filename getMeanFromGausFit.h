/*
 * getMeanFromGausFit.h
 *
 *  Created on: Oct 11, 2017
 *      Author: cligtenb
 */

#ifndef GETMEANFROMGAUSFIT_H_
#define GETMEANFROMGAUSFIT_H_
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

inline double getMeanFromSimpleGausFit( TH1& hist ) {
	TFitResultPtr fitresult=hist.Fit("gaus", "QS");
	if(!fitresult->IsValid()) {std::cerr<< "error: failed to fit histogram "<<hist.GetName()<<std::endl; return 0;}
	return fitresult->Parameter(1);
}
inline double getParameterFromFit( TH1& hist, TF1& fit, int param) {
	fit.SetParameters(hist.GetEntries()/hist.GetNbinsX()*5, hist.GetMean(), hist.GetStdDev() ); //estimates with the correct order of magnitude
	TFitResultPtr fitresult=hist.Fit(&fit, "MSQ"); //More(try to find more than one minimum), Store and Quiet
	if(!fitresult->IsValid()) {
		std::cerr<< "getParameterFromFit error: failed to fit histogram "<<hist.GetName()<<std::endl;
		throw int(1);
	}
	return fitresult->Parameter(param);
}
inline double getMeanFromGausFit( TH1& hist ) {
	TF1 gaus( "myGaus", "[0]*exp(-0.5*((x-[1])/[2])^2)", hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
	int meanParameterNumber=1;
	double mean=0;
	try{
		mean=getParameterFromFit(hist, gaus, meanParameterNumber);
	} catch(const int& error ) {
		if(error==1) {
//			std::cerr<<"retry with gaus default instead"<<std::endl;
//			mean=getMeanFromSimpleGausFit(hist);
			throw error;
		}
	}
	return mean;
}

#endif /* GETMEANFROMGAUSFIT_H_ */
