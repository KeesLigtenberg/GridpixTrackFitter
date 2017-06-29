#include <string>
#include <iostream>

#include "TROOT.h"
#include "TGraph.h"

#include "DetectorConfiguration.h"
#include "TrackFitter.h"

#include "/user/cligtenb/rootmacros/AllCombiner.h"

#if 1 //root?
#include "TrackFitter.cpp"
#include "linearRegressionFit.cpp"
#include "makeNoisyPixelMask.cpp"
#include "ResidualHistogrammer.cpp"
#endif

using namespace std;

const DetectorConfiguration mimosa= {
	6, //planes
	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other!
};


void FitTracks (std::string inputfile) {

	trackFitter telescopeFitter(inputfile, mimosa);

	telescopeFitter.makeMask(5e3);

	const int nRepeatFit=5;

	//initialise alignment parameters
	double recursion[nRepeatFit],
		shiftx[mimosa.nPlanes][nRepeatFit],
		shifty[mimosa.nPlanes][nRepeatFit],
		rotation[mimosa.nPlanes][nRepeatFit];

	for(int i=0; i<nRepeatFit; i++) {
				cout<<"fitting "<<i<<endl;
//				if(i==4) telescopeFitter.displayEvent=true;
				if(i==4) telescopeFitter.makeTrackHistograms=true;

				telescopeFitter.fitTracks("residualHistograms"+to_string(i)+".root"); //
				telescopeFitter.addToShifts( telescopeFitter.getMeans() );
				telescopeFitter.addToAngles( telescopeFitter.getRotations() );

				//store alignment parameters
				recursion[i]=i;
				for(int plane=0; plane<mimosa.nPlanes; ++plane) {
					shiftx[plane][i]=telescopeFitter.getShifts().at(plane).first;
					shifty[plane][i]=telescopeFitter.getShifts().at(plane).second;
					rotation[plane][i]=telescopeFitter.getAngles().at(plane);
				}
	}


	//create and combine graphs
//	std::vector<TGraph*> shiftxGraph,shiftyGraph, rotationGraph;
//	for(int plane=0;plane<mimosa.nPlanes; ++plane) {
//		shiftxGraph.push_back( new TGraph(nRepeatFit, recursion, shiftx[plane]) );
//		shiftyGraph.push_back( new TGraph(nRepeatFit, recursion, shifty[plane]) );
//		rotationGraph.push_back( new TGraph(nRepeatFit, recursion, rotation[plane]) );
//		for(const std::vector<TGraph*>& v: {shiftxGraph, shiftyGraph, rotationGraph}) {
//			v.back()->SetTitle( ("plane "+to_string(plane+1)).c_str() );
//		}
//	}
//
//	AllCombiner<TGraph>  Xcombination("shiftsxCombined", shiftxGraph);
//	AllCombiner<TGraph>  Ycombination("shiftsyCombined", shiftyGraph);
//	AllCombiner<TGraph>  RotCombination("rotationCombined", rotationGraph);
//
//	std::vector<TCanvas*> canv;
//	for(AllCombiner<TGraph>* comb : {&Xcombination, &Ycombination, &RotCombination}) {
//		comb->setStyle(7);
//		canv.push_back( comb->createCombined() );
//	}

}


//redirect to root main function
int main(int argc, const char* argv[]) {
	gROOT->ProcessLine(".L Hit.h+"); //compile hit
	if(argc>1) FitTracks(argv[1]);
}
