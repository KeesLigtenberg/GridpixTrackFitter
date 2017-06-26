#include <string>
#include <iostream>

#include "TROOT.h"
#include "DetectorConfiguration.h"
#include "TrackFitter.h"

using namespace std;

const DetectorConfiguration mimosa= {
	6, //planes
	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other!
};


void FitTracks (std::string inputfile) {

	trackFitter telescopeFitter(inputfile, mimosa);

	telescopeFitter.makeMask();

	const int nRepeatFit=5;
	for(int i=0; i<nRepeatFit; i++) {
				cout<<"fitting "<<i<<endl;
				telescopeFitter.fitTracks("residualHistograms.root"); //"+to_string(i)+"
				telescopeFitter.addToShifts( telescopeFitter.getMeans() );
				telescopeFitter.addToAngles( telescopeFitter.getRotations() );
	}

}


//redirect to root main function
int main(int argc, const char* argv[]) {
	gROOT->ProcessLine(".L Hit.h+"); //compile hit
	if(argc>1) FitTracks(argv[1]);
}
