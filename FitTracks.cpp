#include <string>
#include <iostream>

#include "TROOT.h"
#include "TGraph.h"
#include "TH1.h"

#include "AllCombiner.h"
#include "DetectorConfiguration.h"
#include "TrackFitter.h"

//#include "/user/cligtenb/rootmacros/AllCombiner.h"

#if 1 //root?
#include "TrackFitter.cpp"
#include "linearRegressionFit.cpp"
#include "makeNoisyPixelMask.cpp"
#include "ResidualHistogrammer.cpp"
#endif

using namespace std;

const DetectorConfiguration mimosa= {
	6, //planes
//	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
	{0, 15.8, 31.8, 143.1, 161.55, 179.91 }, //plane positions as measured
	0.0184, //pixelsize in mm
	1153, 577 //row, column, ://one extra because exampleData starts at 1 and our data starts at 0 //TODO: change this to one or the other!
};

//std::vector<std::pair<double, double>>& fixFirstPosition(std::vector<std::pair<double, double>>& means) {
//	for(int i=1; i<means.size(); ++i) {
//		means[i].first+=means[0].first;
//		means[i].second+=means[0].second;
//	}
//	means[0]=std::make_pair(0.,0.);
//	return means;
//}

void FitTracks (std::string inputfile, int nRepeatFit=5) {

	trackFitter telescopeFitter(inputfile, mimosa);

	telescopeFitter.makeMask(5e3);
//	telescopeFitter.setShifts( { {0,0}, {0,0}, {0,0}, {0,0}, {0,1}, {0,0} } );

	//initialise alignment parameters
	double recursion[nRepeatFit],
		shiftx[mimosa.nPlanes][nRepeatFit],
		shifty[mimosa.nPlanes][nRepeatFit],
		rotation[mimosa.nPlanes][nRepeatFit];

	for(int i=0; i<nRepeatFit; i++) {
				cout<<"fitting "<<i<<endl;
				if(i==4) telescopeFitter.displayEvent=true;

				if(i>=3) telescopeFitter.makeTrackHistograms=true;

				if(i<=2) telescopeFitter.selectHitForRefit=[](const PositionHit& h) {return h.plane==1 || h.plane==4;};
				else telescopeFitter.selectHitForRefit=[](const PositionHit& h) {return true;};

				if(i==1) telescopeFitter.constructLineParallelToZ=true;
				else telescopeFitter.constructLineParallelToZ=false;

				//fit tracks!
				telescopeFitter.fitTracks("residualHistograms"+to_string(i)+".root");

				auto means=telescopeFitter.getMeans();
				auto rotations=telescopeFitter.getRotations();

				means[1]=means[4]={0,0};
				rotations[1]=0;

				switch(i) {
				case 0:
					telescopeFitter.addToShifts( means );
					telescopeFitter.setAngles({0,0,0,0,0,0});
					break;
				case 1:
					telescopeFitter.addToAngles( {0,0,0,0,rotations[4],0} );
					break;
				case 2:
					telescopeFitter.addToAngles( rotations );
					break;
				case 3:
					telescopeFitter.setSlopes( telescopeFitter.getSlopes() );
					break;
				default:
					telescopeFitter.addToShifts( means );
					telescopeFitter.addToAngles( rotations );
					break;
				}

				//store alignment parameters
				recursion[i]=i;
				for(int plane=0; plane<mimosa.nPlanes; ++plane) {
					shiftx[plane][i]=telescopeFitter.getShifts().at(plane).first;
					shifty[plane][i]=telescopeFitter.getShifts().at(plane).second;
					rotation[plane][i]=telescopeFitter.getAngles().at(plane);
				}
	}


	//create and combine graphs
	std::vector<TGraph*> shiftxGraph,shiftyGraph, rotationGraph;
	for(int plane=0;plane<mimosa.nPlanes; ++plane) {
		shiftxGraph.push_back( new TGraph(nRepeatFit, recursion, shiftx[plane]) );
		shiftyGraph.push_back( new TGraph(nRepeatFit, recursion, shifty[plane]) );
		rotationGraph.push_back( new TGraph(nRepeatFit, recursion, rotation[plane]) );

		shiftxGraph.back()->SetTitle( ("plane "+to_string(plane+1)+";Iteration;Correction [mm]").c_str() );
		shiftyGraph.back()->SetTitle( ("plane "+to_string(plane+1)+";Iteration;Correction [mm]").c_str() );
		rotationGraph.back()->SetTitle( ("plane "+to_string(plane+1)+";Iteration;Correction [rad.]").c_str() );
	}

	AllCombiner<TGraph>  Xcombination("shiftsxCombined", shiftxGraph);
	AllCombiner<TGraph>  Ycombination("shiftsyCombined", shiftyGraph);
	AllCombiner<TGraph>  RotCombination("rotationCombined", rotationGraph);

	std::vector<TCanvas*> canv;
	for(AllCombiner<TGraph>* comb : {&Xcombination, &Ycombination, &RotCombination}) {
		comb->setStyle(7);
		canv.push_back( comb->createCombined() );
	}

}


//redirect to root main function
int main(int argc, const char* argv[]) {
	gROOT->ProcessLine(".L Hit.h+"); //compile hit
	if(argc>1) FitTracks(argv[1]);
}
