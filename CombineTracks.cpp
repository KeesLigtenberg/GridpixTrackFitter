/*
 * CombineTracks.cpp
 *
 *  Created on: Jul 19, 2017
 *      Author: cligtenb
 */

#include <string>
#include <fstream>

#include "TrackCombiner.cpp"

using namespace std;

//returns correlation factor
void CombineTracks(std::string mimosaInput, std::string timepixInput, bool displayEvent=false) {

	bool doAlignment=false;
	if(!displayEvent and doAlignment) {
		for(int i=0; i<3; ++i) {
			TrackCombiner combiner(mimosaInput, timepixInput);
			cout<<"loading alignment\n";
			combiner.loadAlignment("alignment.dat");
			combiner.openFile("fitResults.root");
			combiner.setDisplayEvent(displayEvent);
			cout<<"processing tracks\n";
			combiner.processTracks();
			cout<<"saving alignment\n";
			combiner.saveAlignment("alignment.dat");
		}
	} else {
		TrackCombiner combiner(mimosaInput, timepixInput);
		combiner.loadAlignment("alignment.dat");
		combiner.openFile(displayEvent ? "tmp.root" : "fitResults.root"); //do not overwrite results if only dipslaying
		combiner.setDisplayEvent(displayEvent);
		combiner.processTracks();
	}
}


