/*
 * CombineTracks.cpp
 *
 *  Created on: Jul 19, 2017
 *      Author: cligtenb
 */

#include <string>

#include "TrackCombiner.cpp"

using namespace std;

//returns correlation factor
void CombineTracks(std::string mimosaInput, std::string timepixInput, int triggerOffset=0,  bool displayEvent=false) {

	TrackCombiner combiner(mimosaInput, timepixInput);
	combiner.openFile("fitResults2.root");
	combiner.processTracks();

}


