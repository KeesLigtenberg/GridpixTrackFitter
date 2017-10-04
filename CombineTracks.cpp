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
void CombineTracks(std::string mimosaInput, std::string timepixInput, bool displayEvent=false) {

	TrackCombiner combiner(mimosaInput, timepixInput);
	if(not displayEvent) combiner.openFile("fitResults.root");
	combiner.setDisplayEvent(displayEvent);
	combiner.processTracks();

}


