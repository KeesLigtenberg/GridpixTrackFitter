//root macro

#include "TimePixEventTreeReader.h"
#include "../trackFitter/getObjectFromFile.h"

void TimePixConvert( std::string inputfile ) {
	TTree* events=getObjectFromFile<TTree>("events", inputfile);
	TimePixTreeReader tr(events);
	inputfile.insert(inputfile.length()-5, "_converted");
	tr.Loop( inputfile );

}
