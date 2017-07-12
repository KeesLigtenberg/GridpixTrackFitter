#ifndef HIT_H
#define HIT_H
#include <vector>

struct Hit {
	unsigned short row, column, charge; //charge=ToT
};

struct TimePixHit : Hit {
	TimePixHit() {};
	TimePixHit(unsigned short row, unsigned short col, unsigned short ToT, int driftTime, unsigned long long ToA)
		: Hit{row, col, ToT}, driftTime(driftTime), ToA(ToA)
		  {};

	int driftTime;
	unsigned long long ToA;
};


#pragma link C++ class Hit+;
#pragma link C++ class std::vector<Hit>+;
#pragma link C++ class std::vector<std::vector<Hit>>+;

#pragma link C++ class TimePixHit+;
#pragma link C++ class std::vector<TimePixHit>+;

#endif

