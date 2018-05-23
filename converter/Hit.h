#ifndef HIT_H
#define HIT_H
#include <vector>

struct Hit {
	Hit() : row(0), column(0), charge(0) {};
	Hit(short row, short column, short charge) : row(row), column(column), charge(charge) {};
	virtual ~Hit() {};
	unsigned short row, column, charge; //charge=ToT
};

struct TimePixHit : Hit {
	TimePixHit() : Hit(), driftTime(0), ToA(0) {};
	TimePixHit(unsigned short row, unsigned short col, unsigned short ToT, int driftTime, unsigned long long ToA)
		: Hit(row, col, ToT), driftTime(driftTime), ToA(ToA)
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

