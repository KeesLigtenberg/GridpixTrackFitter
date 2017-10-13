/*
 * Alignment.h
 *
 *  Created on: Oct 10, 2017
 *      Author: cligtenb
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <vector>
#include <iostream>
#include <array>

#include "TTree.h"
#include "TProfile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "getHistFromTree.h"
#include "PositionHit.h"
#include "getMeanFromGausFit.h"

class TimeWalkCorrector {
public:
	TimeWalkCorrector() : minToT(0.), coefficient(0.), shiftToT(0.), shiftx(0.) {}
	TimeWalkCorrector(double minToT, double coefficient, double shiftt, double shiftx) :
		minToT(minToT), shiftToT(shiftt), shiftx(shiftx), coefficient(coefficient)
	{}

	std::vector<PositionHit>&  correct(std::vector<PositionHit>& spaceHit) const;

	void calculate(TTree* tree);

	void save(std::ostream&) const;
	void load(std::istream&);
private:
	double minToT, shiftToT, shiftx, coefficient;
};


struct RelativeAligner {
	RelativeAligner(const TVector3& shift,const TVector3& timepixCOM, double xangle, double yangle, double zangle ) :
		shift{shift}, timepixCOM{timepixCOM}, angle{xangle,yangle, zangle} {}
	RelativeAligner() : shift(0,0,0), timepixCOM{11,-7,6}, angle{0,0,0} {};

	void calculate(TTree* tree);

	void save(std::ostream&) const;
	void load(std::istream&);

	TVector3 shift, timepixCOM; //com in telescope frame!
	std::array<double, 3> angle; //X, Y, Z
};

struct Alignment {
	Alignment() {};
	Alignment(std::string filename) {
		std::ifstream fin(filename);
		load(fin);
	}

	void load(std::istream& fin) {
		timeWalkCorrection.load(fin);
		relativeAlignment.load(fin);
	}

	void save(std::ostream& fout) const {
		timeWalkCorrection.save(fout);
		relativeAlignment.save(fout);
	}
	void saveToFile(std::string filename) const {
		std::ofstream fout(filename);
		timeWalkCorrection.save(fout);
		relativeAlignment.save(fout);
	}

	//class for keeping all alignment
	TimeWalkCorrector timeWalkCorrection;
	RelativeAligner relativeAlignment;
};


std::vector<PositionHit>& TimeWalkCorrector::correct(std::vector<PositionHit>& spaceHit) const {
	spaceHit.erase(
			std::remove_if(spaceHit.begin(), spaceHit.end(), [this](const PositionHit&h) {return (h.ToT*0.025) <minToT;} ),
			spaceHit.end()
	);
	for(auto& h : spaceHit) {
		h.x=h.x-coefficient/(h.ToT*0.025+shiftToT)+shiftx;
	}
	return spaceHit;
}

void TimeWalkCorrector::save(std::ostream& output) const {
	output<<"TIMEWALKPARAMETERS\n"
		  <<minToT<<" "<<shiftToT<<" "<<shiftx<<" "<<coefficient<<"\n";
}

void TimeWalkCorrector::load(std::istream& input) {
	std::string header;
	input>>header;
	if(header!="TIMEWALKPARAMETERS" or not input.good()) {
		std::cerr<<"failed to read header TIMEWALKPARAMETERS\n"; throw 1;
	}
	input>>minToT>>shiftToT>>shiftx>>coefficient;
	if(not input.good() ) {
		std::cerr<<"failed to read parameters\n";
		throw 1;
	}
}

void TimeWalkCorrector::calculate(TTree* tree) {
		TH1* hist=getHistFromTree(tree, "timepixHits.rx:timepixHits.ToT*0.025", "", "profile", "profgoff");
//		gPad->Update();
//		std::cin.get();

		TProfile* prof = dynamic_cast<TProfile*>(hist);

		if(!prof) {
			std::cerr<<"Could not get profile from tree!\n";
			throw 1;
		}

		auto scoef=std::to_string(coefficient);
		auto sToT=std::to_string(shiftToT);
		//subtract old correction in fit of new one:
		TF1* fun=new TF1("fun", ("[0]/(x+[1])+[2]-"+scoef+"/(x+"+sToT+")").c_str() );

		TFitResultPtr fit=prof->Fit(fun, "SQ", "", minToT, 100);
		if(!fit->IsValid()) { std::cerr<<"Time walk fit failed!\n"; throw 1;}

		std::cout<<"Updated time walk parameters\n"
				 <<"coefficient "<<coefficient<<" -> "<<fit->Parameter(0)<<"\n"
				 <<"shiftToT "<<shiftToT<<" -> "<<fit->Parameter(1)<<"\n";
		coefficient=fit->Parameter(0);
		shiftToT=fit->Parameter(1);
	}

void RelativeAligner::calculate(TTree* tree) {
	//calculate COM
	//fill histogram with average from each track and get average from that, i.e. weight each track, not each hit
	for(int i=0; i<3; i++) {
		const auto x = std::array<std::string,3>{"x", "y", "z"}[i];
		auto xhist=getHistFromTree(tree, "Sum$(timepixHits."+x+")/Length$(timepixHits)", "", x+"hist", "goff" );
		double xmean=getMeanFromGausFit(*xhist);
		timepixCOM[i]=xmean;
	}
	std::cout<<"com="<<timepixCOM[0]<<" "<<timepixCOM[1]<<" "<<timepixCOM[2]<<"\n";
	for(int i=0; i<2;i++) {
		const std::string axis= i ? "Y" : "X";
		//intercepts
		auto zstring=std::to_string( timepixCOM.z() );//mean in Telescope frame -380+~6 =~-374
		auto shiftHist=getHistFromTree(tree, "(telescopeFits[]."+axis+"Z.intercept+telescopeFits[]."+axis+"Z.slope*"+zstring+")-"
			 "(timepixFits[0]."+axis+"Z.intercept+timepixFits[0]."+axis+"Z.slope*"+zstring+")", "", "shiftHist"+axis, "goff" );

//		gPad->Update();
//		std::cin.get();

		double shiftMean=getMeanFromGausFit(*shiftHist);
		shift[i]+=shiftMean;

		std::cout<<"additional shift is "<<shiftMean<<"\n";

		//angles
		auto hist=getHistFromTree(tree, "timepixFits."+axis+"Z.slope", "fabs(timepixFits."+axis+"Z.slope)<0.2", "slopeHist"+axis, "goff" );
		double mean=getMeanFromGausFit(*hist);
		angle[i]+= std::atan( i==0 ? mean : -mean );

		std::cout<<"added rotation is "<<std::atan(mean)<<"\n";
	}
	//Z axis rotation
	auto comx=std::to_string(timepixCOM.x()), comy=std::to_string(timepixCOM.y());
	auto zRotation=getHistFromTree(tree,
			"-dyz*(Sum$(timepixHits.x)/Length$(timepixHits)-"+comx+")"
			"+dxz*(Sum$(timepixHits.y)/Length$(timepixHits)-"+comy+")",
			"", "zRotHist", "goff");
	double zAngle=getMeanFromGausFit(*zRotation);
	angle[2]= zAngle;
	std::cout<<"zrotation is "<<zAngle<<"\n";
}

void RelativeAligner::load(std::istream& input) {
	std::string header;
	input>>header;
	if(header!="RELATIVEALIGNMENT" or not input.good()) {
		std::cerr<<"failed to read header RELATIVEALIGNMENT\n"; throw 1;
	}
	double x,y,z;
	input>>x>>y>>z>>angle[0]>>angle[1]>>angle[2];
	shift.SetXYZ(x,y,z);
	if(not input.good() ) {
		std::cerr<<"failed to read parameters\n";
		throw 1;
	}
}

void RelativeAligner::save(std::ostream& output) const {
	output<<"RELATIVEALIGNMENT\n";
	for(int i=0; i<3; i++) output<<shift[i]<<" ";
	output<<"\n";
	for(int i=0; i<3; i++) output<<angle[i]<<" ";
	output<<"\n";
}



#endif /* ALIGNMENT_H_ */
