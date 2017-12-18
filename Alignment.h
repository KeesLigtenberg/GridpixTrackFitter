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
#include <fstream>
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
	TimeWalkCorrector() : minToT(0.1), coefficient(0.06), shiftToT(0.4), shiftx(0.6) {}
	TimeWalkCorrector(double minToT, double coefficient, double shiftt, double shiftx) :
		minToT(minToT), shiftToT(shiftt), shiftx(shiftx), coefficient(coefficient)
	{}

	std::vector<PositionHit>&  correct(std::vector<PositionHit>& spaceHit) const;

	void calculate(TTree* tree);

	void save(std::ostream&) const;
	void load(std::istream&);
private:
	double minToT, coefficient, shiftToT, shiftx;
};


struct RelativeAligner {
	RelativeAligner(const TVector3& shift,const TVector3& timepixCOM, double xangle, double yangle, double zangle ) :
		shift{shift}, timepixCOM{timepixCOM}, angle{xangle,yangle, zangle} {}
	RelativeAligner() : shift(0,0,0), timepixCOM{11,-7,6}, angle{0,0,0} {};

	void calculate(TTree* tree);

	void save(std::ostream&) const;
	void load(std::istream&);

	TVector3 getCOM() const { return timepixCOM-shift; };//COM in mirrored timepix frame!

	TVector3 shift, timepixCOM; //com in telescope frame!
	std::array<double, 3> angle; //X, Y, Z
};

struct MimosaAligner {
	void save(std::ostream&) const;
	void load(std::istream&);

	std::vector<std::pair<double,double>> shifts;
	std::vector<std::pair<double,double>> COMs;
	std::vector<double> angles;
	std::pair<double,double> slopes;

	int nplanes=6;
};

struct ToTCorrector {
	void load(std::istream&);
	void load(std::string filename);

	std::vector<PositionHit>&  correct(std::vector<PositionHit>& spaceHit) const;

	static const int ncols=256;
	int firstCol, lastCol;
	std::array<double, ncols> correction;
};

struct Alignment {
	//class for keeping all alignment
	Alignment() {};
	Alignment(std::string filename) {
		std::ifstream fin(filename);
		load(fin);
	}

	void load(std::istream& fin) {
		timeWalkCorrection.load(fin);
		relativeAlignment.load(fin);
		mimosa.load(fin);
	}

	void save(std::ostream& fout) const {
		timeWalkCorrection.save(fout);
		relativeAlignment.save(fout);
	}
	void saveToFile(std::string filename) const {
		std::ofstream fout(filename);
		timeWalkCorrection.save(fout);
		relativeAlignment.save(fout);
		mimosa.save(fout);
	}

	TimeWalkCorrector timeWalkCorrection;
	RelativeAligner relativeAlignment;
	MimosaAligner mimosa;
};

//template<class T>
//std::istream& operator>>(std::istream& is, std::pair<T,T> p) {
//	return is>>p.first>>p.second;
//}
void checkHeader( std::istream& input, std::string name) {
	std::string header;
	input>>header;
	if(header!=name or not input.good()) {
		std::cerr<<"failed to read header "<<name<<"\n"; throw 1;
	}
}
std::vector<PositionHit>& TimeWalkCorrector::correct(std::vector<PositionHit>& spaceHit) const {
//	spaceHit.erase(
//			std::remove_if(spaceHit.begin(), spaceHit.end(), [this](const PositionHit&h) {return (h.ToT*0.025) <minToT;} ),
//			spaceHit.end()
//	);
	for(auto& h : spaceHit) {
		if(h.ToT/40.<minToT) {h.flag=-3;};
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
		//require to have passed cuts, except rx cut
		TH1* hist=getHistFromTree(tree, "timepixHits.rx:timepixHits.ToT*0.025", "timepixHits.flag>0 || timepixHits.flag==-1", "profile", "profgoff");
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
	std::string eventcuts="timepixClusterSize/ntimepixHits>0.75 && timepixClusterSize>=30";
	//calculate COM
	//fill histogram with average from each track and get average from that, i.e. weight each track, not each hit
	for(int i=0; i<3; i++) {
		const auto x = std::array<std::string,3>{"x", "y", "z"}[i];
		auto xhist=getHistFromTree(tree, "Sum$(timepixHits."+x+")/Length$(timepixHits)", eventcuts, x+"hist", "goff" );
//		gPad->Update();
//		std::cin.get();
		//todo: correct COM with actual center of mass from only passing hits
		double xmean=getMeanFromGausFit(*xhist);
		timepixCOM[i]=xmean;
	}
	std::cout<<"com="<<timepixCOM[0]<<" "<<timepixCOM[1]<<" "<<timepixCOM[2]<<"\n";
	for(int i=0; i<2;i++) {
		const std::string axis= i ? "Y" : "X";
		//intercepts
		auto zstring=std::to_string( timepixCOM.z() );//mean in Telescope frame -380+~6 =~-374
		auto shiftHist=getHistFromTree(tree, "(telescopeFits[]."+axis+"Z.intercept+telescopeFits[]."+axis+"Z.slope*"+zstring+")-"
			 "(timepixFits[0]."+axis+"Z.intercept+timepixFits[0]."+axis+"Z.slope*"+zstring+")", eventcuts, "shiftHist"+axis, "goff" );
//		gPad->Update();
//		std::cin.get();
		double shiftMean=getMeanFromGausFit(*shiftHist);
		shift[i]+=shiftMean;

		std::cout<<"additional shift is "<<shiftMean<<"\n";

		//angles
		const std::string rotaxis= i ? "X" : "Y";
		auto hist=getHistFromTree(tree, "timepixFits."+rotaxis+"Z.slope", "fabs(timepixFits."+rotaxis+"Z.slope)<0.2 && "+eventcuts, "slopeHist"+rotaxis, "goff" );
//		gPad->Update();
//		std::cin.get();
		double mean=getMeanFromGausFit(*hist);
		angle[i]+= std::atan( i ? -mean : mean );

		std::cout<<"added "<<axis<<" rotation is "<<std::atan(mean)<<"\n";
	}
	//Z axis rotation
	auto comx=std::to_string(timepixCOM.x()), comy=std::to_string(timepixCOM.y());
	std::string avgx="Sum$(timepixHits.x)/Length$(timepixHits)", avgy="Sum$(timepixHits.y)/Length$(timepixHits)";
	auto zRotation=getHistFromTree(tree,
			"(-dyz*("+avgx+"-"+comx+")+dxz*("+avgy+"-"+comy+"))"
			"/(pow("+avgx+"-"+comx+", 2)+pow("+avgy+"-"+comy+",2))",
			"(pow("+avgx+"-"+comx+", 2)+pow("+avgy+"-"+comy+",2))"
			"*(fabs("
			"(-dyz*("+avgx+"-"+comx+")"
			"+dxz*("+avgy+"-"+comy+"))"
			"/(pow("+avgx+"-"+comx+", 2)+pow("+avgy+"-"+comy+",2))"
			"))<1 && "+eventcuts
			, "zRotHist", "goff");
//	std::cout<<"got histogram zrotation\n";
	double zAngle=getMeanFromGausFit(*zRotation);
	angle[2]-=zAngle;
	std::cout<<"added zrotation is "<<zAngle<<"\n";
}

void RelativeAligner::load(std::istream& input) {
	std::string header;
	input>>header;
	if(header!="RELATIVEALIGNMENT" or not input.good()) {
		std::cerr<<"failed to read header RELATIVEALIGNMENT\n"; throw 1;
	}
	double x,y,z;
	double comx, comy, comz;
	input>>x>>y>>z
		 >>comx>>comy>>comz
		 >>angle[0]>>angle[1]>>angle[2];
	shift.SetXYZ(x,y,z);
	timepixCOM.SetXYZ(comx,comy,comz);
	if(not input.good() ) {
		std::cerr<<"failed to read parameters\n";
		throw 1;
	}
}

void RelativeAligner::save(std::ostream& output) const {
	output<<"RELATIVEALIGNMENT\n";
	for(const auto v : {&shift, &timepixCOM} ) {
		for(int i=0; i<3; i++) output<<(*v)[i]<<" ";
		output<<"\n";
	}
	for(int i=0; i<3; i++) output<<angle[i]<<" ";
	output<<"\n";
}

void MimosaAligner::save(std::ostream& os) const {
	os<<"MIMOSAALIGNMENT\n";
	for(const auto& v : {shifts, COMs} ) {
		for(const auto p : v ) {
			os<<p.first<<" "<<p.second<<"\n";
		}
	}
	for(const auto& i : angles) {
		os<<i<<"\n";
	}
	os<<slopes.first<<" "<<slopes.second<<"\n";
}

void MimosaAligner::load(std::istream& input) {
	std::string header;
	input>>header;
	if(header!="MIMOSAALIGNMENT" or not input.good()) {
		std::cerr<<"failed to read header MIMOSAALIGNMENT\n"; throw 1;
	}
	for(auto v : {&shifts, &COMs} ) {
		v->clear();
		for(int i=0; i<nplanes; ++i ) {
			double f, s;
			input>>f>>s;
			v->push_back({f,s});
		}
	}
	angles.clear();
	for(int i=0; i<nplanes; ++i) {
		double a;
		input>>a;
		angles.push_back(a);
	}
	input>>slopes.first>>slopes.second;
}

void ToTCorrector::load(std::istream& in) {
	checkHeader(in, "TOTCORRECTION");
	int nColsInFile=0;
	in>>nColsInFile;
	if(nColsInFile != ncols) {
		std::cerr<<"Found "<<nColsInFile<<" in file, but expected "<<ncols<<"\n";
	}
	in>>firstCol>>lastCol;
	for(int i=0; i<ncols; i++) {
		in>>correction[i];
	}
}

void ToTCorrector::load(std::string filename) {
	std::ifstream in(filename.c_str());
	if(in.good()){
		load(in);
	} else {
		std::cout<<"Did not find ToTCorrection: not using any factor!\n";
		firstCol=0; lastCol=ncols;
		std::fill(correction.begin(), correction.end(), 1.);
	}
}

std::vector<PositionHit>& ToTCorrector::correct(
		std::vector<PositionHit>& spaceHit) const {
	for(auto& h : spaceHit) {
		double factor=correction[ std::min(std::max(h.column,firstCol), lastCol) ];
		h.ToT/=factor;
	}
	return spaceHit;
}

std::vector<FitResult3D> transformFitsToTimepixFrame( const std::vector<FitResult3D>& fits, const RelativeAligner& ra) {
	std::vector<FitResult3D> fitsInTimePixFrame;
	for(auto& f :fits ) //telescopeTPCLines)
		fitsInTimePixFrame.push_back(
				f.makeShifted(-ra.shift)
				 .makeRotated(-ra.angle[2], ra.getCOM(), {0,0,1})
				 .makeRotated(-ra.angle[0], ra.getCOM(), {1,0,0})
				 .makeRotated(-ra.angle[1], ra.getCOM(), {0,1,0})
				 .makeMirrorY() ); //
	return fitsInTimePixFrame;
}

PositionHit& transformHitToTimepixFrame( PositionHit& h, const RelativeAligner& ra) {
	h.SetPosition(h.getPosition() - ra.shift);
	h.RotatePosition(-ra.angle[2], ra.getCOM(), {0,0,1});
	h.RotatePosition(-ra.angle[0], ra.getCOM(), {1,0,0});
	h.RotatePosition(-ra.angle[1], ra.getCOM(), {0,1,0});
	h.y=-h.y;
	return h;
}

std::vector<PositionHit> transformHitsToTimepixFrame( const std::vector<PositionHit>& hits, const RelativeAligner& ra) {
	auto hitsInTimepixFrame=hits;
	for(auto& h: hitsInTimepixFrame) {
		h=transformHitToTimepixFrame(h, ra);
	}
	return hitsInTimepixFrame;
}

#endif /* ALIGNMENT_H_ */
