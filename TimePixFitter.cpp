/*
 * TimePixFitter.cpp
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#include "TimePixFitter.h"

#include "TSystem.h"
#include "TView.h"
#include "TRandom.h"

//#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "getObjectFromFile.h"
#include "makeNoisyPixelMask.h"
#include "transformHits.h"
#include "ResidualHistogrammer.h"
#include "TrackFitter.h"

using namespace std;

TimePixFitter::TimePixFitter(std::string inputfile, const DetectorConfiguration& detector) :
	detector(detector),
	houghTransform(detector.xmin(), detector.xmax(), detector.ymin(), detector.ymax(), 12/*12 xbins*/, 12 /*12 ybins*/ ),
	residualHistograms(nullptr),
	hitsCentre(detector.getCentre() ), //initialise to regular centre of sensor
	averageResidualFromSum{0,0}, //initialise to zero
	rotationZFromSum(0),
	slope1FromSum(0), slope2FromSum(0),
	shift(0,0)
{

	//get tree from file
	file=openFile(inputfile);
	hitTable=getObjectFromFile<TTree>("Hits", file);

	//setup tree for reading
	hitTable->SetBranchAddress("timepix", &rawHits);
	hitTable->SetBranchAddress("triggerNumber", &triggerNumber);
	hitTable->SetBranchAddress("trigger", &timestamp);
	nEvents=hitTable->GetEntriesFast();

}

TimePixFitter::~TimePixFitter() {
	file->Close(); //close and write contents of file
	residualHistograms=nullptr; //make sure to delete residualHistograms before trackHistograms, because they use the same file
	trackHistograms=nullptr; //todo: make this independent of order
	//NumberOfHitsInTrack=nullptr;
}

int TimePixFitter::makeMask(double ntimesThreshold) {
	//first mask noisy pixels
	mask= makeNoisyPixelMask(hitTable,ntimesThreshold,{detector.pixelColumns, detector.pixelRows} );
	return 0; //todo:should return number of pixels masked
}

bool TimePixFitter::passEvent(const std::vector<PositionHit>& spaceHit) const {
	const int nMinHits=30;//require at least 4 planes hit
	return spaceHit.size()>nMinHits;
}

int roundToInt(double x) {return x+0.5;};
std::vector<TimePixHit>& addCrossTalk(std::vector<TimePixHit>& hv, double chance, const TimeWalkCorrector& twc) {
	for(auto& h : hv) {
		if(gRandom->Rndm()<chance) {
			//correct timewalk, copy, undo correction
			const double driftScale=25./4096 /*scale*/ * 0.079 /*um/ns*/;
			TimePixHit ct(h);
			auto originalCharge=h.charge;
			double fraction=gRandom->Rndm();
			ct.charge=roundToInt(ct.charge*(1-fraction));
			h.charge=roundToInt(h.charge*fraction);
			double smearing=gRandom->Gaus()*0.2;
			ct.driftTime+=(twc.getCorrection(ct.charge)-twc.getCorrection(originalCharge))/driftScale;
			h.driftTime+=(twc.getCorrection(h.charge)-twc.getCorrection(originalCharge))/driftScale;
//			cout<<ct.driftTime<<" ("<<ct.charge<<") - "<<h.driftTime<<" ("<<h.charge<<")\n";
			int pm= gRandom->Rndm()<0.5 ? -1 : 1;
			if(gRandom->Rndm()<0.5) {
				if(ct.column!=0 && ct.column!=255) ct.column+=pm;
			} else {
				if(ct.row!=0 && ct.row!=255)  ct.row+=pm;
			}
			hv.push_back(ct);
		}
	}
	return hv;
}

std::vector<PositionHit> TimePixFitter::getSpaceHits() {
	std::vector<PositionHit> spaceHit;
	//apply mask
	if (!mask.empty()) {
		auto maskedHit = applyPixelMask(mask, *rawHits);
		//convert rawHits to positions
		const double driftScale=25./4096 /*scale*/ * 0.079 /*mm/ns*/; //Todo: move parameters to detectorConfiguration
		spaceHit=convertHits( maskedHit, detector.pixelsize, detector.pixelsize, driftScale );
	}
	return spaceHit;
}

std::vector<PositionHit> TimePixFitter::getSpaceHitsWithCrossTalk(const TimeWalkCorrector& twc) {
	std::vector<PositionHit> spaceHit;
	//apply mask
	if (!mask.empty()) {
		auto maskedHit = applyPixelMask(mask, *rawHits);
		maskedHit = addCrossTalk(maskedHit, 0.01, twc);
		//convert rawHits to positions
		const double driftScale=25./4096 /*scale*/ * 0.079 /*mm/ns*/;
		spaceHit=convertHits( maskedHit, detector.pixelsize, detector.pixelsize, driftScale );
	}
	return spaceHit;
}

std::vector<PositionHit>&  TimePixFitter::rotateAndShift(
	std::vector<PositionHit>& spaceHit) {
	//apply translation and rotation
	spaceHit = translateHits(spaceHit, shift);
//	spaceHit = rotateHits(spaceHit, angle, hitsCentre);todo: do proper rotation!
	return spaceHit;
}

int TimePixFitter::getEntry(int iEvent) {
	if(iEvent>=nEvents) return false;
	//get entry
	int nb=hitTable->GetEntry(iEvent);
	if (!rawHits) {
		cout << "Did not find hits in tree!" << endl;
		return false;
	}
	return nb;
}

void TimePixFitter::drawEvent(const vector<PositionHit>& spaceHit,
		const std::vector<FitResult3D>& fits) {
	HoughTransformer::drawCluster(spaceHit, detector);
	//			HoughTransformer::drawClusters(houghClusters, detector);
	for (auto& f : fits)
		f.draw(detector.zmin(), detector.zmax());
	gPad->Update();
}

void TimePixFitter::fitTracks(std::string outputfilename) {

	//open histograms
	residualHistograms=unique_ptr<ResidualHistogrammer>(new ResidualHistogrammer(outputfilename, detector));
	if(makeTrackHistograms) trackHistograms=unique_ptr<TrackHistogrammer>(new TrackHistogrammer(detector) );
	numberOfHitsInTrack=new TH1D("numberOfHitsInTrack", "Number of hits in track;Number of hits in track;Tracks", 50,0,400);
	ToTOfHitsInTrack=new TH1D("ToTOfHitsInTrack", "ToT of hits in track;ToT of hits in track [#mus];Hits", 40,0,2);

	//for calculation of com, means and rotation
	double hitsXSum=0, hitsYSum=0;
	double	residualXSum=0, residualYSum=0;
	double rotationZSum=0, rotationZWeightSum=0;
	int nHits=0;
	double slope1Sum=0, slope2Sum=0;

	//loop over all entries
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%100000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		if(!getEntry(iEvent) ) continue;

		//apply mask and convert
		vector<PositionHit> spaceHit = getSpaceHits();

		//optionally do ToT and timewalk corrections here
		//if no ToT correction, it is best to take slightly larger ybins
		houghTransform.xbins=6;

		//optionally set hit errors here, default is all weights equal. See TrackCombiner for example on how to set weights

		//check event
		if( !passEvent(spaceHit) ) continue;
		++nPassed;

		//put sum x and y here to sum all rawHits including noise

		//apply translation
		spaceHit=rotateAndShift(spaceHit); //just a shift in spite of name
		//hough transform
		auto houghClusters = houghTransform(spaceHit);

		//continue if no clusters are found
		if(!houghClusters.size()) {
			continue;
		}

		//first fit on outer planes and translate inner rawHits
		std::vector<FitResult3D> fits;
		for(auto& hitCluster : houghClusters) {

			const unsigned minHitsForFit=houghTransform.minClusterSize; //cut on 10
			if(hitCluster.size()<minHitsForFit) continue; //need at least 2 points for a fit

			//fit track
			auto fit=regressionFit3d(hitCluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}

			//remove outliers
			auto residuals=calculateResiduals(hitCluster, fit);
			hitCluster=cutOnResiduals(hitCluster, residuals, maxResidual);

			//refit on planes with selection
			HoughTransformer::HitCluster selectedHits;
			std::copy_if(hitCluster.begin(), hitCluster.end(), std::back_inserter(selectedHits), selectHitForRefit ); //default: select (copy) all hits
			if(selectedHits.size()<minHitsForFit) continue;
			if(constructLineParallelToZ) fit = makeLinesParallelToZ( selectedHits.front().x, selectedHits.front().y ); //for allignment purposes
			else fit=regressionFit3d(selectedHits); //refit

			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}

			//sum x and y for calculation of centre of mass from track rawHits, also for alignment
			if(hitCluster.size()) for(auto& h : hitCluster) {
				hitsXSum+=h.x;
				hitsYSum+=h.y;
				++nHits;
			}

			residuals=calculateResiduals(hitCluster, fit);

			residualHistograms->fill(residuals, hitsCentre); //fill residuals, hit centre for rotations
			numberOfHitsInTrack->Fill(hitCluster.size());


			//sum residuals
			if(residuals.size()) for(auto& r : residuals) {
//				cout<<r.x<<" "<<r.y<<endl;
				residualXSum+=r.x;
				residualYSum+=r.y;

				auto& rotationPoint=hitsCentre;
				double xc=rotationPoint.first, yc=rotationPoint.second; //x and y center
				double hx=r.h.x-xc, hy=r.h.y-yc;
				double phi=(hy*r.x-hx*r.y)/(hx*hx+hy*hy);
				double weight=hx*hx+hy*hy;
				rotationZSum+=phi*weight; rotationZWeightSum+=weight;
			}

			//sum fit slope
			slope1Sum+=fit.XZ.slope; slope2Sum+=fit.YZ.slope;

			fits.push_back(fit);
			if(makeTrackHistograms) { trackHistograms->fill(fit); }

			++nClusters;
		}

		if(displayEvent) {
			static TCanvas* canv=new TCanvas("eventDisplay", "Display of event", 600,400);
			canv->cd();
			drawEvent(spaceHit, fits);
			if( trackFitter::processDrawSignals()  ) break;
		}

	}

	std::cout<<"passed: "<<nPassed<<"/"<<nEvents<<" with "<<nClusters<<"\n";

	//Recalculate centre of mass of rawHits for each plane
	if(recalculateCOM) {
		cout<<"centre of chip is "<<detector.planexmax()/2<<", "<<detector.planeymax()/2<<endl;
		hitsCentre={ hitsXSum/nHits, hitsYSum/nHits };
		cout<<"Hits centre of mass is: ("<<hitsCentre.first<<", "<<hitsCentre.second<<")"<<endl;
		recalculateCOM=false;
	}

	cout<<"average residual is "<<residualXSum<<"/"<<nHits<<", "<<residualYSum<<"/"<<nHits<<endl;
	averageResidualFromSum= { residualXSum/nHits, residualYSum/nHits };
	rotationZFromSum= rotationZSum/rotationZWeightSum;


	slope1FromSum=slope1Sum/nClusters;
	slope2FromSum=slope2Sum/nClusters;
	cout<<"slopes are "<<slope1FromSum<<" "<<slope2FromSum<<endl;

//	cin.get();

}

std::pair<double, double> TimePixFitter::getMeans() {
//	auto means= residualHistograms->getMeansOfPlanes();
	auto means= averageResidualFromSum;
	cout<<"shift ("<<means.first<<", "<<means.second<<")"<<endl;
	return means;
}

double TimePixFitter::getRotation() {
//	auto rotation= residualHistograms->getRotationOfPlanes();
	auto rotation= rotationZFromSum;
	cout<<"rotation: "<<rotation<<" = "<<rotation/M_PI*180.<<endl;
	return rotation;
}

void TimePixFitter::setShift(
		const std::pair<double, double> & shiftsIn) {
		shift=shiftsIn;
}



void TimePixFitter::addToShift(
		const std::pair<double, double> & shiftsExtra) {
	shift.first+=shiftsExtra.first;
	shift.second+=shiftsExtra.second;
}

const std::pair<double, double>& TimePixFitter::getShift() const {
	return shift;
}


std::pair<double, double> TimePixFitter::getSlopes() const {
	return {slope1FromSum, slope2FromSum};
}

void TimePixFitter::setSlopes(std::pair<double, double> slopes) {
	houghTransform.angleOfTracksX=slopes.first;
	houghTransform.angleOfTracksY=slopes.second;
}

