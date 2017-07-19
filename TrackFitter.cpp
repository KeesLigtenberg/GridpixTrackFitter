/*
 * trackFitter.cpp
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#include "TrackFitter.h"

#include "TSystem.h"
#include "TView.h"

//#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "getObjectFromFile.h"
#include "linearRegressionFit.h"
#include "makeNoisyPixelMask.h"
#include "transformHits.h"
#include "ResidualHistogrammer.h"

using namespace std;

trackFitter::trackFitter(std::string inputfile, const DetectorConfiguration& detector) :
	detector(detector),
	houghTransform(detector.xmin(), detector.xmax(), detector.ymin(), detector.ymax(), 30/*xbins*/, 15 /*ybins*/ ),
	residualHistograms(nullptr),
	hitsCentre(detector.nPlanes, detector.getCentre() ), //initialise to regular centre of sensor
	averageResidualFromSum(detector.nPlanes, {0,0} ), //initialise to zero
	rotationZFromSum(detector.nPlanes),
	slope1FromSum(0), slope2FromSum(0)
{

	//get tree from file
	file=openFile(inputfile);
	hitTable=getObjectFromFile<TTree>("Hits", file);

	//setup tree for reading
	hitTable->SetBranchAddress("mimosa", &mimosaHit);
	//    unsigned short triggerNumberBegin, triggerNumberEnd;
//	hitTable->SetBranchAddress("triggerNumberBegin", &triggerNumberBegin);
//	hitTable->SetBranchAddress("triggerNumberEnd", &triggerNumberEnd);

}

trackFitter::~trackFitter() {
	file->Close();
	residualHistograms=nullptr; //make sure to close residualHistograms before track, because they use the same file
	trackHistograms=nullptr; //todo: combine both in a proper way;
}

int trackFitter::makeMask(double ntimesThreshold) {
	//first mask noisy pixels
	for(int iplane=0; iplane<detector.nPlanes; iplane++) {
		mask.emplace_back( makeNoisyPixelMask(hitTable,iplane,ntimesThreshold,{detector.pixelColumns, detector.pixelRows} ) );
	}
	return 0; //should return number of pixels masked
}

bool trackFitter::passEvent(std::vector<std::vector<PositionHit> > spaceHit) {
	int nPlanesHit=0;
	int nTotalHits=0;
	for(auto& v:spaceHit) {
		int size=v.size();
		if(size) {
			nPlanesHit++;
			nTotalHits+=size;
		}
	}
	const int nMinPlanesHit=3;//require at least 4 planes hit
	return nPlanesHit>nMinPlanesHit;
}

std::vector<std::vector<PositionHit> > trackFitter::getSpaceHits() {
	std::vector<std::vector<PositionHit> > spaceHit;
	//apply mask
	if (!mask.empty()) {
		auto itmask = mask.begin();
		for (unsigned plane = 0; plane < mimosaHit->size(); plane++) {
			auto maskedHit = applyPixelMask(*itmask++, mimosaHit->at(plane));
			//convert hits to positions
			spaceHit.push_back(
					convertHits(maskedHit, detector.planePosition[plane],
							detector.pixelsize, detector.pixelsize, plane));
		}
	}
	return spaceHit;
}

std::vector<std::vector<PositionHit> >&  trackFitter::rotateAndShift(
		std::vector<std::vector<PositionHit> >& spaceHit) {
	//apply translation and rotation
	if (!shifts.empty())
		spaceHit = translateHits(spaceHit, shifts);
	if (!angles.empty())
		spaceHit = rotateHits(spaceHit, angles, hitsCentre);

}

int trackFitter::getEntry(int iEvent) {
	//get entry
	int nb=hitTable->GetEntry(iEvent);
	if (!mimosaHit) {
		cout << "Did not find mimosaHit!" << endl;
		return false;
	}
	return nb;
}

void trackFitter::fitTracks(std::string outputfilename) {
	residualHistograms=unique_ptr<ResidualHistogrammer>(new ResidualHistogrammer(outputfilename, detector));
	if(makeTrackHistograms) trackHistograms=unique_ptr<TrackHistogrammer>(new TrackHistogrammer(detector) );

	//for calculation of com, means and rotation
	std::vector<double> hitsXSum(detector.nPlanes), hitsYSum(detector.nPlanes);
	std::vector<double>	residualXSum(detector.nPlanes), residualYSum(detector.nPlanes);
	std::vector<double> rotationZSum(detector.nPlanes), rotationZWeightSum(detector.nPlanes);
	std::vector<int> nHits(detector.nPlanes);
	double slope1Sum=0, slope2Sum=0;

	//loop over all entries
	const long long nEvents=hitTable->GetEntriesFast(); //std::min( (long long) 2,);
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%100000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		if(!getEntry(iEvent) ) continue;

		//apply mask and convert
		std::vector<std::vector<PositionHit> > spaceHit = getSpaceHits();

		//check event
		if( !passEvent(spaceHit) ) continue;
		++nPassed;

		//sum x and y here to sum all hits including noise

		//apply translation and rotation
		spaceHit=rotateAndShift(spaceHit);
		//hough transform
		auto houghClusters = houghTransform(spaceHit);

		//require at least nmin planes to be hit
		const int nMinPlanesHit=4;
		houghClusters.remove_if([](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
		if(!houghClusters.size()) {
			continue;
		}

		//first fit on outer planes and translate inner hits
		std::vector<SimpleFitResult> fits;
		for(auto& hitCluster : houghClusters) {

			if(hitCluster.size()<2 || hitCluster.getNPlanesHit()<=1) continue;

			//fit track
			auto fit=linearRegressionFit(hitCluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}

			//remove outliers
			auto residuals=calculateResiduals(hitCluster, fit);
			hitCluster=cutOnResiduals(hitCluster, residuals, maxResidual);

			//refit on planes with selection
			HoughTransformer::HitCluster selectedHits;
			std::copy_if(hitCluster.begin(), hitCluster.end(), std::back_inserter(selectedHits), selectHitForRefit );
			if(selectedHits.size()<2 || selectedHits.recalculateNPlanesHit()<=1) continue;
			if(constructLineParallelToZ) fit = SimpleFitResult{0,selectedHits.front().x, 0, selectedHits.front().y, 0,0,0,0};
			else fit=linearRegressionFit(selectedHits);

			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; cin.get(); continue;	}

			//sum x and y for calculation of centre of mass from track hits
			if(hitCluster.size()) for(auto& h : hitCluster) {
				hitsXSum[h.plane]+=h.x;
				hitsYSum[h.plane]+=h.y;
				++nHits[h.plane];
			}

			residuals=calculateResiduals(hitCluster, fit);

			residualHistograms->fill(residuals, hitsCentre);

			//sum residuals
			if(residuals.size()) for(auto& r : residuals) {
//				cout<<r.x<<" "<<r.y<<endl;
				residualXSum[r.h.plane]+=r.x;
				residualYSum[r.h.plane]+=r.y;

				auto& rotationPoint=hitsCentre[r.h.plane];
				double xc=rotationPoint.first, yc=rotationPoint.second; //x and y center
				double hx=r.h.x-xc, hy=r.h.y-yc;
				double phi=(hy*r.x-hx*r.y)/(hx*hx+hy*hy);
				double weight=hx*hx+hy*hy;
				rotationZSum[r.h.plane]+=phi*weight; rotationZWeightSum[r.h.plane]+=weight;
			}

			//sum fit slope
			slope1Sum+=fit.slope1; slope2Sum+=fit.slope2;

			fits.push_back(fit);
			if(makeTrackHistograms) { trackHistograms->fill(fit); }

			++nClusters;
		}

		if(displayEvent) {
			HoughTransformer::drawClusters(spaceHit, detector);
//			HoughTransformer::drawClusters(houghClusters, detector);
			for(auto& f : fits) f.draw(0, detector.planePosition.back());
			gPad->Update();
			auto signal=std::cin.get();
			if(signal=='q') break;
			else if(signal=='l') {
				while(!gSystem->ProcessEvents()) {
				   gSystem->Sleep(50);
				}
				break;
			} else if(signal=='w') {
				//rotate and write as animated gif!
				double phiView=40;
				for(int thetaView=0; thetaView<360; thetaView+=2 ) {
				   gPad->GetView()->RotateView(thetaView,phiView);
				   gPad->Modified();
				   gPad->Update();
				   gPad->Print( thetaView==358 ? "eventAnimation.gif++5++" : "eventAnimation.gif+5");
				}
			}
		}

	}

	std::cout<<"passed: "<<nPassed<<"/"<<nEvents<<" with "<<nClusters<<"\n";

	//Recalculate centre of mass of hits for each plane
	if(recalculateCOM) {
		cout<<"centre of chip is "<<detector.planexmax()/2<<", "<<detector.planeymax()/2<<endl;
		for(int i=0; i<detector.nPlanes; ++i) {
			hitsCentre[i]={ hitsXSum[i]/nHits[i], hitsYSum[i]/nHits[i] };
			cout<<"Hits centre of mass is: ("<<hitsCentre[i].first<<", "<<hitsCentre[i].second<<")"<<endl;
		}
		recalculateCOM=false;
	}

	for(int i=0; i<detector.nPlanes; i++) {
		cout<<"average residual is "<<residualXSum[i]<<"/"<<nHits[i]<<", "<<residualYSum[i]<<"/"<<nHits[i]<<endl;
		averageResidualFromSum[i]= { residualXSum[i]/nHits[i], residualYSum[i]/nHits[i] };
		rotationZFromSum[i]= rotationZSum[i]/rotationZWeightSum[i];
	}

	slope1FromSum=slope1Sum/nClusters;
	slope2FromSum=slope2Sum/nClusters;
	cout<<"slopes are "<<slope1FromSum<<" "<<slope2FromSum<<endl;

//	cin.get();

}

std::vector<std::pair<double, double> > trackFitter::getMeans() {
//	auto means= residualHistograms->getMeansOfPlanes();
	auto means= averageResidualFromSum;
	for(auto& m : means ) cout<<"shift ("<<m.first<<", "<<m.second<<")"<<endl;
	return means;
}

std::vector<double> trackFitter::getRotations() {
//	auto rotation= residualHistograms->getRotationOfPlanes();
	auto rotation= rotationZFromSum;
	for(auto& r : rotation ) cout<<"rotation: "<<r<<" = "<<r/M_PI*180.<<endl;
	return rotation;
}

void trackFitter::setShifts(
		const std::vector<std::pair<double, double> >& shiftsIn) {
		shifts=shiftsIn;
}

void trackFitter::setAngles(const std::vector<double>& anglesIn) {
	angles=anglesIn;
}

void trackFitter::addToShifts(
		const std::vector<std::pair<double, double> >& shiftsExtra) {
	if(shifts.empty()) return setShifts(shiftsExtra);

	auto itExtra= shiftsExtra.begin();
	for(auto& s : shifts) {
		s.first+=itExtra->first;
		s.second+=itExtra->second;
		++itExtra;
	}
}

void trackFitter::addToAngles(const std::vector<double>& anglesExtra) {
	if(angles.empty()) return setAngles(anglesExtra);

	auto itExtra=anglesExtra.begin();
	for(auto& a : angles) {
		a+=*itExtra++;
	}
}

const std::vector<std::pair<double, double> >& trackFitter::getShifts() const {
	return shifts;
}

const std::vector<double>& trackFitter::getAngles() const {
	return angles;
}

std::pair<double, double> trackFitter::getSlopes() const {
	return {slope1FromSum, slope2FromSum};
}

void trackFitter::setSlopes(std::pair<double, double> slopes) {
	houghTransform.angleOfTracksX=slopes.first;
	houghTransform.angleOfTracksY=slopes.second;
}
