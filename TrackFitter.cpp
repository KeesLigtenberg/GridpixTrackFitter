/*
 * trackFitter.cpp
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#include "TrackFitter.h"

#include "TSystem.h"
#include "TView.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "linearRegressionFit.h"
#include "makeNoisyPixelMask.h"
#include "transformHits.h"
#include "ResidualHistogrammer.h"

using namespace std;

trackFitter::trackFitter(std::string inputfile, const DetectorConfiguration& detector) :
	detector(detector),
	houghTransform(detector.planexmax(), detector.planeymax(), 50/*xbins*/, 25 /*ybins*/ ),
	residualHistograms(nullptr)
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
	residualHistograms=nullptr; //make sure to close residualHistograms before track
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

void trackFitter::fitTracks(std::string outputfilename) {
	residualHistograms=unique_ptr<ResidualHistogrammer>(new ResidualHistogrammer(outputfilename, detector));
	if(makeTrackHistograms) trackHistograms=unique_ptr<TrackHistogrammer>(new TrackHistogrammer(detector) );

	//loop over all entries
	const long long nEvents=hitTable->GetEntriesFast(); //std::min( (long long) 2,);
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%1000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<std::endl;

		//get entry
		hitTable->GetEntry(iEvent);

		std::vector<std::vector<PositionHit> > spaceHit;

		//apply mask
		if(!mask.empty()) {
			auto itmask=mask.begin();
			for(unsigned plane=0; plane<mimosaHit->size(); plane++ ) {
				auto maskedHit = applyPixelMask( *itmask++ , mimosaHit->at(plane) ) ;
				//convert hits to positions
				spaceHit.push_back( convertHits(maskedHit, detector.planePosition[plane], detector.pixelsize, detector.pixelsize, plane) );
			}
		}

		if( !passEvent(spaceHit) ) continue;
		++nPassed;

		//apply translation and rotation
		if(!shifts.empty()) spaceHit=translateHits(spaceHit,shifts);
		if(!angles.empty()) spaceHit=rotateHits(spaceHit, angles);

		//hough transform
		auto houghClusters = houghTransform(spaceHit);

		//require at least 5/6 planes to be hit
		const int nMinPlanesHit=5;
		houghClusters.remove_if([](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
		if(!houghClusters.size()) {
			continue;
		}

//		std::cout<<"found "<<houghClusters.size()<<" clusters of size ";
//		for(auto& cl : houghClusters) {
//			std::cout<<cl.clusterSize<<", ";
//		}
//		std::cout<<std::endl;
//		HoughTransformer::drawClusters(houghClusters, mimosa);


		//first fit on outer planes and translate inner hits
		std::vector<SimpleFitResult> fits;
		for(auto& hitCluster : houghClusters) {

			//remove all hits except first and last one
//			auto hitsOnOuterPlanes=hitCluster;
//			const int firstPlane=0, lastPlane=5;
//			hitsOnOuterPlanes.remove_if([](PositionHit& hit){return hit.plane != firstPlane && hit.plane != lastPlane;});

			//fit track
			auto fit=linearRegressionFit(hitCluster);
//			fit.draw(0, mimosa.planePosition.back() );

			//remove outliers
			auto residuals=calculateResiduals(hitCluster, fit);
			const double maxResidual=0.2;
			hitCluster=cutOnResiduals(hitCluster, residuals, maxResidual);

			++nClusters;
			//refit on outer planes again
			fit=linearRegressionFit(hitCluster);
			residuals=calculateResiduals(hitCluster, fit);
			residualHistograms->fill(residuals);
			fits.push_back(fit);

			if(makeTrackHistograms) { trackHistograms->fill(fit); }
		}

		if(displayEvent) {
//			HoughTransformer::drawClusters(spaceHit, detector);
			HoughTransformer::drawClusters(houghClusters, detector);
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


}

std::vector<std::pair<double, double> > trackFitter::getMeans() {
	auto means= residualHistograms->getMeansOfPlanes();
	for(auto& m : means ) cout<<"shift ("<<m.first<<", "<<m.second<<")"<<endl;
	return means;
}

std::vector<double> trackFitter::getRotations() {
	auto rotation= residualHistograms->getRotationOfPlanes();
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
