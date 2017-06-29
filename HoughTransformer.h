#ifndef GETHOUGHTRANSF_H
#define GETHOUGHTRANSF_H

#include <list>
#include <functional>
#include <memory>
#include <iostream>

#include "TStyle.h"
#include "TTree.h"
#include "TH2.h"
#include "TPad.h"

#include "PositionHit.h"
#include "Hit.h"
#include "DetectorConfiguration.h"

struct HoughTransformer {

	constexpr static int nPlanes=6;

	HoughTransformer(double xmax, double ymax, int xbins, int ybins) :
		xmax(xmax), ymax(ymax), xbins(xbins), ybins(ybins) {};

	struct HitCluster : std::list<PositionHit> {
		int clusterSize=0, planeHit[nPlanes] = {0} ;
		void add( PositionHit h, int plane ) {
			push_back(h);
			++clusterSize;
			++planeHit[plane];
		};
		void clear() {
			 std::list<PositionHit>::clear();
			 clusterSize=0;
		};
		void mergeWith(HitCluster& o) {
			for(int i=0; i<nPlanes; i++) planeHit[i] += o.planeHit[i];
			clusterSize+=o.clusterSize;
			o.clusterSize=0;
			splice(begin(), o);
		}
		int getNPlanesHit() const {
			int n=0;
			for(int i=0; i<nPlanes; i++) if(planeHit[i]) ++n;
			return n;
		}
		int recalculateNPlanesHit() {
			for(int i=0; i<nPlanes;i++) planeHit[i]=0;
			for( PositionHit& h : (*this) ) {
				++planeHit[h.plane];
			}
			return getNPlanesHit();
		}
	};

	template<class T >
	static void drawClusters(const T& clusters, const DetectorConfiguration& detector);

	double xmax, ymax;
	int xbins, ybins;

	
	std::list<HitCluster> operator() ( const std::vector< std::vector<PositionHit> >& hv, int minClusterSize=5, int minCandidateSize=3  ) {

		//construct grid
		std::vector< std::vector< std::unique_ptr<HitCluster> > > houghGrid( xbins );
		for(auto& v : houghGrid) v.resize(ybins);

//		TH2D graphicHistogram("graphicHistogram", "Histogram of hough transform", xbins,0,xbins, ybins,0,ybins );
		for( int plane=0; plane<nPlanes; plane++ ) {
			for(auto& h : hv[plane] ) {
				int binx= h.x/xmax*xbins;
				int biny= h.y/ymax*ybins;
				if( binx>=xbins ) binx=xbins-1;
				if( biny>=ybins ) biny=ybins-1;
				if( binx<0 ) binx=0;
				if( biny<0 ) biny=0;

				if(!houghGrid.at(binx).at(biny)) houghGrid.at(binx).at(biny)=std::unique_ptr<HitCluster>(new HitCluster());
				houghGrid.at(binx).at(biny)->add(h, plane);

//				graphicHistogram.Fill(binx,biny);
			}
		}

		//draw histogram
//		std::cout<<"drawing histogram of hough transform!"<<std::endl;
//		graphicHistogram.Draw("colz");
//		gPad->Update();
//		if(std::cin.get()=='q') { throw graphicHistogram; //serious abuse of language
//		}

		//get grid positions and sort by size
		std::list< std::tuple<int, int, int> > gridPositions;
		for(int i=0; i<xbins; i++) {
			for(int j=0; j<ybins; j++) {
				if(houghGrid.at(i).at(j)) {
					int size=houghGrid.at(i).at(j)->clusterSize;
					if(size >= minCandidateSize) {
						gridPositions.emplace_back(size, i, j);
					}
				}
			}
		}
		gridPositions.sort();

		//merge neighbouring bins, starting with largest
		std::list<HitCluster> foundClusters;
		for(auto it=gridPositions.rbegin();	it!=gridPositions.rend(); it++) {
			int size, i, j;
			std::tie(size,i,j)=*it;

			if(!houghGrid.at(i).at(j)) continue;
			auto& currentCluster=*houghGrid.at(i).at(j);
			if(!currentCluster.clusterSize) continue;

			for(int dx=-1; dx<=1; dx++) for(int dy=-1; dy<=1; dy++) {
				if(i+dx<0 or i+dx >= xbins or j+dy<0 or j+dy >= ybins or (!dx and !dy) ) continue; //outside grid
				if(! houghGrid.at(i+dx).at(j+dy) ) continue; //no hits in bin
				currentCluster.mergeWith( *houghGrid.at(i+dx).at(j+dy) );
			}

			if(currentCluster.clusterSize >= minClusterSize) {
				foundClusters.push_back( std::move(currentCluster) );
			}
			houghGrid.at(i).at(j).reset(nullptr);
		}

		return foundClusters;
	}
	
};

template<class T > //=std::list<HitCluster>
inline void HoughTransformer::drawClusters(const T& clusters, const DetectorConfiguration& detector) {
	TTree pointTree;
	PositionHit h(0,0,0);
	int cluster=1;
	pointTree.Branch("h", &h);
	pointTree.Branch("cluster", &cluster);

	for(auto& iClus : clusters) {
		for(auto& iHit : iClus ) {
			h=iHit;
			pointTree.Fill();
		}
		++cluster;
	}

	gStyle->SetMarkerStyle(20);
	pointTree.Draw("h.z:h.y:h.x:cluster*10", "", "*");
	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	const double xmax=detector.planexmax(), ymax=detector.planeymax();
	axisObject->GetXaxis()->SetLimits(0,xmax);
	axisObject->GetYaxis()->SetLimits(0,ymax);
	axisObject->DrawClone();
	gPad->Update();
}




#endif
