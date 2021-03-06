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
#include "TCanvas.h"
#include "TView.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"

#include "PositionHit.h"
#include "Hit.h"
#include "DetectorConfiguration.h"

struct HoughTransformer {

	constexpr static int nPlanes=6; //todo: make dynamic

	HoughTransformer( double xmin, double xmax, double ymin, double ymax, int xbins, int ybins ) :
		xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), xbins(xbins), ybins(ybins) {};

	virtual ~HoughTransformer() {};

	struct HitCluster : std::list<PositionHit> { //todo: make member variable instead of inheritance (bad practice inheriting from stl classes)
		int clusterSize=0, planeHit[nPlanes] = {} ;
		void add( const PositionHit& h ) {
			push_back(h);
			++clusterSize;
			++planeHit[h.plane];
		};
		void add( const PositionHit& h, int plane ) {
			push_back(h);
			++clusterSize;
			++planeHit[plane];
		};
		void clear() {
			 std::list<PositionHit>::clear();
			 clusterSize=0;
			 for(int i=0; i<nPlanes; i++) planeHit[i]=0;
		};
		void mergeWith(HitCluster& o) {
			for(int i=0; i<nPlanes; i++) planeHit[i] += o.planeHit[i];
			clusterSize+=o.clusterSize;
			splice(begin(), o);
			o.clear();
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
		TVector3 getAveragePosition(bool rejectFlagged=false) const {
			TVector3 sum(0,0,0);
			int n=0;
			for(const auto& hit : *this) {
				if(rejectFlagged && hit.flag<0) continue;
				sum+=hit.getPosition();
				n++;
			}
			sum*=(1./n);
			return sum;
		}
		int getNHitsUnflagged() {
			int n=0;
			for(const auto&h : (*this) ) {
				if(h.flag>0) ++n;
			}
			return n;
		}
	};

	template<class T >
	static void drawClusters(const T& clusters, const DetectorConfiguration& detector);
	template<class T >
	static void drawCluster(const T& cluster, const DetectorConfiguration& detector);

	double xmin, xmax, ymin, ymax;
	int xbins, ybins;

	int minClusterSize=5, minCandidateSize=3;
	double angleOfTracksX=0., angleOfTracksY=0.;
	

	void changeStat(TH1&h) {
		gPad->Update();
		TPaveStats *st = (TPaveStats*)h.FindObject("stats");
		st->SetX1NDC(0.12); st->SetY1NDC(0.83);
		st->SetX2NDC(0.32); st->SetY2NDC(0.88);
		st->SetLineWidth(0);
		st->Draw();
		gPad->Update();
	}

	//telescope
	//vector of vectors, this means multiple planes
	virtual std::list<HitCluster> operator() ( const std::vector< std::vector<PositionHit> >& hv  ) {

		//construct grid
		std::vector< std::vector< std::unique_ptr<HitCluster> > > houghGrid( xbins );
		for(auto& v : houghGrid) v.resize(ybins);

		constexpr bool DrawHistogram=false;
		std::unique_ptr<TH2D> graphicHistogram{
				DrawHistogram ?
				new TH2D("graphicHistogram", "Histogram of hough transform;x bin;y bin", xbins,0,xbins, ybins,0,ybins):
				nullptr };
		static TCanvas* canv=new TCanvas("houghTelCanv", "Canvas with cluster histogram", 600,400);
		for( int plane=0; plane<nPlanes; plane++ ) {
			for(auto& h : hv[plane] ) {
				int binx= (h.x-xmin-angleOfTracksX*h.z)/(xmax-xmin)*xbins;
				int biny= (h.y-angleOfTracksY*h.z)/ymax*ybins;
				if( binx>=xbins ) binx=xbins-1;
				if( biny>=ybins ) biny=ybins-1;
				if( binx<0 ) binx=0;
				if( biny<0 ) biny=0;

				if(!houghGrid.at(binx).at(biny)) houghGrid.at(binx).at(biny)=std::unique_ptr<HitCluster>( new HitCluster() );
				houghGrid.at(binx).at(biny)->add(h, plane);

				if(DrawHistogram) graphicHistogram->Fill(binx,biny);
			}
		}

		//draw histogram
		if(DrawHistogram) {
			std::cout<<"drawing histogram of hough transform!"<<std::endl;
			canv->cd();
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat("e");
			graphicHistogram->Draw("colz");
			changeStat(*graphicHistogram);
			canv->SetTicks(1,1);
			char c=std::cin.get();
			if(c=='w') {gPad->Print("telescopeCluster.pdf");}
			if(c=='q') { throw graphicHistogram;} //abuse of throw mechanism
		}

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
			auto& currentCluster= *houghGrid.at(i).at(j);
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
	

	//timepix
	//just a vector, this is one plane
	std::list<HoughTransformer::HitCluster> operator() ( const std::vector<PositionHit>& hv ) {

		//construct grid
		std::vector< std::vector< std::unique_ptr<HitCluster> > > houghGrid( xbins ); //Houghgrid[z][y]
		for(auto& v : houghGrid) v.resize(ybins);

		//set drawhistogram to true to make a histogram of the hough-like transform.
		constexpr bool DrawHistogram=false;
		static TCanvas* canv=nullptr;
		if(canv) canv->Clear();
		std::unique_ptr<TH2D> graphicHistogram{
				DrawHistogram ?
				new TH2D("graphicHistogram", "Histogram of hough transform;x bin;y bin", xbins,0,xbins, ybins,0,ybins ):
				nullptr};
		for(auto& h : hv ) {
			int binx= (h.x-xmin-angleOfTracksX*h.z)/(xmax-xmin)*xbins;
			int biny= (h.y-angleOfTracksY*h.z)/ymax*ybins;
			if( binx>=xbins ) binx=xbins-1;
			if( biny>=ybins ) biny=ybins-1;
			if( binx<0 ) binx=0;
			if( biny<0 ) biny=0;

			if(!houghGrid.at(binx).at(biny)) houghGrid.at(binx).at(biny)=std::unique_ptr<HitCluster>( new HitCluster() );
			houghGrid.at(binx).at(biny)->add(h);

			if(DrawHistogram) graphicHistogram->Fill(binx,biny);
		}

		//draw histogram
		if(DrawHistogram) {
			std::cout<<"drawing histogram of hough transform!"<<std::endl;
			if(!canv) canv=new TCanvas("houghCanv", "Canvas with cluster histogram", 600,400);
			canv->cd();
			graphicHistogram->Draw("colz");
			changeStat(*graphicHistogram);
			char c=std::cin.get();
			canv->SetTicks(1,1);
			gPad->Update();
			if(c=='w') {gPad->Print("timepixCluster.pdf");}
			if(c=='q') { throw 1;} //abuse of throw mechanism
		}


		struct ClusterCandidate {
			int size, x, y;
			ClusterCandidate(int size, int x, int y) : size(size),x(x),y(y) {};
		};

		//get grid positions
		std::list< ClusterCandidate > ClusterCandidatePositions;
		for(int i=0; i<xbins; i++) {
			for(int j=0; j<ybins; j++) {
				if(houghGrid.at(i).at(j)) {
					int size=houghGrid.at(i).at(j)->clusterSize;
					if(size >= minCandidateSize) {
						ClusterCandidatePositions.emplace_back(size, i, j);
					}
				}
			}
		}
		//sort by size
		ClusterCandidatePositions.sort([](const ClusterCandidate& a, const ClusterCandidate& b) {return a.size>b.size;});

		//merge neighbouring bins, starting with largest
		std::list<HitCluster> foundClusters;
		for(auto it=ClusterCandidatePositions.rbegin();	it!=ClusterCandidatePositions.rend(); it++) {
			int i=it->x, j=it->y;

			if(!houghGrid.at(i).at(j)) continue;
			auto& currentCluster= *houghGrid.at(i).at(j);
			if(!currentCluster.clusterSize) continue; //already merged with another

			//merge with neighbours
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
	pointTree.Draw("h.z:h.y:h.x:cluster*10", "", "*colz");
	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	const double xmax=detector.planexmax(), ymax=detector.planeymax();
	axisObject->GetXaxis()->SetLimits(0,xmax);
	axisObject->GetYaxis()->SetLimits(0,ymax);
	axisObject->DrawClone();
	gPad->Update();
}

template<class T > //=std::list<HitCluster>
inline void HoughTransformer::drawCluster(const T& cluster, const DetectorConfiguration& detector) {
	TTree pointTree;
	PositionHit h(1E4,1E4,1E4, 0,0,0, 0);
	pointTree.Branch("h", &h);
	pointTree.Fill(); //fill one with zero ToT to set scale

	for(auto& iHit : cluster ) {
		h=iHit;
		if(h.flag>0)
			pointTree.Fill();
	}

	gStyle->SetMarkerStyle(20);
	gStyle->SetOptTitle(0);
	double totAxis=1.6;
	pointTree.Draw( ("h.x:h.y:h.z:TMath::Min(h.ToT*0.025, "+std::to_string(totAxis)+")").c_str() , "", "*colz"); //ToT to microseconds
	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	auto zaxis=axisObject->GetZaxis();
	zaxis->SetLimits(detector.xmin(),detector.xmax());
	zaxis->SetTitle("z-axis (drift direction) [mm]") ;
	auto yaxis=axisObject->GetYaxis();
	yaxis->SetTitle("y-axis [mm]");
	yaxis->SetLimits(detector.ymin(),detector.ymax());
	auto xaxis=axisObject->GetXaxis();
	xaxis->SetLimits(detector.zmin(),detector.zmax());
	xaxis->SetTitle("x-axis (beam direction) [mm]");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.3);
	for(auto axis : {xaxis, yaxis,zaxis} ) {
		axis->SetTitleSize(0.05);
		axis->SetLabelSize(0.05);
	}
//	axisObject->SetMaximum(totAxis);
//	axisObject->SetMinimum(0);
	axisObject->Draw("colz");
//	std::cout<<"give angles!"<<std::endl;
//	std::cin>>theta>>phi;
	gPad->SetMargin(0.1,0.175,0.15,0.1);//l r b t
	gPad->Update();


	TPaletteAxis* palette= dynamic_cast<TPaletteAxis*>(gPad->GetPrimitive("palette"));
	if(!palette) throw "could not find paletteAxis!";
	palette->SetX1NDC(0.85);
	palette->SetX2NDC(0.90);
	palette->SetY2NDC(0.74);
	palette->SetY1NDC(0.1);
	//draw TPaveText over Palette axis title
	auto paletteAxisLabel = new TPaveLabel(0.96,0.1,1,0.75, "ToT [#mus]", "NDC");
	paletteAxisLabel->SetFillColor(kWhite);
	paletteAxisLabel->SetBorderSize(0);
	paletteAxisLabel->SetTextAngle(90);
	paletteAxisLabel->SetTextSize(0.08);
	paletteAxisLabel->SetTextFont(42);
	paletteAxisLabel->SetTextAlign(kHAlignCenter+kVAlignCenter);
	paletteAxisLabel->Draw();


	double theta=70,phi=60;
	gPad->GetView()->RotateView(theta, phi);

	TLegend* legend= new TLegend( 0.6, 0.8, 0.95,0.95 );
	legend->SetName("eventDisplayLegend");
	legend->AddEntry(axisObject, "Timepix hits", "p");
	axisObject->SetLineColor(kOrange+7);
	axisObject->SetLineWidth(2);
	legend->AddEntry(axisObject, "Track fit", "l");
	legend->Draw();

	gPad->Update();

}


#endif
