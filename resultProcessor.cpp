#include "resultProcessor.h"

#include <deque>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile2D.h"

#include "testBeamSetup.h"

resultProcessor::resultProcessor(TTree *tree) :
	fChain(0),
	alignment("alignment.dat")
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("fitResults.root");
      if (!f || !f->IsOpen()) {
    	  f=TFile::Open("fitResults.root", "READ");
      }
      if (!f || !f->IsOpen()) {
         std::cerr<<"could not find file!"<<std::endl;
         throw 1;
      }
      f->GetObject("fitResults",tree);

   }
   Init(tree);
}

resultProcessor::~resultProcessor()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t resultProcessor::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   auto nb=fChain->GetEntry(entry);

   //calculate derived quantities
   timepixFrameFits=transformFitsToTimepixFrame(*timepixFits, alignment.relativeAlignment);

   return nb;
}
Long64_t resultProcessor::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void resultProcessor::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   timepixClusterSize = 0;
   timepixHits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
//   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("telescopeFits", &telescopeFits, &b_telescopeFits);
   fChain->SetBranchAddress("timepixFits", &timepixFits, &b_timepixFits);
   fChain->SetBranchAddress("ntimepixHits", &ntimepixHits, &b_ntimepixHits);
   fChain->SetBranchAddress("ntelescopeHits", &ntelescopeHits, &b_ntelescopeHits);
   fChain->SetBranchAddress("timepixClusterSize", &timepixClusterSize, &b_timepixClusterSize);
   fChain->SetBranchAddress("timepixHits", &timepixHits, &b_timepixHits);
   fChain->SetBranchAddress("dxz", &dxz, &b_dxz);
   fChain->SetBranchAddress("dyz", &dyz, &b_dyz);
   fChain->SetBranchAddress("nresiduals", &nresiduals, &b_nresiduals);
   fChain->SetBranchAddress("nfitted", &nfitted, &b_nfitted);
   Notify();
}

Bool_t resultProcessor::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void resultProcessor::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t resultProcessor::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   if(timepixClusterSize->front()/double(ntimepixHits)<0.75) return -1;
   if(timepixClusterSize->front()<30) return -1;
//		   or timepixClusterSize->front()>300

   return 1;
}

bool isInsideArea(int x, int y) {
	if(x<=16 or x>=240) {return false;}
	if(y< 80-50./225*x) return false;
	return true;
}

TProfile2D* removeBinsWithFewerEntries(TProfile2D* h, int minEntries){
	for(int x=1; x<=h->GetNbinsX(); x++) {
		for(int y=1; y<=h->GetNbinsY(); y++) {
			int bin=h->GetBin(x,y);
			if( h->GetBinEntries(bin) < minEntries
//					or !isInsideArea( h->GetXaxis()->GetBinCenter(x), h->GetYaxis()->GetBinCenter(y) )
				) {
				h->SetBinEntries(bin, 0);
			}
		}
	}
	return h;
}


//do frequency on an unweighted histogram, i.e. all entreis should have same weight
void getFrequencyHistogram(TProfile2D* original, double systematicError, double min=-0.1, double max=0.1, int nBins=40, double entryWeight=1.0) {
	TH1D* frequencyHist=new TH1D(
		(original->GetName()+std::string("freq")).c_str(),
		(original->GetTitle()+std::string("frequency; frequency;")+original->GetXaxis()->GetTitle() ).c_str(),
		nBins, min, max);
	TH1D* frequencyPull=new TH1D(
		(original->GetName()+std::string("freqPull")).c_str(),
		(original->GetTitle()+std::string("frequencyPull; frequency;")+original->GetXaxis()->GetTitle() ).c_str(),
		100, -5, 5);
	for(int i=1;i<=original->GetNbinsX();i++) {
		for(int j=1; j<=original->GetNbinsY();j++) {
			int bin=original->GetBin(i,j);
			if(!isInsideArea( original->GetXaxis()->GetBinCenter(i), original->GetYaxis()->GetBinCenter(j) )) continue;
			if( original->GetBinEntries(bin) > 0 ) {
				double binContent= original->GetBinContent(i,j); //possible rounding error here
				double binError=original->GetBinError(i,j);
				double error=sqrt(binError*binError+systematicError*systematicError);
				frequencyHist->Fill(binContent);
				frequencyPull->Fill(binContent/error);
			}
		}
	}
}

std::vector<int> rebinHits(const std::vector<int>& hits, int startbin, int endbin, int nBinsPerNewBin) {
	if(!(endbin-startbin)%nBinsPerNewBin) {
		std::cerr<<"bins should exactly match range!\n"; throw 1;
	}
	std::vector<int> newBins((endbin-startbin)/nBinsPerNewBin);
	for(int i = startbin; i<endbin; i++ ) {
		newBins[(i-startbin)/nBinsPerNewBin]+=hits[i];
	}
	return newBins;
}

double truncatedMean(std::vector<int>& binned, double fraction=0.7) {
	  std::sort(binned.begin(), binned.end());
	  int sumFirstN=binned.size()*fraction;
	  double sum=std::accumulate(binned.begin(), binned.begin() + sumFirstN , 0);
	  double mean=sum/sumFirstN;
	  return mean;
}
template <class RAcontainer=std::vector<int> >
double truncatedSum(RAcontainer& binned, double fraction=0.7) {
	  std::sort(binned.begin(), binned.end());
	  int sumFirstN=binned.size()*fraction;
	  double sum=std::accumulate(binned.begin(), binned.begin() + sumFirstN , 0);
	  return sum;
}

bool isZero(double x,double threshold=1E-30) { return fabs(x)<threshold; }
struct crosstalkCalculator {
	  std::array<std::array<double, 256>, 256> ToTmap{{}};
	  void Fill(int col, int row, double ToT) {
		  ToTmap[col][row]=ToT;
	  }
	  bool isFilled(int i, int j) {
		  return not isZero(ToTmap[i][j]);
	  }
	  int getNNeighbours(int i, int j) {
		  int nNeighbours=0;
		  for(int a=-1; a<=1; a++) {
			  for(int b=-1; b<=1; b++) {
				  if(a==0 and b==0) continue;
				  if(i+a<0 or j+b<0 or i+a>=256 or j+b>=256) continue; //out of grid
				  if( isFilled(i+a,j+b) ) ++nNeighbours;
			  }
		  }
		  return nNeighbours;
	  }
	  bool isIsolated(int i, int j) {
		  return not getNNeighbours(i,j);
	  }
	  bool isPair(int i, int j) { //is part of pair
		  if(not isFilled(i,j)) return false;
		  if(getNNeighbours(i,j)!=1) return false;
		  //has one neighbour, find it and check if also one neighbour.
		  for(int a=-1; a<=1; a++) {
			  for(int b=-1; b<=1; b++) {
				  if(a==0 and b==0) continue;
				  if(i+a<0 or j+b<0 or i+a>=256 or j+b>=256) continue; //out of grid
				  if( isFilled(i+a,j+b) ) {//neighbour of hit under inspection
					  if(getNNeighbours(i+a,j+b)!=1) return false;
					  if( std::abs(a)+std::abs(b)!=1  ) return false; //should be a direct neighbour! no diagonals
				  }
			  }
		  }
		  return true;
	  }
	  std::vector<double> getIsolatedToTs() {
		  std::vector<double> ToTs;
		  for(int i=0; i<256; i++) {
			  for(int j=0; j<256; j++) {
				  if( isFilled(i,j) and isIsolated(i,j) ) {
					  	ToTs.push_back(ToTmap[i][j]);
				  }
			  }
		  }
//		  std::cout<<"isolated ToTs: "<<ToTs.size()<<"\n";
		  return ToTs;
	  }
	  std::vector<double> getPairToTs() {
		  std::vector<double> ToTs;
		  for(int i=0; i<256; i++) {
			  for(int j=0; j<256; j++) {
				  if( isPair(i,j) ) {
					  ToTs.push_back(ToTmap[i][j]);
				  }
			  }
		  }
		  return ToTs;
	  }
};

#pragma link C++ class std::vector<int>+;

void resultProcessor::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   TFile output("histograms.root", "RECREATE");

   TH2D hitmap("hitmap", "map of hits in tracks", 256,0,256,256,0,256 );
   TProfile ToTByCol("ToTByCol", "ToT by column", 256,0,256, 0, 4);
   const int nbins=64;
   TProfile2D deformationsyExp("deformationsyExp", "profile of y residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsxExp("deformationsxExp", "profile of x residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsy("deformationsy", "profile of y residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsx("deformationsx", "profile of x residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TH2D diffusionx("diffusionx", "x residuals as a function of drift distance", 50,4,20,40,-2,2);
   TH2D diffusiony("diffusiony", "y residuals as a function of drift distance", 50,4,20,40,-2,2);

   TH1D trackLength("trackLength", "length of track in tpc", 40,13,16);
   TH1D truncatedSumHits("truncatedSumHits", "mean per n bins", 100, 0, 100);
   TH1D sumHits("sumHits", "mean per n bins", 200, 0, 200);
   TH1D chargePerBin("chargePerBin", "chargePerBin", 400,0,400);
   TH1D aggravatedMeanHits("aggravatedHits", "mean per n bins", 400, 0, 400);

   TH1D isolatedToT("isolatedToT","ToT of isolated hits", 40,0,2);
   TH1D pairToT("pairToT","ToT of isolated hit-pairs", 40,0,2);

   std::vector<int> hitsAlongTrack(512);
   std::deque<int> aggravatedHits;
   TTree binnedHitsTree("binnedHitsTree", "tree with hits per track");
   binnedHitsTree.Branch("hits", &hitsAlongTrack );

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	  if(!(jentry%10000)) std::cout<<"entry "<<jentry<<"/"<<nentries<<"\n";
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      auto nb = GetEntry(jentry);  nbytes += nb;
       if (Cut(ientry) < 0) continue;

       trackLength.Fill(timepixFrameFits.at(0).getTrackLength( timePixChip.zmin(), timePixChip.zmax() ));

 	  std::fill(hitsAlongTrack.begin(), hitsAlongTrack.end(), 0);
 	  crosstalkCalculator crosstalk;


      for(auto& h : timepixHits->front() ) {
    	  if(h.ToT*0.025<0.15 or h.flag<0) continue;

    	  hitmap.Fill( h.col, h.row );
    	  crosstalk.Fill(h.col,h.row,h.ToT*0.025);

    	  double hryp=h.ry/cos(timepixFits->front().YZ.slope);
    	  deformationsxExp.Fill( h.col-h.rz/.055, h.row+h.ry/.055, h.rx );
    	  deformationsyExp.Fill( h.col-h.rz/.055, h.row+h.ry/.055, h.ry/cos(timepixFits->front().YZ.slope) ); //todo: check where XZ slope enters
    	  deformationsx.Fill( h.col, h.row+hryp/0.055, h.rx );
    	  deformationsy.Fill( h.col, h.row+hryp/0.055, h.ry/cos(timepixFits->front().YZ.slope) ); //todo: check where XZ slope enters

    	  diffusiony.Fill(h.x-h.rx, h.ry);
    	  diffusionx.Fill(h.x-h.rx, h.rx);

    	  ToTByCol.Fill(h.col, h.ToT*0.025);

    	  //calculate xy
    	  double xInTimepixFrame=0.055*h.col;
    	  double yInTimepixFrame=0.055*h.row;

    	  {
    		  //check with algebra
//        	  auto pos=h.createPositionHit();
//        	  pos=transformHitToTimepixFrame(pos, alignment.relativeAlignment);
//    		  auto& f=timepixFrameFits.at(0);
//    		  double z0=timePixChip.zmin(), y0=f.yAt(z0);
//    		  double l=( (pos.z-z0)+(pos.y-y0)*f.YZ.slope )/sqrt(1+pow(f.YZ.slope,2));
//    		  double bin=l/(sqrt(1+pow(f.YZ.slope,2))*.055)-0.5;
//    		  std::cout<<bin<<" - "<<h.col-h.rz/0.055<<"    ("<<bin-(h.col-h.rz/0.055)<<")\n";

    		  ++hitsAlongTrack[2*(h.col-h.rz/.055)];
    	  }
      }

      //sum hits for dEdx: fill corresponding tree and make histograms
      binnedHitsTree.Fill();
      {
    	  //view bins
//    	  TH1D hitsAlongTrackHist("hitsAlongTrackHist", "hits along track", hitsAlongTrack.size(), 0, hitsAlongTrack.size() );
//    	  std::vector<double> asDouble(hitsAlongTrack.begin(), hitsAlongTrack.end());
//    	  hitsAlongTrackHist.SetContent(asDouble.data());
//    	  hitsAlongTrackHist.Draw();
//    	  gPad->Update();

		  auto rebinnedHits=rebinHits(hitsAlongTrack, 32, 472, 40);
		  for(const auto& i : rebinnedHits) {
			  chargePerBin.Fill(i);
		  }
//		  auto rebinnedHits=rebinHits(hitsAlongTrack, 16, 236, 10);
    	  double mean=truncatedSum(rebinnedHits,0.7);
		  truncatedSumHits.Fill(mean);
		  sumHits.Fill( truncatedSum(rebinnedHits,1.) );

		  aggravatedHits.insert(aggravatedHits.end(), rebinnedHits.begin(), rebinnedHits.end() );
		  if(aggravatedHits.size() >= 110 ) { //10 events with 11 bins
//		  if(aggravatedHits.size() >= 220 ) { //10 events with 22 bins
			  aggravatedMeanHits.Fill( truncatedSum(aggravatedHits, 0.7) );
			  aggravatedHits.clear();
		  }

//		  std::cout<<"total mean "<<truncatedMean(rebinnedHitsAlongTrack,1.)<<" truncated "<<mean<<"\n";
//    	  if(std::cin.get()=='q') break;;
      }

      //check pair ToT
      {
    	  auto isolatedToTs=crosstalk.getIsolatedToTs();
    	  for(auto& t : isolatedToTs) isolatedToT.Fill(t);
    	  auto pairToTs=crosstalk.getPairToTs();
    	  for(auto& t : pairToTs) pairToT.Fill(t);
      }

   }
   gStyle->SetPalette(1);


   std::vector<double> errorVector={0.005, .020, 0.005, .020};
   auto error=errorVector.begin();
   for(auto d : {&deformationsy, &deformationsx, &deformationsyExp, &deformationsxExp} ) {
	   removeBinsWithFewerEntries(d, 100);
	   getFrequencyHistogram(d, *error++);
	   d->SetMinimum(-0.1);
	   d->SetMaximum(0.1);
   }

   output.Write();
}
