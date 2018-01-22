#include "resultProcessor.h"

#include <deque>
#include <utility>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile2D.h"
#include "TRandom.h"

#include "testBeamSetup.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "goodArea.h"
//#include "/user/cligtenb/rootmacros/makeGraphsFromFits.h"

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

namespace {

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

	TProfile2D* setMinMax(TProfile2D* h, double min, double max){
		for(int x=1; x<=h->GetNbinsX(); x++) {
			for(int y=1; y<=h->GetNbinsY(); y++) {
				int bin=h->GetBin(x,y);
				double bc=h->GetBinContent(bin);
				if( bc < min ) { h->SetBinContent(bin, min); h->SetBinEntries(bin,1);}
				else if(bc > max) { h->SetBinContent(bin,max); h->SetBinEntries(bin,1);}
			}
		}
		h->SetMinimum(min);
		h->SetMaximum(max);
		return h;
	}


	//do frequency on an unweighted histogram, i.e. all entreis should have same weight
	void getFrequencyHistogram(TProfile2D* original, double systematicError, double min=-0.2, double max=0.2, int nBins=80, double entryWeight=1.0) {
		TH1D* frequencyHist=new TH1D(
			(original->GetName()+std::string("freq")).c_str(),
			(original->GetTitle()+std::string(" frequency;")+original->GetZaxis()->GetTitle()+";" ).c_str(),
			nBins, min, max);
		TH1D* frequencyHistAll=new TH1D(
			(original->GetName()+std::string("freqAll")).c_str(),
			(original->GetTitle()+std::string(" frequency (with hits outside area);")+original->GetZaxis()->GetTitle()+";" ).c_str(),
			nBins, min, max);
//		TH1D* frequencyPull=new TH1D(
//			(original->GetName()+std::string("freqPull")).c_str(),
//			(original->GetTitle()+std::string(" frequencyPull;")+original->GetZaxis()->GetTitle()+";Frequency" ).c_str(),
//			100, -5, 5);
		for(int i=1;i<=original->GetNbinsX();i++) {
			for(int j=1; j<=original->GetNbinsY();j++) {
				int bin=original->GetBin(i,j);
				if( original->GetBinEntries(bin) <= 0 ) continue;
				double binContent= original->GetBinContent(i,j); //possible rounding error here
				double binError=original->GetBinError(i,j);
				double error=sqrt(binError*binError+systematicError*systematicError);
				frequencyHistAll->Fill(binContent);
				if(!goodArea::isInsideArea( original->GetXaxis()->GetBinCenter(i), original->GetYaxis()->GetBinCenter(j) )) continue;
				frequencyHist->Fill(binContent);
//				frequencyPull->Fill(binContent/error);

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
		  std::array<std::array<char, 256>, 256> flagmap{{}};
		  void Fill(int col, int row, double ToT) {
			  ToTmap[col][row]=ToT;
			  flagmap[col][row]=1;
		  }
		  bool isFilled(int i, int j) {
			  auto& t=flagmap[i][j];
			  return t>0;
		  }
		  using Pattern=std::vector<std::vector<char>>; //not bool, because vector<bool> is not vector<char>
		  Pattern transpose(const Pattern& pat) {
			  int csize=pat.size(), rsize=pat.at(0).size();
			  Pattern n(rsize, std::vector<char>(csize) );
			  for(int a=0; a<csize; a++) {
				  for(int b=0; b<rsize; b++) {
					  n[b][a]=pat[a][b];
				  }
			  }
			  return n;
		  }
		  //match pattern at i,j
		  bool matchPattern(int i, int j, Pattern pat) {
			  int csize=pat.size(), rsize=pat.at(0).size();
			  if(i+csize>=256 or j+rsize>=256) return false;
			  for(int a=0; a<csize; a++) {
				  for(int b=0; b<rsize; b++) {
					  bool isf=isFilled(i+a, j+b), isp=pat[a][b];
					  if( (isf && !isp) or (!isf && isp) ) {
						  return false;
					  }
				  }
			  }
			  return true;
		  }
		  //assumes pattern was already matched
		  std::vector<double> getToTsFromPattern(int i, int j, Pattern pat) {
			  std::vector<double> ToTs;
			  int csize=pat.size(), rsize=pat.at(0).size();
			  if(i+csize>=256 or j+rsize>=256) return ToTs;
			  for(int a=0; a<csize; a++) {
				  for(int b=0; b<rsize; b++) {
					  if(pat[a][b]) {
						  auto& t=ToTmap[i+a][j+b];
						  ToTs.push_back( t );
						  flagmap[1+a][j+b]=-1; //flag as visited by setting negative
					  }
				  }
			  }
			  return ToTs;
		  }
		  const std::array<Pattern,2> isolated = {{
				  { {0,0,0}, {0,1,0}, {0,0,0}, {0,0,0} },
				  { {0,0,0}, {0,0,0}, {0,1,0}, {0,0,0} },
		  }};
		  std::vector<double> getIsolatedToTs() {
			  std::vector<double> ToTs;
			  for(int i=0; i<256; i++) {
				  for(int j=0; j<256; j++) {
					  for(const auto& p : isolated) { //loop over isolated patterns
						  if( matchPattern(i,j,p) ) {
	//						  std::cout<<"pattern matched ";
							  auto ToT=getToTsFromPattern(i,j,p);
							  if(ToT.size()) ToTs.push_back(ToT.at(0)); //should be just one
	//						  else std::cout<<" but no TOT!!\n";
							  break;
						  }
						  auto pT=transpose(p);
						  if( matchPattern(i,j,pT) ) {
							  auto ToT=getToTsFromPattern(i,j,pT);
							  if(ToT.size()) ToTs.push_back(ToT.at(0)); //should be just one
							  break;
						  }
					  }
				  }
			  }
	//		  std::cout<<"isolated ToTs: "<<ToTs.size()<<"\n";
			  return ToTs;
		  }
		  const Pattern pair = {{
				  {0,0,0}, {0,1,0}, {0,1,0}, {0,0,0}
		  }};
//		  const Pattern controlPair = {{
//				  std::vector<char>{1},std::vector<char>{0},std::vector<char>{1}
//		  }};
		  const Pattern controlPair = {{
				  {0,0,0}, {0,1,0}, {0,0,0}, {0,1,0}
		  }};
		  std::vector<std::pair<double,double>> getPairToTs(const Pattern& pattern) {
			  std::vector<std::pair<double,double>> ToTs;
			  for(int i=0; i<256; i++) {
				  for(int j=0; j<256; j++) {
					  if( matchPattern(i,j,pattern) ) {
						  auto ToT=getToTsFromPattern(i,j,pattern);
						  if(ToT.size()==2) ToTs.push_back( std::minmax(ToT.at(0), ToT.at(1)) );
						  else if(ToT.size()!=0) std::cerr<<"expected 2 values!\n"; //should be two!
						  break;
					  }
					  auto pTransp=transpose(pattern);
					  if( matchPattern(i,j,pTransp) ) {
						  auto ToT=getToTsFromPattern(i,j,pTransp);
						  if(ToT.size()==2) ToTs.push_back( std::minmax(ToT.at(0), ToT.at(1)) );
						  else if(ToT.size()!=0) std::cerr<<"expected 2 values, size was "<<ToT.size()<<"\n";//should be two!
						  break;
					  }
				  }
			  }
			  return ToTs;
		  }
	};

}

#pragma link C++ class std::vector<int>+;

void resultProcessor::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   TFile output("histograms.root", "RECREATE");

   TH2D hitmap("hitmap", "map of hits in tracks", 256,0,256,256,0,256 );
   TProfile ToTByCol("ToTByCol", "ToT by column", 256,0,256, 0, 4);
   const int nbins=256;
   TProfile2D deformationsyExp("deformationsyExp", "profile of y residuals;Column;Row;y-residual [mm]", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsxExp("deformationsxExp", "profile of z residuals;Column;Row;z-residual [mm]", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsy("deformationsy", "profile of y residuals;Column;Row;y-residual [mm]", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsx("deformationsx", "profile of z residuals;Column;Row;z-residual [mm]", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TH2D diffusionx("diffusionx", "z residuals as a function of drift distance;Drift distance [mm];z-residual [mm]", 50,4,24,40,-2,2);
   TH2D diffusiony("diffusiony", "y residuals as a function of drift distance;Drift distance [mm];y-residual [mm]", 50,4,24,40,-2,2);

   TH2D timewalk("timewalk", "x residual by ToT;ToT [#mus]; x-residual [mm]", 100,0,2.5, 200,-5,5);
   TH2D timewalkCorrected("timewalkCorrected", "x residual by ToT;ToT [#mus]; x-residual [mm]", 100,0,2.5, 200,-5,5);

   TH1D trackLength("trackLength", "length of track in tpc", 40,13,16);
   TH1D truncatedSumHits("truncatedSumHits", "mean per n bins", 100, 0, 100);
   TH1D sumHits("sumHits", "mean per n bins", 200, 0, 200);
   TH1D chargePerBin("chargePerBin", "chargePerBin", 400,0,400);
   TH1D aggravatedMeanHits("aggravatedHits", "mean per n bins", 400, 0, 400);

   TH1D isolatedToT("isolatedToT","ToT of isolated hits;ToT [#mu s];Entries", 40,0,2);
   TH1D pairToTMin("pairToTMin","ToT of isolated hit-pairs;ToT [#mu s];Entries", 40,0,2);
   TH1D pairToTMax("pairToTMax","ToT of isolated hit-pairs;ToT [#mu s];Entries", 40,0,2);
   TH1D pairZDiff("pairZDiff", "Z-diff of isolated hit-pairs;#Delta z [mm];Entries", 30,0,3.515625);
   TH1D pairToTMinControl("pairToTMinControl","ToT of isolated hits;ToT [#mu s];Entries", 40,0,2);
   TH1D pairToTMaxControl("pairToTMaxControl","ToT of isolated hits;ToT [#mu s];Entries", 40,0,2);
   TH1D pairZDiffControl("pairZDiffControl", "Z-diff of isolated hit-pairs;#Delta z [mm];Entries", 30,0,3.515625);

   std::vector<int> hitsAlongTrack(512);
   std::deque<int> aggravatedHits;
   TTree binnedHitsTree("binnedHitsTree", "tree with hits per track");
   binnedHitsTree.Branch("hits", &hitsAlongTrack );

   TH1* hitHist=getHistFromTree(fChain,"timepixClusterSize", "", "nhitsHist");
   double actualNumberOfHits=getMeanFromSimpleGausFit(*hitHist);
   double targetNumberOfHits=69.541; //from 330 V run333
   constexpr bool dropHits=false;
   if(dropHits) std::cout<<"actual number of hits :"<<actualNumberOfHits<<", hits will be dropped to reach target "<<targetNumberOfHits<<"hits \n";

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0;
		   jentry<std::min(nentries,10000LL);//nentries; //
		   jentry++) {
	  if(!(jentry%100))
		  std::cout<<"entry "<<jentry<<"/"<<nentries<<"\n";
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      auto nb = GetEntry(jentry);  nbytes += nb;
       if (Cut(ientry) < 0) continue;

       trackLength.Fill(timepixFrameFits.at(0).getTrackLength( timePixChip.zmin(), timePixChip.zmax() ));

 	  std::fill(hitsAlongTrack.begin(), hitsAlongTrack.end(), 0);
 	  crosstalkCalculator crosstalk, crosstalkZ;


      for(auto& h : timepixHits->front() ) {
    	  crosstalk.Fill(h.col,h.row,h.ToT*0.025);
    	  double dxTW=alignment.timeWalkCorrection.getCorrection(h.ToT);//time walk correction
    	  if(h.ToT*0.025>0.15) crosstalkZ.Fill(h.col,h.row,h.x+dxTW);

    	  if(h.flag==-3 or h.flag>=-1) {
        	  double dxTW=alignment.timeWalkCorrection.getCorrection(h.ToT);//time walk correction
    		  timewalk.Fill(h.ToT*0.025,h.rx+dxTW);
    		  timewalkCorrected.Fill(h.ToT*0.025,h.rx);
    	  }

    	  if(h.ToT*0.025<0.15 || fabs(h.rx)>2 ) continue;

    	  double hryp=-h.ry/cos(timepixFits->front().YZ.slope);
    	  diffusiony.Fill(h.x-h.rx, hryp);
    	  diffusionx.Fill(h.x-h.rx, h.rx);

    	  hitmap.Fill( h.col, h.row );

    	  if(h.flag<0) continue;

    	  //have a chance to drop hits
    	  if(dropHits and gRandom->Rndm() > targetNumberOfHits/actualNumberOfHits) continue;


    	  deformationsxExp.Fill( h.col-h.rz/.055, h.row+h.ry/.055, h.rx );
    	  deformationsyExp.Fill( h.col-h.rz/.055, h.row+h.ry/.055, hryp ); //todo: check where XZ slope enters
    	  deformationsx.Fill( h.col, h.row, h.rx );
    	  deformationsy.Fill( h.col, h.row, hryp ); //todo: check where XZ slope enters


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
      bool dodEdX=false;
      if(dodEdX){
          binnedHitsTree.Fill();
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
      constexpr bool doTopologicalMatching=true;
      if(doTopologicalMatching){
    	  auto isolatedToTs=crosstalk.getIsolatedToTs();
    	  crosstalkZ.getIsolatedToTs();
    	  for(auto& t : isolatedToTs) isolatedToT.Fill(t);
    	  auto pairToTs=crosstalk.getPairToTs(crosstalk.pair);
    	  for(const auto& t : pairToTs) {
    		  pairToTMin.Fill(t.first);
    		  pairToTMax.Fill(t.second);
    	  }
    	  auto pairZs=crosstalkZ.getPairToTs(crosstalkZ.pair);
//    	  if(pairToTs.size()!=pairZs.size()) {std::cerr<<"ToT vec "<<pairToTs.size()<<" is different length than Z vec "<<pairZs.size()<<" in "<<jentry<<"\n";}
    	  for(const auto& z : pairZs) {
    		  pairZDiff.Fill(z.second-z.first);
    	  }
    	  auto pairControlToTs=crosstalk.getPairToTs(crosstalk.controlPair);
    	  for(auto& t : pairControlToTs) {
    		  pairToTMinControl.Fill(t.first);
    		  pairToTMaxControl.Fill(t.second);
    	  }
    	  auto pairControlZs=crosstalkZ.getPairToTs(crosstalkZ.controlPair);
    	  for(const auto& z : pairControlZs) {
    		  pairZDiffControl.Fill( z.second-z.first );
    	  }
      }

   }
   gStyle->SetPalette(kRainBow);


   std::vector<double> errorVector={0.005, .020, 0.005, .020};
   std::vector<double> minmaxVector={0.1,0.2,0.1,0.2};
   auto error=errorVector.begin(), minmax=minmaxVector.begin();
   for(auto d : {&deformationsy, &deformationsx, &deformationsyExp, &deformationsxExp} ) {
	   d=removeBinsWithFewerEntries(d, 1000);
	   getFrequencyHistogram(d, *error++);
	   double max=*minmax++;
	   setMinMax(d,-max,max);
   }

   output.Write();
}
