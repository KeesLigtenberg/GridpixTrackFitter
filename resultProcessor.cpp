#include "resultProcessor.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile2D.h"

resultProcessor::resultProcessor(TTree *tree) : fChain(0)
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
   return fChain->GetEntry(entry);
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
   if(timepixClusterSize->front()/double(ntimepixHits)<0.60) return -1;
   if(timepixClusterSize->front()<30 or timepixClusterSize->front()>300) return -1;

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
			if( h->GetBinEntries(bin) < minEntries or !isInsideArea( h->GetXaxis()->GetBinCenter(x), h->GetYaxis()->GetBinCenter(y) )) {
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
void resultProcessor::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   TFile output("histograms.root", "RECREATE");

   TProfile ToTByCol("ToTByCol", "ToT by column", 256,0,256, 0, 4);
   const int nbins=64;
   TProfile2D deformationsyExp("deformationsyExp", "profile of y residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsxExp("deformationsxExp", "profile of x residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsy("deformationsy", "profile of y residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TProfile2D deformationsx("deformationsx", "profile of x residuals", nbins, 0, 256, nbins, 0, 256, -1, 1);
   TH2D diffusionx("diffusionx", "x residuals as a function of drift distance", 50,4,20,40,-2,2);
   TH2D diffusiony("diffusiony", "y residuals as a function of drift distance", 50,4,20,40,-2,2);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       if (Cut(ientry) < 0) continue;

      for(auto& h : timepixHits->front() ) {
    	  if(h.ToT*0.025<0.15 or h.flag<0) continue;
    	  double hryp=h.ry/cos(timepixFits->front().YZ.slope);
    	  deformationsxExp.Fill( h.col-h.rz/.055, h.row+h.ry/.055, h.rx );
    	  deformationsyExp.Fill( h.col-h.rz/.055, h.row+h.ry/.055, h.ry/cos(timepixFits->front().YZ.slope) ); //todo: check where XZ slope enters
    	  deformationsx.Fill( h.col, h.row+hryp/0.055, h.rx );
    	  deformationsy.Fill( h.col, h.row+hryp/0.055, h.ry/cos(timepixFits->front().YZ.slope) ); //todo: check where XZ slope enters

    	  diffusiony.Fill(h.x-h.rx, h.ry);
    	  diffusionx.Fill(h.x-h.rx, h.rx);

    	  ToTByCol.Fill(h.col, h.ToT*0.025);
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
