//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 11 09:12:19 2017 by ROOT version 6.08/06
// from TTree hits/wat
// found on file: W0015_H04-170411-114133-55.root
//////////////////////////////////////////////////////////

#ifndef TimePixTreeReader_h
#define TimePixTreeReader_h

#include <iostream>
#include <string>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "../trackFitter/getHistFromTree.h"
#include "..trackFitter/getMeanFromGausFit.h"

#include "../trackFitter/Hit.h"


struct HitTree {

	HitTree() :
		tree( "Hits", "Hits by event" ),
		timepixHits()
	{
		tree.Branch("triggerNumber", &trigger_number);
		tree.Branch("trigger",&trigger);
		tree.Branch("timepix", &timepixHits);
	}

	void Fill() {
		tree.Fill();
		timepixHits.clear();
	}

	TTree tree;

	Long64_t trigger_number=-1;
//	UShort_t triggerNumber;
	ULong64_t trigger=0;
	std::vector<TimePixHit> timepixHits;

};

// Header file for the classes stored in the TTree if any.

class TimePixTreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UShort_t        no_hits;
   static constexpr int maxNumberOfHitsPerEvent=32000;
   UChar_t         Col[maxNumberOfHitsPerEvent];   //[no_hits]
   UChar_t         Row[maxNumberOfHitsPerEvent];   //[no_hits]
   ULong64_t       ToA[maxNumberOfHitsPerEvent];   //[no_hits]
   UShort_t        ToT[maxNumberOfHitsPerEvent];   //[no_hits]
   ULong64_t       trigger_time;
   ULong64_t       trigger_number;

   Long64_t        Drift_time[maxNumberOfHitsPerEvent];

   // List of branches
   TBranch        *b_no_hits;   //!
   TBranch        *b_Col;   //!
   TBranch        *b_Row;   //!
   TBranch        *b_ToA;   //!
   TBranch        *b_ToT;   //!
   TBranch        *b_trigger_time;   //!
   TBranch        *b_trigger_number;   //!

   TimePixTreeReader(TTree *tree);
   virtual ~TimePixTreeReader();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outputfile);
   virtual void     Show(Long64_t entry = -1);

private:
   double findTriggerOffset(double estimate, double range);
};


using namespace std;

TimePixTreeReader::TimePixTreeReader(TTree *tree) : fChain(0)
{
   Init(tree);
}

TimePixTreeReader::~TimePixTreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TimePixTreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TimePixTreeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void TimePixTreeReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("no_hits", &no_hits, &b_no_hits);
   fChain->SetBranchAddress("Col", Col, &b_Col);
   fChain->SetBranchAddress("Row", Row, &b_Row);
   fChain->SetBranchAddress("ToA", ToA, &b_ToA);
   fChain->SetBranchAddress("ToT", ToT, &b_ToT);
   fChain->SetBranchAddress("trigger_time", &trigger_time, &b_trigger_time);
   fChain->SetBranchAddress("trigger_number", &trigger_number, &b_trigger_number);

}

void TimePixTreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

double TimePixTreeReader::findTriggerOffset( double estimate, double range ) {
	using std::to_string;
	//plot drift time with cut on range
	TH1* hist=getHistFromTree(*fChain, "(ToA-trigger_time)*25./4096-"+to_string(estimate), "fabs( (ToA-trigger_time)*25./4096-"+to_string(estimate)+")<"+to_string(range));
	double correction=getMeanFromGausFit(*hist);
//	delete hist;
	cin.get();
	return estimate+correction;
}

void TimePixTreeReader::Loop(std::string outputfile)
{
   if (fChain == 0) return;
   cout<<"Looping over tree"<<endl;

   Long64_t nEntry = fChain->GetEntriesFast();

//   const double offSetEstimate=-207300, //-124500,
//		   offsetSearchRange=3e3; //ns!
//   double triggerOffset=findTriggerOffset(offSetEstimate, offsetSearchRange);
//   long long driftOffset=triggerOffset*4096/25.; //from ns->int64
   long long driftOffset=0;


   TFile* outfile=new TFile(outputfile.c_str(), "RECREATE");
   if(!outfile) return;
   HitTree outputTree;

   Long64_t nbytes = 0, nb = 0;
   bool firstEvent=true;
   ULong64_t firstToA=0;
   for( Long64_t iEntry=0; iEntry<nEntry; iEntry++) {
	   fChain->GetEntry(iEntry);
	   if(iEntry and !(iEntry%100000) ) cout<<"entry "<<iEntry<<endl;

	   outputTree.trigger=trigger_time;
	   outputTree.trigger_number=trigger_number;

	   for (int jhit=0; jhit<no_hits;jhit++) {

		  //add hits in event
		  Drift_time[jhit]=ToA[jhit]-trigger_time;
		  Drift_time[jhit]-=driftOffset;
		  const double maxDriftTime=1e3; //ns
		  if(fabs(Drift_time[jhit]*25/4096.)>maxDriftTime) continue;

//		  if(firstEvent) {
//			  firstToA=ToA[jhit];
//			  std::cout<<"first ToA "<<firstToA<<std::endl;
//			  firstEvent=false;
//		  }
//		  ToA[jhit]-=firstToA;

		  TimePixHit hit{ Row[jhit],  Col[jhit], ToT[jhit], int( Drift_time[jhit] ), ToA[jhit] };
		  outputTree.timepixHits.push_back( hit );
	   }

		outputTree.Fill();
	}

	outfile->Write();
	outfile->Close();

	cout<<"found "<<trigger_number<<" triggers in "<<nEntry<<" entries"<<endl;
}

#endif

