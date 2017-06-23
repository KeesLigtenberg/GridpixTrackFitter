//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 14 09:55:46 2017 by ROOT version 6.08/02
// from TTree Hits/Hits by event
// found on file: /data/lepcol/common/BonnTelescopeTestJune17/combinedConverted.root
//////////////////////////////////////////////////////////

#ifndef tmp_h
#define tmp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class tmp {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        eventNumber;
   UShort_t        triggerNumberBegin;
   UShort_t        triggerNumberEnd;
   UInt_t          timestamp;
   vector<vector<Hit> > *mimosa;

   // List of branches
   TBranch        *b_eventNumber;   //!
   TBranch        *b_triggerNumberBegin;   //!
   TBranch        *b_triggerNumberEnd;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_mimosa;   //!

   tmp(TTree *tree=0);
   virtual ~tmp();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tmp_cxx
tmp::tmp(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/lepcol/common/BonnTelescopeTestJune17/combinedConverted.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/lepcol/common/BonnTelescopeTestJune17/combinedConverted.root");
      }
      f->GetObject("Hits",tree);

   }
   Init(tree);
}

tmp::~tmp()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tmp::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tmp::LoadTree(Long64_t entry)
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

void tmp::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mimosa = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("triggerNumberBegin", &triggerNumberBegin, &b_triggerNumberBegin);
   fChain->SetBranchAddress("triggerNumberEnd", &triggerNumberEnd, &b_triggerNumberEnd);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("mimosa", &mimosa, &b_mimosa);
   Notify();
}

Bool_t tmp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tmp::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tmp::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tmp_cxx
