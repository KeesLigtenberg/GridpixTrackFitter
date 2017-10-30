//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct 29 20:39:27 2017 by ROOT version 6.08/02
// from TTree fitResults/Tree with telescope and timepix fit results
// found on file: fitResults.root
//////////////////////////////////////////////////////////

#ifndef resultProcessor_h
#define resultProcessor_h

#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "HitEntry.h"
#include "linearRegressionFit.cpp"

class resultProcessor {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxtelescopeFits = 1;
   static const Int_t kMaxtimepixFits = 1;

   // Declaration of leaf types
   std::vector<FitResult3D> *telescopeFits;
   std::vector<FitResult3D> *timepixFits;

   Int_t           ntimepixHits;
   Int_t           ntelescopeHits;
   std::vector<int>     *timepixClusterSize;
   std::vector<std::vector<HitEntry> > *timepixHits;
   Double_t        dxz;
   Double_t        dyz;
   Int_t           nresiduals;
   Int_t           nfitted;

   // List of branches
   TBranch        *b_telescopeFits;
   TBranch        *b_timepixFits;

   TBranch        *b_ntimepixHits;   //!
   TBranch        *b_ntelescopeHits;   //!
   TBranch        *b_timepixClusterSize;   //!
   TBranch        *b_timepixHits;   //!
   TBranch        *b_dxz;   //!
   TBranch        *b_dyz;   //!
   TBranch        *b_nresiduals;   //!
   TBranch        *b_nfitted;   //!

   resultProcessor(TTree *tree);
   virtual ~resultProcessor();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

