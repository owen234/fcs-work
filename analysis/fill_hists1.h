//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 31 12:51:21 2023 by ROOT version 6.24/07
// from TChain fwd/
//////////////////////////////////////////////////////////

#ifndef fill_hists1_h
#define fill_hists1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "c++/v1/vector"
#include "c++/v1/vector"

class fill_hists1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           rcN;
   vector<float>   *rcPt;
   vector<float>   *rcEta;
   vector<float>   *rcPhi;
   vector<int>     *rcNumFST;
   vector<int>     *rcNumFTT;
   vector<float>   *rcProjEcalx;
   vector<float>   *rcProjEcaly;
   vector<float>   *rcProjEcalz;
   vector<float>   *rcProjHcalx;
   vector<float>   *rcProjHcaly;
   vector<float>   *rcProjHcalz;
   vector<float>   *rcProjEcalPx;
   vector<float>   *rcProjEcalPy;
   vector<float>   *rcProjEcalPz;
   Int_t           fcs_rec_ecalN;
   vector<float>   *fcs_rec_ecalX;
   vector<float>   *fcs_rec_ecalY;
   vector<float>   *fcs_rec_ecalZ;
   vector<float>   *fcs_rec_ecalLX;
   vector<float>   *fcs_rec_ecalLY;
   vector<float>   *fcs_rec_ecalE;
   vector<int>     *fcs_rec_ecalId;
   vector<int>     *fcs_rec_ecalDet;
   Int_t           fcs_rec_hcalN;
   vector<float>   *fcs_rec_hcalX;
   vector<float>   *fcs_rec_hcalY;
   vector<float>   *fcs_rec_hcalZ;
   vector<float>   *fcs_rec_hcalLX;
   vector<float>   *fcs_rec_hcalLY;
   vector<float>   *fcs_rec_hcalE;
   vector<int>     *fcs_rec_hcalId;
   vector<int>     *fcs_rec_hcalDet;

   // List of branches
   TBranch        *b_rcN;   //!
   TBranch        *b_rcPt;   //!
   TBranch        *b_rcEta;   //!
   TBranch        *b_rcPhi;   //!
   TBranch        *b_rcNumFST;   //!
   TBranch        *b_rcNumFTT;   //!
   TBranch        *b_rcProjEcalx;   //!
   TBranch        *b_rcProjEcaly;   //!
   TBranch        *b_rcProjEcalz;   //!
   TBranch        *b_rcProjHcalx;   //!
   TBranch        *b_rcProjHcaly;   //!
   TBranch        *b_rcProjHcalz;   //!
   TBranch        *b_rcProjEcalPx;   //!
   TBranch        *b_rcProjEcalPy;   //!
   TBranch        *b_rcProjEcalPz;   //!
   TBranch        *b_fcs_rec_ecalN;   //!
   TBranch        *b_fcs_rec_ecalX;   //!
   TBranch        *b_fcs_rec_ecalY;   //!
   TBranch        *b_fcs_rec_ecalZ;   //!
   TBranch        *b_fcs_rec_ecalLX;   //!
   TBranch        *b_fcs_rec_ecalLY;   //!
   TBranch        *b_fcs_rec_ecalE;   //!
   TBranch        *b_fcs_rec_ecalId;   //!
   TBranch        *b_fcs_rec_ecalDet;   //!
   TBranch        *b_fcs_rec_hcalN;   //!
   TBranch        *b_fcs_rec_hcalX;   //!
   TBranch        *b_fcs_rec_hcalY;   //!
   TBranch        *b_fcs_rec_hcalZ;   //!
   TBranch        *b_fcs_rec_hcalLX;   //!
   TBranch        *b_fcs_rec_hcalLY;   //!
   TBranch        *b_fcs_rec_hcalE;   //!
   TBranch        *b_fcs_rec_hcalId;   //!
   TBranch        *b_fcs_rec_hcalDet;   //!

   fill_hists1(TTree *tree=0);
   virtual ~fill_hists1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop( bool verb=false, int first_event = 0 );
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fill_hists1_cxx
fill_hists1::fill_hists1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

      TChain * chain = new TChain("mipanalysis","");
      chain->Add("mip-analysis-ttree.root");
      tree = chain;

   }
   Init(tree);
}

fill_hists1::~fill_hists1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fill_hists1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fill_hists1::LoadTree(Long64_t entry)
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

void fill_hists1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   rcPt = 0;
   rcEta = 0;
   rcPhi = 0;
   rcNumFST = 0;
   rcNumFTT = 0;
   rcProjEcalx = 0;
   rcProjEcaly = 0;
   rcProjEcalz = 0;
   rcProjHcalx = 0;
   rcProjHcaly = 0;
   rcProjHcalz = 0;
   rcProjEcalPx = 0;
   rcProjEcalPy = 0;
   rcProjEcalPz = 0;
   fcs_rec_ecalX = 0;
   fcs_rec_ecalY = 0;
   fcs_rec_ecalZ = 0;
   fcs_rec_ecalLX = 0;
   fcs_rec_ecalLY = 0;
   fcs_rec_ecalE = 0;
   fcs_rec_ecalId = 0;
   fcs_rec_ecalDet = 0;
   fcs_rec_hcalX = 0;
   fcs_rec_hcalY = 0;
   fcs_rec_hcalZ = 0;
   fcs_rec_hcalLX = 0;
   fcs_rec_hcalLY = 0;
   fcs_rec_hcalE = 0;
   fcs_rec_hcalId = 0;
   fcs_rec_hcalDet = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rcN", &rcN, &b_rcN);
   fChain->SetBranchAddress("rcPt", &rcPt, &b_rcPt);
   fChain->SetBranchAddress("rcEta", &rcEta, &b_rcEta);
   fChain->SetBranchAddress("rcPhi", &rcPhi, &b_rcPhi);
   fChain->SetBranchAddress("rcNumFST", &rcNumFST, &b_rcNumFST);
   fChain->SetBranchAddress("rcNumFTT", &rcNumFTT, &b_rcNumFTT);
   fChain->SetBranchAddress("rcProjEcalx", &rcProjEcalx, &b_rcProjEcalx);
   fChain->SetBranchAddress("rcProjEcaly", &rcProjEcaly, &b_rcProjEcaly);
   fChain->SetBranchAddress("rcProjEcalz", &rcProjEcalz, &b_rcProjEcalz);
   fChain->SetBranchAddress("rcProjHcalx", &rcProjHcalx, &b_rcProjHcalx);
   fChain->SetBranchAddress("rcProjHcaly", &rcProjHcaly, &b_rcProjHcaly);
   fChain->SetBranchAddress("rcProjHcalz", &rcProjHcalz, &b_rcProjHcalz);
   fChain->SetBranchAddress("rcProjEcalPx", &rcProjEcalPx, &b_rcProjEcalPx);
   fChain->SetBranchAddress("rcProjEcalPy", &rcProjEcalPy, &b_rcProjEcalPy);
   fChain->SetBranchAddress("rcProjEcalPz", &rcProjEcalPz, &b_rcProjEcalPz);
   fChain->SetBranchAddress("fcs_rec_ecalN", &fcs_rec_ecalN, &b_fcs_rec_ecalN);
   fChain->SetBranchAddress("fcs_rec_ecalX", &fcs_rec_ecalX, &b_fcs_rec_ecalX);
   fChain->SetBranchAddress("fcs_rec_ecalY", &fcs_rec_ecalY, &b_fcs_rec_ecalY);
   fChain->SetBranchAddress("fcs_rec_ecalZ", &fcs_rec_ecalZ, &b_fcs_rec_ecalZ);
   fChain->SetBranchAddress("fcs_rec_ecalLX", &fcs_rec_ecalLX, &b_fcs_rec_ecalLX);
   fChain->SetBranchAddress("fcs_rec_ecalLY", &fcs_rec_ecalLY, &b_fcs_rec_ecalLY);
   fChain->SetBranchAddress("fcs_rec_ecalE", &fcs_rec_ecalE, &b_fcs_rec_ecalE);
   fChain->SetBranchAddress("fcs_rec_ecalId", &fcs_rec_ecalId, &b_fcs_rec_ecalId);
   fChain->SetBranchAddress("fcs_rec_ecalDet", &fcs_rec_ecalDet, &b_fcs_rec_ecalDet);
   fChain->SetBranchAddress("fcs_rec_hcalN", &fcs_rec_hcalN, &b_fcs_rec_hcalN);
   fChain->SetBranchAddress("fcs_rec_hcalX", &fcs_rec_hcalX, &b_fcs_rec_hcalX);
   fChain->SetBranchAddress("fcs_rec_hcalY", &fcs_rec_hcalY, &b_fcs_rec_hcalY);
   fChain->SetBranchAddress("fcs_rec_hcalZ", &fcs_rec_hcalZ, &b_fcs_rec_hcalZ);
   fChain->SetBranchAddress("fcs_rec_hcalLX", &fcs_rec_hcalLX, &b_fcs_rec_hcalLX);
   fChain->SetBranchAddress("fcs_rec_hcalLY", &fcs_rec_hcalLY, &b_fcs_rec_hcalLY);
   fChain->SetBranchAddress("fcs_rec_hcalE", &fcs_rec_hcalE, &b_fcs_rec_hcalE);
   fChain->SetBranchAddress("fcs_rec_hcalId", &fcs_rec_hcalId, &b_fcs_rec_hcalId);
   fChain->SetBranchAddress("fcs_rec_hcalDet", &fcs_rec_hcalDet, &b_fcs_rec_hcalDet);
   Notify();
}

Bool_t fill_hists1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fill_hists1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fill_hists1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fill_hists1_cxx
