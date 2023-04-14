//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 13 17:21:11 2023 by ROOT version 6.24/07
// from TChain fwd/
//////////////////////////////////////////////////////////

#ifndef plots2_h
#define plots2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "c++/v1/vector"
#include "c++/v1/vector"
#include "c++/v1/vector"

class plots2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *rc_chi2_energy_ecal;
   vector<float>   *rc_chi2_energy_hcal;
   vector<float>   *rc_chi2_energy;
   vector<float>   *rc_chi2_iso1_ecal;
   vector<float>   *rc_chi2_iso1_hcal;
   vector<float>   *rc_chi2_iso1;
   vector<float>   *rc_chi2_iso2_ecal;
   vector<float>   *rc_chi2_iso2_hcal;
   vector<float>   *rc_chi2_iso2;
   vector<float>   *rc_dx;
   vector<float>   *rc_dy;
   vector<float>   *rc_dx_gradient;
   vector<float>   *rc_dy_gradient;
   vector<float>   *rc_energy_ecal;
   vector<float>   *rc_energy_hcal;
   vector<int>     *rc_nhit_ecal;
   vector<int>     *rc_nhit_hcal;
   vector<float>   *rc_nm1_chi2_energy_ecal;
   vector<float>   *rc_nm1_chi2_energy_hcal;
   vector<float>   *rc_nm1_chi2_energy;
   vector<float>   *rc_nm1_chi2_iso1_ecal;
   vector<float>   *rc_nm1_chi2_iso1_hcal;
   vector<float>   *rc_nm1_chi2_iso1;
   vector<float>   *rc_nm1_chi2_iso2_ecal;
   vector<float>   *rc_nm1_chi2_iso2_hcal;
   vector<float>   *rc_nm1_chi2_iso2;
   vector<float>   *rc_nm1_dx;
   vector<float>   *rc_nm1_dy;
   vector<float>   *rc_nm1_dx_gradient;
   vector<float>   *rc_nm1_dy_gradient;
   vector<float>   *rc_nm1_energy_ecal;
   vector<float>   *rc_nm1_energy_hcal;
   vector<int>     *rc_nm1_nhit_ecal;
   vector<int>     *rc_nm1_nhit_hcal;
   vector<float>   *rc_nm1_hcal_L1;
   vector<float>   *rc_nm1_hcal_E1;
   vector<float>   *rc_nm1_hcal_L2;
   vector<float>   *rc_nm1_hcal_E2;
   vector<int>     *rc_nm1_hcal_cell_id;
   vector<int>     *rc_nm1_hcal_cell_row;
   vector<int>     *rc_nm1_hcal_cell_col;
   vector<int>     *rc_nm1_hcal_cell_ns;
   vector<float>   *rc_Pt;
   vector<float>   *rc_Eta;
   vector<float>   *rc_Phi;
   vector<int>     *rc_mcPid;
   vector<int>     *rc_rci;
   vector<bool>    *rc_ismu;
   vector<float>   *rc_P;
   vector<float>   *rc_Charge;
   vector<float>   *rc_mcPt;
   vector<float>   *rc_mcEta;
   vector<float>   *rc_mcPhi;
   vector<float>   *rc_mcdR;
   vector<float>   *rc_mcP;
   vector<float>   *rc_mcCharge;

   // List of branches
   TBranch        *b_rc_chi2_energy_ecal;   //!
   TBranch        *b_rc_chi2_energy_hcal;   //!
   TBranch        *b_rc_chi2_energy;   //!
   TBranch        *b_rc_chi2_iso1_ecal;   //!
   TBranch        *b_rc_chi2_iso1_hcal;   //!
   TBranch        *b_rc_chi2_iso1;   //!
   TBranch        *b_rc_chi2_iso2_ecal;   //!
   TBranch        *b_rc_chi2_iso2_hcal;   //!
   TBranch        *b_rc_chi2_iso2;   //!
   TBranch        *b_rc_dx;   //!
   TBranch        *b_rc_dy;   //!
   TBranch        *b_rc_dx_gradient;   //!
   TBranch        *b_rc_dy_gradient;   //!
   TBranch        *b_rc_energy_ecal;   //!
   TBranch        *b_rc_energy_hcal;   //!
   TBranch        *b_rc_nhit_ecal;   //!
   TBranch        *b_rc_nhit_hcal;   //!
   TBranch        *b_rc_nm1_chi2_energy_ecal;   //!
   TBranch        *b_rc_nm1_chi2_energy_hcal;   //!
   TBranch        *b_rc_nm1_chi2_energy;   //!
   TBranch        *b_rc_nm1_chi2_iso1_ecal;   //!
   TBranch        *b_rc_nm1_chi2_iso1_hcal;   //!
   TBranch        *b_rc_nm1_chi2_iso1;   //!
   TBranch        *b_rc_nm1_chi2_iso2_ecal;   //!
   TBranch        *b_rc_nm1_chi2_iso2_hcal;   //!
   TBranch        *b_rc_nm1_chi2_iso2;   //!
   TBranch        *b_rc_nm1_dx;   //!
   TBranch        *b_rc_nm1_dy;   //!
   TBranch        *b_rc_nm1_dx_gradient;   //!
   TBranch        *b_rc_nm1_dy_gradient;   //!
   TBranch        *b_rc_nm1_energy_ecal;   //!
   TBranch        *b_rc_nm1_energy_hcal;   //!
   TBranch        *b_rc_nm1_nhit_ecal;   //!
   TBranch        *b_rc_nm1_nhit_hcal;   //!
   TBranch        *b_rc_nm1_hcal_L1;   //!
   TBranch        *b_rc_nm1_hcal_E1;   //!
   TBranch        *b_rc_nm1_hcal_L2;   //!
   TBranch        *b_rc_nm1_hcal_E2;   //!
   TBranch        *b_rc_nm1_hcal_cell_id;   //!
   TBranch        *b_rc_nm1_hcal_cell_row;   //!
   TBranch        *b_rc_nm1_hcal_cell_col;   //!
   TBranch        *b_rc_nm1_hcal_cell_ns;   //!
   TBranch        *b_rc_Pt;   //!
   TBranch        *b_rc_Eta;   //!
   TBranch        *b_rc_Phi;   //!
   TBranch        *b_rc_mcPid;   //!
   TBranch        *b_rc_rci;   //!
   TBranch        *b_rc_ismu;   //!
   TBranch        *b_rc_P;   //!
   TBranch        *b_rc_Charge;   //!
   TBranch        *b_rc_mcPt;   //!
   TBranch        *b_rc_mcEta;   //!
   TBranch        *b_rc_mcPhi;   //!
   TBranch        *b_rc_mcdR;   //!
   TBranch        *b_rc_mcP;   //!
   TBranch        *b_rc_mcCharge;   //!

   plots2(TTree *tree=0);
   virtual ~plots2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plots2_cxx
plots2::plots2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("fwd",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("fwd","");
      //chain->Add("mb/analysis-tree-job-*.root/fwd");
      chain->Add("mu-filter-b/analysis-tree-job-*.root/fwd");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

plots2::~plots2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t plots2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plots2::LoadTree(Long64_t entry)
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

void plots2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   rc_chi2_energy_ecal = 0;
   rc_chi2_energy_hcal = 0;
   rc_chi2_energy = 0;
   rc_chi2_iso1_ecal = 0;
   rc_chi2_iso1_hcal = 0;
   rc_chi2_iso1 = 0;
   rc_chi2_iso2_ecal = 0;
   rc_chi2_iso2_hcal = 0;
   rc_chi2_iso2 = 0;
   rc_dx = 0;
   rc_dy = 0;
   rc_dx_gradient = 0;
   rc_dy_gradient = 0;
   rc_energy_ecal = 0;
   rc_energy_hcal = 0;
   rc_nhit_ecal = 0;
   rc_nhit_hcal = 0;
   rc_nm1_chi2_energy_ecal = 0;
   rc_nm1_chi2_energy_hcal = 0;
   rc_nm1_chi2_energy = 0;
   rc_nm1_chi2_iso1_ecal = 0;
   rc_nm1_chi2_iso1_hcal = 0;
   rc_nm1_chi2_iso1 = 0;
   rc_nm1_chi2_iso2_ecal = 0;
   rc_nm1_chi2_iso2_hcal = 0;
   rc_nm1_chi2_iso2 = 0;
   rc_nm1_dx = 0;
   rc_nm1_dy = 0;
   rc_nm1_dx_gradient = 0;
   rc_nm1_dy_gradient = 0;
   rc_nm1_energy_ecal = 0;
   rc_nm1_energy_hcal = 0;
   rc_nm1_nhit_ecal = 0;
   rc_nm1_nhit_hcal = 0;
   rc_nm1_hcal_L1 = 0;
   rc_nm1_hcal_E1 = 0;
   rc_nm1_hcal_L2 = 0;
   rc_nm1_hcal_E2 = 0;
   rc_nm1_hcal_cell_id = 0;
   rc_nm1_hcal_cell_row = 0;
   rc_nm1_hcal_cell_col = 0;
   rc_nm1_hcal_cell_ns = 0;
   rc_Pt = 0;
   rc_Eta = 0;
   rc_Phi = 0;
   rc_mcPid = 0;
   rc_rci = 0;
   rc_ismu = 0;
   rc_P = 0;
   rc_Charge = 0;
   rc_mcPt = 0;
   rc_mcEta = 0;
   rc_mcPhi = 0;
   rc_mcdR = 0;
   rc_mcP = 0;
   rc_mcCharge = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rc_chi2_energy_ecal", &rc_chi2_energy_ecal, &b_rc_chi2_energy_ecal);
   fChain->SetBranchAddress("rc_chi2_energy_hcal", &rc_chi2_energy_hcal, &b_rc_chi2_energy_hcal);
   fChain->SetBranchAddress("rc_chi2_energy", &rc_chi2_energy, &b_rc_chi2_energy);
   fChain->SetBranchAddress("rc_chi2_iso1_ecal", &rc_chi2_iso1_ecal, &b_rc_chi2_iso1_ecal);
   fChain->SetBranchAddress("rc_chi2_iso1_hcal", &rc_chi2_iso1_hcal, &b_rc_chi2_iso1_hcal);
   fChain->SetBranchAddress("rc_chi2_iso1", &rc_chi2_iso1, &b_rc_chi2_iso1);
   fChain->SetBranchAddress("rc_chi2_iso2_ecal", &rc_chi2_iso2_ecal, &b_rc_chi2_iso2_ecal);
   fChain->SetBranchAddress("rc_chi2_iso2_hcal", &rc_chi2_iso2_hcal, &b_rc_chi2_iso2_hcal);
   fChain->SetBranchAddress("rc_chi2_iso2", &rc_chi2_iso2, &b_rc_chi2_iso2);
   fChain->SetBranchAddress("rc_dx", &rc_dx, &b_rc_dx);
   fChain->SetBranchAddress("rc_dy", &rc_dy, &b_rc_dy);
   fChain->SetBranchAddress("rc_dx_gradient", &rc_dx_gradient, &b_rc_dx_gradient);
   fChain->SetBranchAddress("rc_dy_gradient", &rc_dy_gradient, &b_rc_dy_gradient);
   fChain->SetBranchAddress("rc_energy_ecal", &rc_energy_ecal, &b_rc_energy_ecal);
   fChain->SetBranchAddress("rc_energy_hcal", &rc_energy_hcal, &b_rc_energy_hcal);
   fChain->SetBranchAddress("rc_nhit_ecal", &rc_nhit_ecal, &b_rc_nhit_ecal);
   fChain->SetBranchAddress("rc_nhit_hcal", &rc_nhit_hcal, &b_rc_nhit_hcal);
   fChain->SetBranchAddress("rc_nm1_chi2_energy_ecal", &rc_nm1_chi2_energy_ecal, &b_rc_nm1_chi2_energy_ecal);
   fChain->SetBranchAddress("rc_nm1_chi2_energy_hcal", &rc_nm1_chi2_energy_hcal, &b_rc_nm1_chi2_energy_hcal);
   fChain->SetBranchAddress("rc_nm1_chi2_energy", &rc_nm1_chi2_energy, &b_rc_nm1_chi2_energy);
   fChain->SetBranchAddress("rc_nm1_chi2_iso1_ecal", &rc_nm1_chi2_iso1_ecal, &b_rc_nm1_chi2_iso1_ecal);
   fChain->SetBranchAddress("rc_nm1_chi2_iso1_hcal", &rc_nm1_chi2_iso1_hcal, &b_rc_nm1_chi2_iso1_hcal);
   fChain->SetBranchAddress("rc_nm1_chi2_iso1", &rc_nm1_chi2_iso1, &b_rc_nm1_chi2_iso1);
   fChain->SetBranchAddress("rc_nm1_chi2_iso2_ecal", &rc_nm1_chi2_iso2_ecal, &b_rc_nm1_chi2_iso2_ecal);
   fChain->SetBranchAddress("rc_nm1_chi2_iso2_hcal", &rc_nm1_chi2_iso2_hcal, &b_rc_nm1_chi2_iso2_hcal);
   fChain->SetBranchAddress("rc_nm1_chi2_iso2", &rc_nm1_chi2_iso2, &b_rc_nm1_chi2_iso2);
   fChain->SetBranchAddress("rc_nm1_dx", &rc_nm1_dx, &b_rc_nm1_dx);
   fChain->SetBranchAddress("rc_nm1_dy", &rc_nm1_dy, &b_rc_nm1_dy);
   fChain->SetBranchAddress("rc_nm1_dx_gradient", &rc_nm1_dx_gradient, &b_rc_nm1_dx_gradient);
   fChain->SetBranchAddress("rc_nm1_dy_gradient", &rc_nm1_dy_gradient, &b_rc_nm1_dy_gradient);
   fChain->SetBranchAddress("rc_nm1_energy_ecal", &rc_nm1_energy_ecal, &b_rc_nm1_energy_ecal);
   fChain->SetBranchAddress("rc_nm1_energy_hcal", &rc_nm1_energy_hcal, &b_rc_nm1_energy_hcal);
   fChain->SetBranchAddress("rc_nm1_nhit_ecal", &rc_nm1_nhit_ecal, &b_rc_nm1_nhit_ecal);
   fChain->SetBranchAddress("rc_nm1_nhit_hcal", &rc_nm1_nhit_hcal, &b_rc_nm1_nhit_hcal);
   fChain->SetBranchAddress("rc_nm1_hcal_L1", &rc_nm1_hcal_L1, &b_rc_nm1_hcal_L1);
   fChain->SetBranchAddress("rc_nm1_hcal_E1", &rc_nm1_hcal_E1, &b_rc_nm1_hcal_E1);
   fChain->SetBranchAddress("rc_nm1_hcal_L2", &rc_nm1_hcal_L2, &b_rc_nm1_hcal_L2);
   fChain->SetBranchAddress("rc_nm1_hcal_E2", &rc_nm1_hcal_E2, &b_rc_nm1_hcal_E2);
   fChain->SetBranchAddress("rc_nm1_hcal_cell_id", &rc_nm1_hcal_cell_id, &b_rc_nm1_hcal_cell_id);
   fChain->SetBranchAddress("rc_nm1_hcal_cell_row", &rc_nm1_hcal_cell_row, &b_rc_nm1_hcal_cell_row);
   fChain->SetBranchAddress("rc_nm1_hcal_cell_col", &rc_nm1_hcal_cell_col, &b_rc_nm1_hcal_cell_col);
   fChain->SetBranchAddress("rc_nm1_hcal_cell_ns", &rc_nm1_hcal_cell_ns, &b_rc_nm1_hcal_cell_ns);
   fChain->SetBranchAddress("rc_Pt", &rc_Pt, &b_rc_Pt);
   fChain->SetBranchAddress("rc_Eta", &rc_Eta, &b_rc_Eta);
   fChain->SetBranchAddress("rc_Phi", &rc_Phi, &b_rc_Phi);
   fChain->SetBranchAddress("rc_mcPid", &rc_mcPid, &b_rc_mcPid);
   fChain->SetBranchAddress("rc_rci", &rc_rci, &b_rc_rci);
   fChain->SetBranchAddress("rc_ismu", &rc_ismu, &b_rc_ismu);
   fChain->SetBranchAddress("rc_P", &rc_P, &b_rc_P);
   fChain->SetBranchAddress("rc_Charge", &rc_Charge, &b_rc_Charge);
   fChain->SetBranchAddress("rc_mcPt", &rc_mcPt, &b_rc_mcPt);
   fChain->SetBranchAddress("rc_mcEta", &rc_mcEta, &b_rc_mcEta);
   fChain->SetBranchAddress("rc_mcPhi", &rc_mcPhi, &b_rc_mcPhi);
   fChain->SetBranchAddress("rc_mcdR", &rc_mcdR, &b_rc_mcdR);
   fChain->SetBranchAddress("rc_mcP", &rc_mcP, &b_rc_mcP);
   fChain->SetBranchAddress("rc_mcCharge", &rc_mcCharge, &b_rc_mcCharge);
   Notify();
}

Bool_t plots2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plots2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plots2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plots2_cxx
