//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr  9 12:34:36 2023 by ROOT version 6.24/07
// from TChain fwd/
//////////////////////////////////////////////////////////

#ifndef mufrac_h
#define mufrac_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "c++/v1/vector"
#include "c++/v1/vector"
#include "c++/v1/vector"

class mufrac {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           fttN;
   vector<float>   *fttX;
   vector<float>   *fttY;
   vector<float>   *fttZ;
   vector<float>   *fttMcpx;
   vector<float>   *fttMcpy;
   vector<float>   *fttMcpz;
   vector<int>     *fttTrackId;
   vector<int>     *fttVolumeId;
   vector<float>   *fttPt;
   vector<int>     *fttVertexId;
   Int_t           fstN;
   vector<float>   *fstX;
   vector<float>   *fstY;
   vector<float>   *fstZ;
   vector<int>     *fstTrackId;
   Int_t           mcN;
   vector<float>   *mcPt;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<int>     *mcCharge;
   vector<int>     *mcPid;
   vector<int>     *mcVertexId;
   Int_t           vmcN;
   vector<float>   *vmcX;
   vector<float>   *vmcY;
   vector<float>   *vmcZ;
   Int_t           vrcN;
   vector<float>   *vrcX;
   vector<float>   *vrcY;
   vector<float>   *vrcZ;
   Int_t           rcN;
   vector<float>   *rcPt;
   vector<float>   *rcEta;
   vector<float>   *rcPhi;
   vector<int>     *rcCharge;
   vector<int>     *rcTrackId;
   vector<int>     *rcTrackIdT;
   vector<int>     *rcNumFST;
   vector<int>     *rcNumFTT;
   vector<int>     *rcNumPV;
   vector<float>   *rcQuality;
   Int_t           thdN;
   vector<float>   *thdX;
   vector<float>   *thdY;
   vector<float>   *thaX;
   vector<float>   *thaY;
   vector<float>   *thaZ;
   Int_t           tprojN;
   vector<int>     *tprojIdD;
   vector<int>     *tprojIdT;
   vector<float>   *tprojX;
   vector<float>   *tprojY;
   vector<float>   *tprojZ;
   vector<float>   *tprojPx;
   vector<float>   *tprojPy;
   vector<float>   *tprojPz;
   Int_t           fcsN;
   vector<float>   *fcsX;
   vector<float>   *fcsY;
   vector<float>   *fcsZ;
   vector<float>   *fcsE;
   Int_t           fcs_mc_ecalN;
   vector<float>   *fcs_mc_ecalX;
   vector<float>   *fcs_mc_ecalY;
   vector<float>   *fcs_mc_ecalZ;
   vector<float>   *fcs_mc_ecalE;
   vector<int>     *fcs_mc_ecalId;
   vector<int>     *fcs_mc_ecalVid;
   Int_t           fcs_mc_hcalN;
   vector<float>   *fcs_mc_hcalX;
   vector<float>   *fcs_mc_hcalY;
   vector<float>   *fcs_mc_hcalZ;
   vector<float>   *fcs_mc_hcalE;
   vector<int>     *fcs_mc_hcalId;
   vector<int>     *fcs_mc_hcalVid;
   Int_t           fcs_rec_ecalN;
   vector<float>   *fcs_rec_ecalX;
   vector<float>   *fcs_rec_ecalY;
   vector<float>   *fcs_rec_ecalZ;
   vector<float>   *fcs_rec_ecalLX;
   vector<float>   *fcs_rec_ecalLY;
   vector<float>   *fcs_rec_ecalE;
   vector<int>     *fcs_rec_ecalId;
   Int_t           fcs_rec_hcalN;
   vector<float>   *fcs_rec_hcalX;
   vector<float>   *fcs_rec_hcalY;
   vector<float>   *fcs_rec_hcalZ;
   vector<float>   *fcs_rec_hcalLX;
   vector<float>   *fcs_rec_hcalLY;
   vector<float>   *fcs_rec_hcalE;
   vector<int>     *fcs_rec_hcalId;
   vector<float>   *Crit2_RZRatio;
   vector<int>     *Crit2_RZRatio_trackIds;
   vector<float>   *Crit2_RZRatio_x1;
   vector<float>   *Crit2_RZRatio_y1;
   vector<float>   *Crit2_RZRatio_z1;
   vector<float>   *Crit2_RZRatio_x2;
   vector<float>   *Crit2_RZRatio_y2;
   vector<float>   *Crit2_RZRatio_z2;
   vector<int>     *Crit2_RZRatio_h1;
   vector<int>     *Crit2_RZRatio_h2;
   vector<int>     *Crit2_RZRatio_h3;
   vector<float>   *Crit2_DeltaPhi;
   vector<int>     *Crit2_DeltaPhi_trackIds;
   vector<float>   *Crit2_DeltaRho;
   vector<int>     *Crit2_DeltaRho_trackIds;
   vector<float>   *Crit2_StraightTrackRatio;
   vector<int>     *Crit2_StraightTrackRatio_trackIds;
   vector<float>   *Crit3_3DAngle;
   vector<int>     *Crit3_3DAngle_trackIds;
   vector<float>   *Crit3_PT;
   vector<int>     *Crit3_PT_trackIds;
   vector<float>   *Crit3_ChangeRZRatio;
   vector<int>     *Crit3_ChangeRZRatio_trackIds;
   vector<float>   *Crit3_2DAngle;
   vector<int>     *Crit3_2DAngle_trackIds;
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

   // List of branches
   TBranch        *b_fttN;   //!
   TBranch        *b_fttX;   //!
   TBranch        *b_fttY;   //!
   TBranch        *b_fttZ;   //!
   TBranch        *b_fttMcpx;   //!
   TBranch        *b_fttMcpy;   //!
   TBranch        *b_fttMcpz;   //!
   TBranch        *b_fttTrackId;   //!
   TBranch        *b_fttVolumeId;   //!
   TBranch        *b_fttPt;   //!
   TBranch        *b_fttVertexId;   //!
   TBranch        *b_fstN;   //!
   TBranch        *b_fstX;   //!
   TBranch        *b_fstY;   //!
   TBranch        *b_fstZ;   //!
   TBranch        *b_fstTrackId;   //!
   TBranch        *b_mcN;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcCharge;   //!
   TBranch        *b_mcPid;   //!
   TBranch        *b_mcVertexId;   //!
   TBranch        *b_vmcN;   //!
   TBranch        *b_vmcX;   //!
   TBranch        *b_vmcY;   //!
   TBranch        *b_vmcZ;   //!
   TBranch        *b_vrcN;   //!
   TBranch        *b_vrcX;   //!
   TBranch        *b_vrcY;   //!
   TBranch        *b_vrcZ;   //!
   TBranch        *b_rcN;   //!
   TBranch        *b_rcPt;   //!
   TBranch        *b_rcEta;   //!
   TBranch        *b_rcPhi;   //!
   TBranch        *b_rcCharge;   //!
   TBranch        *b_rcTrackId;   //!
   TBranch        *b_rcTrackIdT;   //!
   TBranch        *b_rcNumFST;   //!
   TBranch        *b_rcNumFTT;   //!
   TBranch        *b_rcNumPV;   //!
   TBranch        *b_rcQuality;   //!
   TBranch        *b_thdN;   //!
   TBranch        *b_thdX;   //!
   TBranch        *b_thdY;   //!
   TBranch        *b_thaX;   //!
   TBranch        *b_thaY;   //!
   TBranch        *b_thaZ;   //!
   TBranch        *b_tprojN;   //!
   TBranch        *b_tprojIdD;   //!
   TBranch        *b_tprojIdT;   //!
   TBranch        *b_tprojX;   //!
   TBranch        *b_tprojY;   //!
   TBranch        *b_tprojZ;   //!
   TBranch        *b_tprojPx;   //!
   TBranch        *b_tprojPy;   //!
   TBranch        *b_tprojPz;   //!
   TBranch        *b_fcsN;   //!
   TBranch        *b_fcsX;   //!
   TBranch        *b_fcsY;   //!
   TBranch        *b_fcsZ;   //!
   TBranch        *b_fcsE;   //!
   TBranch        *b_fcs_mc_ecalN;   //!
   TBranch        *b_fcs_mc_ecalX;   //!
   TBranch        *b_fcs_mc_ecalY;   //!
   TBranch        *b_fcs_mc_ecalZ;   //!
   TBranch        *b_fcs_mc_ecalE;   //!
   TBranch        *b_fcs_mc_ecalId;   //!
   TBranch        *b_fcs_mc_ecalVid;   //!
   TBranch        *b_fcs_mc_hcalN;   //!
   TBranch        *b_fcs_mc_hcalX;   //!
   TBranch        *b_fcs_mc_hcalY;   //!
   TBranch        *b_fcs_mc_hcalZ;   //!
   TBranch        *b_fcs_mc_hcalE;   //!
   TBranch        *b_fcs_mc_hcalId;   //!
   TBranch        *b_fcs_mc_hcalVid;   //!
   TBranch        *b_fcs_rec_ecalN;   //!
   TBranch        *b_fcs_rec_ecalX;   //!
   TBranch        *b_fcs_rec_ecalY;   //!
   TBranch        *b_fcs_rec_ecalZ;   //!
   TBranch        *b_fcs_rec_ecalLX;   //!
   TBranch        *b_fcs_rec_ecalLY;   //!
   TBranch        *b_fcs_rec_ecalE;   //!
   TBranch        *b_fcs_rec_ecalId;   //!
   TBranch        *b_fcs_rec_hcalN;   //!
   TBranch        *b_fcs_rec_hcalX;   //!
   TBranch        *b_fcs_rec_hcalY;   //!
   TBranch        *b_fcs_rec_hcalZ;   //!
   TBranch        *b_fcs_rec_hcalLX;   //!
   TBranch        *b_fcs_rec_hcalLY;   //!
   TBranch        *b_fcs_rec_hcalE;   //!
   TBranch        *b_fcs_rec_hcalId;   //!
   TBranch        *b_Crit2_RZRatio;   //!
   TBranch        *b_Crit2_RZRatio_trackIds;   //!
   TBranch        *b_Crit2_RZRatio_x1;   //!
   TBranch        *b_Crit2_RZRatio_y1;   //!
   TBranch        *b_Crit2_RZRatio_z1;   //!
   TBranch        *b_Crit2_RZRatio_x2;   //!
   TBranch        *b_Crit2_RZRatio_y2;   //!
   TBranch        *b_Crit2_RZRatio_z2;   //!
   TBranch        *b_Crit2_RZRatio_h1;   //!
   TBranch        *b_Crit2_RZRatio_h2;   //!
   TBranch        *b_Crit2_RZRatio_h3;   //!
   TBranch        *b_Crit2_DeltaPhi;   //!
   TBranch        *b_Crit2_DeltaPhi_trackIds;   //!
   TBranch        *b_Crit2_DeltaRho;   //!
   TBranch        *b_Crit2_DeltaRho_trackIds;   //!
   TBranch        *b_Crit2_StraightTrackRatio;   //!
   TBranch        *b_Crit2_StraightTrackRatio_trackIds;   //!
   TBranch        *b_Crit3_3DAngle;   //!
   TBranch        *b_Crit3_3DAngle_trackIds;   //!
   TBranch        *b_Crit3_PT;   //!
   TBranch        *b_Crit3_PT_trackIds;   //!
   TBranch        *b_Crit3_ChangeRZRatio;   //!
   TBranch        *b_Crit3_ChangeRZRatio_trackIds;   //!
   TBranch        *b_Crit3_2DAngle;   //!
   TBranch        *b_Crit3_2DAngle_trackIds;   //!
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

   mufrac(TTree *tree=0);
   virtual ~mufrac();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mufrac_cxx
mufrac::mufrac(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../analysis-tree/analysis-tree-mb-mufilter-4b.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../analysis-tree/analysis-tree-mb-mufilter-4b.root");
      }
      f->GetObject("fwd",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("fwd","");
      chain->Add("../analysis-tree/analysis-tree-mb-1a.root/fwd");
      chain->Add("../analysis-tree/analysis-tree-mb-4a.root/fwd");
      chain->Add("../analysis-tree/analysis-tree-mb-4b.root/fwd");
      chain->Add("../analysis-tree/analysis-tree-mb-4c.root/fwd");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

mufrac::~mufrac()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mufrac::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mufrac::LoadTree(Long64_t entry)
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

void mufrac::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   fttX = 0;
   fttY = 0;
   fttZ = 0;
   fttMcpx = 0;
   fttMcpy = 0;
   fttMcpz = 0;
   fttTrackId = 0;
   fttVolumeId = 0;
   fttPt = 0;
   fttVertexId = 0;
   fstX = 0;
   fstY = 0;
   fstZ = 0;
   fstTrackId = 0;
   mcPt = 0;
   mcEta = 0;
   mcPhi = 0;
   mcCharge = 0;
   mcPid = 0;
   mcVertexId = 0;
   vmcX = 0;
   vmcY = 0;
   vmcZ = 0;
   vrcX = 0;
   vrcY = 0;
   vrcZ = 0;
   rcPt = 0;
   rcEta = 0;
   rcPhi = 0;
   rcCharge = 0;
   rcTrackId = 0;
   rcTrackIdT = 0;
   rcNumFST = 0;
   rcNumFTT = 0;
   rcNumPV = 0;
   rcQuality = 0;
   thdX = 0;
   thdY = 0;
   thaX = 0;
   thaY = 0;
   thaZ = 0;
   tprojIdD = 0;
   tprojIdT = 0;
   tprojX = 0;
   tprojY = 0;
   tprojZ = 0;
   tprojPx = 0;
   tprojPy = 0;
   tprojPz = 0;
   fcsX = 0;
   fcsY = 0;
   fcsZ = 0;
   fcsE = 0;
   fcs_mc_ecalX = 0;
   fcs_mc_ecalY = 0;
   fcs_mc_ecalZ = 0;
   fcs_mc_ecalE = 0;
   fcs_mc_ecalId = 0;
   fcs_mc_ecalVid = 0;
   fcs_mc_hcalX = 0;
   fcs_mc_hcalY = 0;
   fcs_mc_hcalZ = 0;
   fcs_mc_hcalE = 0;
   fcs_mc_hcalId = 0;
   fcs_mc_hcalVid = 0;
   fcs_rec_ecalX = 0;
   fcs_rec_ecalY = 0;
   fcs_rec_ecalZ = 0;
   fcs_rec_ecalLX = 0;
   fcs_rec_ecalLY = 0;
   fcs_rec_ecalE = 0;
   fcs_rec_ecalId = 0;
   fcs_rec_hcalX = 0;
   fcs_rec_hcalY = 0;
   fcs_rec_hcalZ = 0;
   fcs_rec_hcalLX = 0;
   fcs_rec_hcalLY = 0;
   fcs_rec_hcalE = 0;
   fcs_rec_hcalId = 0;
   Crit2_RZRatio = 0;
   Crit2_RZRatio_trackIds = 0;
   Crit2_RZRatio_x1 = 0;
   Crit2_RZRatio_y1 = 0;
   Crit2_RZRatio_z1 = 0;
   Crit2_RZRatio_x2 = 0;
   Crit2_RZRatio_y2 = 0;
   Crit2_RZRatio_z2 = 0;
   Crit2_RZRatio_h1 = 0;
   Crit2_RZRatio_h2 = 0;
   Crit2_RZRatio_h3 = 0;
   Crit2_DeltaPhi = 0;
   Crit2_DeltaPhi_trackIds = 0;
   Crit2_DeltaRho = 0;
   Crit2_DeltaRho_trackIds = 0;
   Crit2_StraightTrackRatio = 0;
   Crit2_StraightTrackRatio_trackIds = 0;
   Crit3_3DAngle = 0;
   Crit3_3DAngle_trackIds = 0;
   Crit3_PT = 0;
   Crit3_PT_trackIds = 0;
   Crit3_ChangeRZRatio = 0;
   Crit3_ChangeRZRatio_trackIds = 0;
   Crit3_2DAngle = 0;
   Crit3_2DAngle_trackIds = 0;
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
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fttN", &fttN, &b_fttN);
   fChain->SetBranchAddress("fttX", &fttX, &b_fttX);
   fChain->SetBranchAddress("fttY", &fttY, &b_fttY);
   fChain->SetBranchAddress("fttZ", &fttZ, &b_fttZ);
   fChain->SetBranchAddress("fttMcpx", &fttMcpx, &b_fttMcpx);
   fChain->SetBranchAddress("fttMcpy", &fttMcpy, &b_fttMcpy);
   fChain->SetBranchAddress("fttMcpz", &fttMcpz, &b_fttMcpz);
   fChain->SetBranchAddress("fttTrackId", &fttTrackId, &b_fttTrackId);
   fChain->SetBranchAddress("fttVolumeId", &fttVolumeId, &b_fttVolumeId);
   fChain->SetBranchAddress("fttPt", &fttPt, &b_fttPt);
   fChain->SetBranchAddress("fttVertexId", &fttVertexId, &b_fttVertexId);
   fChain->SetBranchAddress("fstN", &fstN, &b_fstN);
   fChain->SetBranchAddress("fstX", &fstX, &b_fstX);
   fChain->SetBranchAddress("fstY", &fstY, &b_fstY);
   fChain->SetBranchAddress("fstZ", &fstZ, &b_fstZ);
   fChain->SetBranchAddress("fstTrackId", &fstTrackId, &b_fstTrackId);
   fChain->SetBranchAddress("mcN", &mcN, &b_mcN);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcCharge", &mcCharge, &b_mcCharge);
   fChain->SetBranchAddress("mcPid", &mcPid, &b_mcPid);
   fChain->SetBranchAddress("mcVertexId", &mcVertexId, &b_mcVertexId);
   fChain->SetBranchAddress("vmcN", &vmcN, &b_vmcN);
   fChain->SetBranchAddress("vmcX", &vmcX, &b_vmcX);
   fChain->SetBranchAddress("vmcY", &vmcY, &b_vmcY);
   fChain->SetBranchAddress("vmcZ", &vmcZ, &b_vmcZ);
   fChain->SetBranchAddress("vrcN", &vrcN, &b_vrcN);
   fChain->SetBranchAddress("vrcX", &vrcX, &b_vrcX);
   fChain->SetBranchAddress("vrcY", &vrcY, &b_vrcY);
   fChain->SetBranchAddress("vrcZ", &vrcZ, &b_vrcZ);
   fChain->SetBranchAddress("rcN", &rcN, &b_rcN);
   fChain->SetBranchAddress("rcPt", &rcPt, &b_rcPt);
   fChain->SetBranchAddress("rcEta", &rcEta, &b_rcEta);
   fChain->SetBranchAddress("rcPhi", &rcPhi, &b_rcPhi);
   fChain->SetBranchAddress("rcCharge", &rcCharge, &b_rcCharge);
   fChain->SetBranchAddress("rcTrackId", &rcTrackId, &b_rcTrackId);
   fChain->SetBranchAddress("rcTrackIdT", &rcTrackIdT, &b_rcTrackIdT);
   fChain->SetBranchAddress("rcNumFST", &rcNumFST, &b_rcNumFST);
   fChain->SetBranchAddress("rcNumFTT", &rcNumFTT, &b_rcNumFTT);
   fChain->SetBranchAddress("rcNumPV", &rcNumPV, &b_rcNumPV);
   fChain->SetBranchAddress("rcQuality", &rcQuality, &b_rcQuality);
   fChain->SetBranchAddress("thdN", &thdN, &b_thdN);
   fChain->SetBranchAddress("thdX", &thdX, &b_thdX);
   fChain->SetBranchAddress("thdY", &thdY, &b_thdY);
   fChain->SetBranchAddress("thaX", &thaX, &b_thaX);
   fChain->SetBranchAddress("thaY", &thaY, &b_thaY);
   fChain->SetBranchAddress("thaZ", &thaZ, &b_thaZ);
   fChain->SetBranchAddress("tprojN", &tprojN, &b_tprojN);
   fChain->SetBranchAddress("tprojIdD", &tprojIdD, &b_tprojIdD);
   fChain->SetBranchAddress("tprojIdT", &tprojIdT, &b_tprojIdT);
   fChain->SetBranchAddress("tprojX", &tprojX, &b_tprojX);
   fChain->SetBranchAddress("tprojY", &tprojY, &b_tprojY);
   fChain->SetBranchAddress("tprojZ", &tprojZ, &b_tprojZ);
   fChain->SetBranchAddress("tprojPx", &tprojPx, &b_tprojPx);
   fChain->SetBranchAddress("tprojPy", &tprojPy, &b_tprojPy);
   fChain->SetBranchAddress("tprojPz", &tprojPz, &b_tprojPz);
   fChain->SetBranchAddress("fcsN", &fcsN, &b_fcsN);
   fChain->SetBranchAddress("fcsX", &fcsX, &b_fcsX);
   fChain->SetBranchAddress("fcsY", &fcsY, &b_fcsY);
   fChain->SetBranchAddress("fcsZ", &fcsZ, &b_fcsZ);
   fChain->SetBranchAddress("fcsE", &fcsE, &b_fcsE);
   fChain->SetBranchAddress("fcs_mc_ecalN", &fcs_mc_ecalN, &b_fcs_mc_ecalN);
   fChain->SetBranchAddress("fcs_mc_ecalX", &fcs_mc_ecalX, &b_fcs_mc_ecalX);
   fChain->SetBranchAddress("fcs_mc_ecalY", &fcs_mc_ecalY, &b_fcs_mc_ecalY);
   fChain->SetBranchAddress("fcs_mc_ecalZ", &fcs_mc_ecalZ, &b_fcs_mc_ecalZ);
   fChain->SetBranchAddress("fcs_mc_ecalE", &fcs_mc_ecalE, &b_fcs_mc_ecalE);
   fChain->SetBranchAddress("fcs_mc_ecalId", &fcs_mc_ecalId, &b_fcs_mc_ecalId);
   fChain->SetBranchAddress("fcs_mc_ecalVid", &fcs_mc_ecalVid, &b_fcs_mc_ecalVid);
   fChain->SetBranchAddress("fcs_mc_hcalN", &fcs_mc_hcalN, &b_fcs_mc_hcalN);
   fChain->SetBranchAddress("fcs_mc_hcalX", &fcs_mc_hcalX, &b_fcs_mc_hcalX);
   fChain->SetBranchAddress("fcs_mc_hcalY", &fcs_mc_hcalY, &b_fcs_mc_hcalY);
   fChain->SetBranchAddress("fcs_mc_hcalZ", &fcs_mc_hcalZ, &b_fcs_mc_hcalZ);
   fChain->SetBranchAddress("fcs_mc_hcalE", &fcs_mc_hcalE, &b_fcs_mc_hcalE);
   fChain->SetBranchAddress("fcs_mc_hcalId", &fcs_mc_hcalId, &b_fcs_mc_hcalId);
   fChain->SetBranchAddress("fcs_mc_hcalVid", &fcs_mc_hcalVid, &b_fcs_mc_hcalVid);
   fChain->SetBranchAddress("fcs_rec_ecalN", &fcs_rec_ecalN, &b_fcs_rec_ecalN);
   fChain->SetBranchAddress("fcs_rec_ecalX", &fcs_rec_ecalX, &b_fcs_rec_ecalX);
   fChain->SetBranchAddress("fcs_rec_ecalY", &fcs_rec_ecalY, &b_fcs_rec_ecalY);
   fChain->SetBranchAddress("fcs_rec_ecalZ", &fcs_rec_ecalZ, &b_fcs_rec_ecalZ);
   fChain->SetBranchAddress("fcs_rec_ecalLX", &fcs_rec_ecalLX, &b_fcs_rec_ecalLX);
   fChain->SetBranchAddress("fcs_rec_ecalLY", &fcs_rec_ecalLY, &b_fcs_rec_ecalLY);
   fChain->SetBranchAddress("fcs_rec_ecalE", &fcs_rec_ecalE, &b_fcs_rec_ecalE);
   fChain->SetBranchAddress("fcs_rec_ecalId", &fcs_rec_ecalId, &b_fcs_rec_ecalId);
   fChain->SetBranchAddress("fcs_rec_hcalN", &fcs_rec_hcalN, &b_fcs_rec_hcalN);
   fChain->SetBranchAddress("fcs_rec_hcalX", &fcs_rec_hcalX, &b_fcs_rec_hcalX);
   fChain->SetBranchAddress("fcs_rec_hcalY", &fcs_rec_hcalY, &b_fcs_rec_hcalY);
   fChain->SetBranchAddress("fcs_rec_hcalZ", &fcs_rec_hcalZ, &b_fcs_rec_hcalZ);
   fChain->SetBranchAddress("fcs_rec_hcalLX", &fcs_rec_hcalLX, &b_fcs_rec_hcalLX);
   fChain->SetBranchAddress("fcs_rec_hcalLY", &fcs_rec_hcalLY, &b_fcs_rec_hcalLY);
   fChain->SetBranchAddress("fcs_rec_hcalE", &fcs_rec_hcalE, &b_fcs_rec_hcalE);
   fChain->SetBranchAddress("fcs_rec_hcalId", &fcs_rec_hcalId, &b_fcs_rec_hcalId);
   fChain->SetBranchAddress("Crit2_RZRatio", &Crit2_RZRatio, &b_Crit2_RZRatio);
   fChain->SetBranchAddress("Crit2_RZRatio_trackIds", &Crit2_RZRatio_trackIds, &b_Crit2_RZRatio_trackIds);
   fChain->SetBranchAddress("Crit2_RZRatio_x1", &Crit2_RZRatio_x1, &b_Crit2_RZRatio_x1);
   fChain->SetBranchAddress("Crit2_RZRatio_y1", &Crit2_RZRatio_y1, &b_Crit2_RZRatio_y1);
   fChain->SetBranchAddress("Crit2_RZRatio_z1", &Crit2_RZRatio_z1, &b_Crit2_RZRatio_z1);
   fChain->SetBranchAddress("Crit2_RZRatio_x2", &Crit2_RZRatio_x2, &b_Crit2_RZRatio_x2);
   fChain->SetBranchAddress("Crit2_RZRatio_y2", &Crit2_RZRatio_y2, &b_Crit2_RZRatio_y2);
   fChain->SetBranchAddress("Crit2_RZRatio_z2", &Crit2_RZRatio_z2, &b_Crit2_RZRatio_z2);
   fChain->SetBranchAddress("Crit2_RZRatio_h1", &Crit2_RZRatio_h1, &b_Crit2_RZRatio_h1);
   fChain->SetBranchAddress("Crit2_RZRatio_h2", &Crit2_RZRatio_h2, &b_Crit2_RZRatio_h2);
   fChain->SetBranchAddress("Crit2_RZRatio_h3", &Crit2_RZRatio_h3, &b_Crit2_RZRatio_h3);
   fChain->SetBranchAddress("Crit2_DeltaPhi", &Crit2_DeltaPhi, &b_Crit2_DeltaPhi);
   fChain->SetBranchAddress("Crit2_DeltaPhi_trackIds", &Crit2_DeltaPhi_trackIds, &b_Crit2_DeltaPhi_trackIds);
   fChain->SetBranchAddress("Crit2_DeltaRho", &Crit2_DeltaRho, &b_Crit2_DeltaRho);
   fChain->SetBranchAddress("Crit2_DeltaRho_trackIds", &Crit2_DeltaRho_trackIds, &b_Crit2_DeltaRho_trackIds);
   fChain->SetBranchAddress("Crit2_StraightTrackRatio", &Crit2_StraightTrackRatio, &b_Crit2_StraightTrackRatio);
   fChain->SetBranchAddress("Crit2_StraightTrackRatio_trackIds", &Crit2_StraightTrackRatio_trackIds, &b_Crit2_StraightTrackRatio_trackIds);
   fChain->SetBranchAddress("Crit3_3DAngle", &Crit3_3DAngle, &b_Crit3_3DAngle);
   fChain->SetBranchAddress("Crit3_3DAngle_trackIds", &Crit3_3DAngle_trackIds, &b_Crit3_3DAngle_trackIds);
   fChain->SetBranchAddress("Crit3_PT", &Crit3_PT, &b_Crit3_PT);
   fChain->SetBranchAddress("Crit3_PT_trackIds", &Crit3_PT_trackIds, &b_Crit3_PT_trackIds);
   fChain->SetBranchAddress("Crit3_ChangeRZRatio", &Crit3_ChangeRZRatio, &b_Crit3_ChangeRZRatio);
   fChain->SetBranchAddress("Crit3_ChangeRZRatio_trackIds", &Crit3_ChangeRZRatio_trackIds, &b_Crit3_ChangeRZRatio_trackIds);
   fChain->SetBranchAddress("Crit3_2DAngle", &Crit3_2DAngle, &b_Crit3_2DAngle);
   fChain->SetBranchAddress("Crit3_2DAngle_trackIds", &Crit3_2DAngle_trackIds, &b_Crit3_2DAngle_trackIds);
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
   Notify();
}

Bool_t mufrac::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mufrac::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mufrac::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mufrac_cxx
