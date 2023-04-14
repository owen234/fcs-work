
#include "histio.c"
#include "utils.c"

   void draw_scatter1( const char* infilepat = "../analysis-tree/analysis-tree-mb-mufilter-4*.root" ) {

      gDirectory -> Delete( "h*" ) ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadLeftMargin(0.15) ;
      gStyle -> SetPadRightMargin(0.05) ;

      TChain ch("fwd") ;

      //ch.Add( infilepat ) ;
      ch.Add( "../analysis-tree/analysis-tree-mb-1a.root" ) ;
      ch.Add( "../analysis-tree/analysis-tree-mb-4a.root" ) ;
      ch.Add( "../analysis-tree/analysis-tree-mb-4b.root" ) ;
      ch.Add( "../analysis-tree/analysis-tree-mb-4c.root" ) ;

      ch.SetMarkerStyle(20) ;
      ch.SetMarkerSize(0.3) ;

      int nentries = ch.GetEntries() ;
      printf("\n\n N entries in %s :  %d\n\n", infilepat, nentries ) ;

      TH2F* h_e_vs_l = new TH2F( "h_e_vs_l", "", 500, 0., 100., 500, 0., 2. ) ;
      h_e_vs_l -> SetXTitle( "HCAL cell track path length (cm)" ) ;
      h_e_vs_l -> SetYTitle( "HCAL cell energy (GeV)" ) ;
      h_e_vs_l -> SetTitleOffset( 1.4, "y") ;

      TText* tt = new TText() ;
      float tx(0.15) ;
      float ty(0.92) ;



//    TCanvas* can1 = get_canvas( "can1", "HCAL E vs L", 50, 50, 1400, 1000 ) ;
//    can1 -> cd() ;

//    can1 -> Clear() ;
//    can1 -> Divide(4,2) ;

//    ch.SetMarkerColor(4) ;

//    can1 -> cd(1) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "muons" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && rc_ismu", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

//    can1 -> cd(2) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "muons, iso cut" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && rc_nm1_chi2_iso1_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e2 && rc_ismu", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

//    can1 -> cd(3) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "muons, energy chi2 cut" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && rc_nm1_hcal_L2>0 && rc_nm1_chi2_energy_ecal<1e1 && rc_nm1_chi2_energy_hcal<1e1 && rc_ismu", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

//    can1 -> cd(4) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "muons, iso and energy chi2 cuts" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && rc_nm1_hcal_L2>0 && rc_nm1_chi2_energy_ecal<1e1 && rc_nm1_chi2_energy_hcal<1e1 && rc_ismu && rc_nm1_chi2_iso1_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e2", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;


//    ch.SetMarkerColor(2) ;

//    can1 -> cd(5) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "Not muons" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && !rc_ismu", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

//    can1 -> cd(6) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "Not muons, iso cut" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && rc_nm1_chi2_iso1_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e2 && !rc_ismu", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

//    can1 -> cd(7) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "Not muons, energy chi2 cut" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && rc_nm1_hcal_L2>0 && rc_nm1_chi2_energy_ecal<1e1 && rc_nm1_chi2_energy_hcal<1e1 && !rc_ismu", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

//    can1 -> cd(8) ;
//    h_e_vs_l -> Draw() ;
//    tt->DrawTextNDC( tx, ty, "Not muons, iso and energy chi2 cuts" ) ;
//    ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "rc_nm1_chi2_energy_ecal>0 && rc_nm1_chi2_energy_hcal>0 && rc_nm1_nhit_hcal>1 && rc_nm1_hcal_L2>0 && rc_nm1_chi2_energy_ecal<1e1 && rc_nm1_chi2_energy_hcal<1e1 && rc_nm1_chi2_iso1_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e2 && !rc_ismu", "same" ) ;
//    can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;




   //-----------------

      TCanvas* can2 = get_canvas( "can2", "HCAL E vs L", 100, 100, 1400, 1000 ) ;
      can2 -> cd() ;

      can2 -> Clear() ;
      can2 -> Divide(4,2) ;

      can2 -> cd(1) ;
      h_e_vs_l -> Draw() ;
      ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "!rc_ismu && rc_nhit_hcal>1 && rc_energy_ecal>0.23 && rc_energy_ecal<0.34 && rc_nm1_chi2_iso1_ecal<1e0 && rc_nm1_chi2_energy_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e0 ", "same" ) ;
      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      can2 -> cd(2) ;
      h_e_vs_l -> Draw() ;
      ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "!rc_ismu && rc_nhit_hcal>1 && rc_energy_ecal>0.23 && rc_energy_ecal<0.34 && rc_nm1_chi2_iso1_ecal<1e0 && rc_nm1_chi2_energy_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e0 && rc_nm1_chi2_energy_hcal<1e1 ", "same" ) ;
      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      can2 -> cd(3) ;
      h_e_vs_l -> Draw() ;
      ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "!rc_ismu && rc_nhit_hcal>1 && rc_energy_ecal>0.23 && rc_energy_ecal<0.34 && rc_nm1_chi2_iso1_ecal<1e0 && rc_nm1_chi2_energy_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e0 && rc_energy_hcal > 1.1 && rc_energy_hcal < 1.5", "same" ) ;
      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      can2 -> cd(4) ;
      h_e_vs_l -> Draw() ;
      ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "!rc_ismu && rc_nhit_hcal>1 && rc_energy_ecal>0.23 && rc_energy_ecal<0.34 && rc_nm1_chi2_iso1_ecal<1e0 && rc_nm1_chi2_energy_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e0 && rc_energy_hcal > 1.1 && rc_energy_hcal < 1.5 && rc_nm1_chi2_energy_hcal<1e1", "same" ) ;
      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      can2 -> cd(5) ;
      h_e_vs_l -> Draw() ;
      ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "!rc_ismu && rc_nhit_hcal>1 && rc_energy_ecal>0.23 && rc_energy_ecal<0.34 && rc_nm1_chi2_iso1_ecal<1e0 && rc_nm1_chi2_energy_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e0 && rc_energy_hcal > 1.1 && rc_energy_hcal < 1.5 && rc_nm1_chi2_energy_hcal<1e1 && rc_nm1_hcal_E2>0.2", "same" ) ;
      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      can2 -> cd(6) ;
      h_e_vs_l -> Draw() ;
      ch.Draw( "rc_nm1_hcal_E1:rc_nm1_hcal_L1", "!rc_ismu && rc_nhit_hcal>1 && rc_energy_ecal>0.23 && rc_energy_ecal<0.34 && rc_nm1_chi2_iso1_ecal<1e0 && rc_nm1_chi2_energy_ecal<1e2 && rc_nm1_chi2_iso1_hcal<1e0 && rc_energy_hcal > 0.1 && rc_energy_hcal < 10.5 && rc_nm1_chi2_energy_hcal<1e1 && rc_nm1_hcal_E2>0.2", "same" ) ;
      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;



   }
















