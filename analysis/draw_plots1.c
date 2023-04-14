

#include "histio.c"
#include "utils.c"


   void draw_plots1( const char* input_file = "plots.root" ) {

      gDirectory -> Delete( "h*" ) ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetStatW(0.60) ;
      gStyle -> SetStatH(0.25) ;
      gStyle -> SetStatY(0.92) ;

      loadHist( input_file ) ;

      int canwid = 900 ;
      int can_dxy = 50 ;
      int can_w_per_plot = 400 ;
      int can_h_per_plot = 500 ;

      int canx0(0) ;
      int cany0(0) ;

      int ny ;
      int nx ;

      int ci ;

      TCanvas* can ;

    //------------------------------------

      nx = 3 ;
      ny = 2 ;

      gStyle -> SetOptStat("emr") ;

      canx0 += can_dxy ;
      cany0 += can_dxy ;

      TCanvas* can1 = get_canvas( "can1", "chi2 plots, energy", canx0, cany0, nx*can_w_per_plot, ny*can_h_per_plot ) ;

      can = can1 ;

      ci = 1 ;

      can -> Clear() ;
      can -> cd() ;
      can -> Divide( 3, ny ) ;


      TH1F* h_chi2_energy_ecal_mu = get_hist( "h_chi2_energy_ecal_mu" ) ;
      TH1F* h_chi2_energy_hcal_mu = get_hist( "h_chi2_energy_hcal_mu" ) ;
      TH1F* h_chi2_energy_mu = get_hist( "h_chi2_energy_mu" ) ;

      TH1F* h_chi2_energy_ecal_notmu = get_hist( "h_chi2_energy_ecal_notmu" ) ;
      TH1F* h_chi2_energy_hcal_notmu = get_hist( "h_chi2_energy_hcal_notmu" ) ;
      TH1F* h_chi2_energy_notmu = get_hist( "h_chi2_energy_notmu" ) ;

      h_chi2_energy_ecal_mu -> SetLineColor(4) ;
      h_chi2_energy_hcal_mu -> SetLineColor(4) ;
      h_chi2_energy_mu -> SetLineColor(4) ;

      h_chi2_energy_ecal_notmu -> SetLineColor(2) ;
      h_chi2_energy_hcal_notmu -> SetLineColor(2) ;
      h_chi2_energy_notmu -> SetLineColor(2) ;


      can -> cd( ci++ ) ;
      h_chi2_energy_ecal_mu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_energy_hcal_mu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_energy_mu -> Draw() ;
      gPad -> SetGridx(1) ;


      can -> cd( ci++ ) ;
      h_chi2_energy_ecal_notmu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_energy_hcal_notmu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_energy_notmu -> Draw() ;
      gPad -> SetGridx(1) ;


      can -> Update() ;
      can -> Draw() ;
      gSystem -> ProcessEvents() ;


    //------------------------------------

      nx = 3 ;
      ny = 2 ;

      gStyle -> SetOptStat("emr") ;

      canx0 += can_dxy ;
      cany0 += can_dxy ;

      TCanvas* can1b = get_canvas( "can1b", "chi2 plots, isolation", canx0, cany0, nx*can_w_per_plot, ny*can_h_per_plot ) ;

      can = can1b ;

      ci = 1 ;

      can -> Clear() ;
      can -> cd() ;
      can -> Divide( 3, ny ) ;


      TH1F* h_chi2_iso_ecal_mu = get_hist( "h_chi2_iso_ecal_mu" ) ;
      TH1F* h_chi2_iso_hcal_mu = get_hist( "h_chi2_iso_hcal_mu" ) ;
      TH1F* h_chi2_iso_mu = get_hist( "h_chi2_iso_mu" ) ;

      TH1F* h_chi2_iso_ecal_notmu = get_hist( "h_chi2_iso_ecal_notmu" ) ;
      TH1F* h_chi2_iso_hcal_notmu = get_hist( "h_chi2_iso_hcal_notmu" ) ;
      TH1F* h_chi2_iso_notmu = get_hist( "h_chi2_iso_notmu" ) ;

      h_chi2_iso_ecal_mu -> SetLineColor(4) ;
      h_chi2_iso_hcal_mu -> SetLineColor(4) ;
      h_chi2_iso_mu -> SetLineColor(4) ;

      h_chi2_iso_ecal_notmu -> SetLineColor(2) ;
      h_chi2_iso_hcal_notmu -> SetLineColor(2) ;
      h_chi2_iso_notmu -> SetLineColor(2) ;


      can -> cd( ci++ ) ;
      h_chi2_iso_ecal_mu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_iso_hcal_mu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_iso_mu -> Draw() ;
      gPad -> SetGridx(1) ;


      can -> cd( ci++ ) ;
      h_chi2_iso_ecal_notmu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_iso_hcal_notmu -> Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_chi2_iso_notmu -> Draw() ;
      gPad -> SetGridx(1) ;


      can -> Update() ;
      can -> Draw() ;
      gSystem -> ProcessEvents() ;


    //------------------------------------



    //------------------------------------

      nx = 4 ;
      ny = 2 ;

      gStyle -> SetOptStat("emr") ;

      canx0 += can_dxy ;
      cany0 += can_dxy ;

      TCanvas* can2 = get_canvas( "can2", "Track pt, eta, Nhit", canx0, cany0, nx*can_w_per_plot, ny*can_h_per_plot ) ;


      can = can2 ;

      can -> Clear() ;
      can -> cd() ;
      can -> Divide( nx, ny ) ;

      ci = 1 ;

      TH1F* h_rc_pt_mu = get_hist( "h_rc_pt_mu" ) ;
      TH1F* h_rc_eta_mu = get_hist( "h_rc_eta_mu" ) ;

      TH1F* h_rc_pt_notmu = get_hist( "h_rc_pt_notmu" ) ;
      TH1F* h_rc_eta_notmu = get_hist( "h_rc_eta_notmu" ) ;

      TH2F* h_nhits_ecal_vs_eta_mu = get_hist2d( "h_nhits_ecal_vs_eta_mu" ) ;
      TH2F* h_nhits_hcal_vs_eta_mu = get_hist2d( "h_nhits_hcal_vs_eta_mu" ) ;

      TH2F* h_nhits_ecal_vs_eta_notmu = get_hist2d( "h_nhits_ecal_vs_eta_notmu" ) ;
      TH2F* h_nhits_hcal_vs_eta_notmu = get_hist2d( "h_nhits_hcal_vs_eta_notmu" ) ;

      can -> cd( ci++ ) ;
      h_rc_pt_mu -> Draw() ;
      gPad -> SetLogy(1) ;

      can -> cd( ci++ ) ;
      h_rc_eta_mu -> Draw() ;

      can -> cd( ci++ ) ;
      h_nhits_ecal_vs_eta_mu -> Draw("box" ) ;

      can -> cd( ci++ ) ;
      h_nhits_hcal_vs_eta_mu -> Draw("box" ) ;




      can -> cd( ci++ ) ;
      h_rc_pt_notmu -> Draw() ;
      gPad -> SetLogy(1) ;

      can -> cd( ci++ ) ;
      h_rc_eta_notmu -> Draw() ;

      can -> cd( ci++ ) ;
      h_nhits_ecal_vs_eta_notmu -> Draw("box" ) ;

      can -> cd( ci++ ) ;
      h_nhits_hcal_vs_eta_notmu -> Draw("box" ) ;









    //------------------------------------



      nx = 2 ;
      ny = 2 ;

      gStyle -> SetOptStat("emr") ;

      canx0 += can_dxy ;
      cany0 += can_dxy ;

      TCanvas* can3 = get_canvas( "can3", "Energy", canx0, cany0, nx*can_w_per_plot, ny*can_h_per_plot ) ;


      can = can3 ;

      can -> Clear() ;
      can -> cd() ;
      can -> Divide( nx, ny ) ;

      ci = 1 ;

      TH1F* h_energy_ecal_mu = get_hist( "h_energy_ecal_mu" ) ;
      TH1F* h_energy_ecal_notmu = get_hist( "h_energy_ecal_notmu" ) ;

      TH1F* h_energy_hcal_mu = get_hist( "h_energy_hcal_mu" ) ;
      TH1F* h_energy_hcal_notmu = get_hist( "h_energy_hcal_notmu" ) ;

      can -> cd( ci++ ) ;
      h_energy_ecal_mu->Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_energy_hcal_mu->Draw() ;
      gPad -> SetGridx(1) ;


      can -> cd( ci++ ) ;
      h_energy_ecal_notmu->Draw() ;
      gPad -> SetGridx(1) ;

      can -> cd( ci++ ) ;
      h_energy_hcal_notmu->Draw() ;
      gPad -> SetGridx(1) ;




    //------------------------------------



      nx = 2 ;
      ny = 2 ;

      gStyle -> SetOptStat("emr") ;

      canx0 += can_dxy ;
      cany0 += can_dxy ;

      TCanvas* can4 = get_canvas( "can4", "N-1, E/L plots", canx0, cany0, nx*can_w_per_plot, ny*can_h_per_plot ) ;


      can = can4 ;

      can -> Clear() ;
      can -> cd() ;
      can -> Divide( nx, ny ) ;

      ci = 1 ;

      TH1F* h_hcal_nm1_E_over_L_mu = get_hist( "h_hcal_nm1_E_over_L_mu" ) ;
      TH1F* h_hcal_nm1_E_over_L_notmu = get_hist( "h_hcal_nm1_E_over_L_notmu" ) ;

      TH1F* h_hcal_nm1_E_over_L_mu_isoloose = get_hist( "h_hcal_nm1_E_over_L_mu_isoloose" ) ;
      TH1F* h_hcal_nm1_E_over_L_notmu_isoloose = get_hist( "h_hcal_nm1_E_over_L_notmu_isoloose" ) ;

      TH1F* h_hcal_nm1_E_over_L_mu_isotight = get_hist( "h_hcal_nm1_E_over_L_mu_isotight" ) ;
      TH1F* h_hcal_nm1_E_over_L_notmu_isotight = get_hist( "h_hcal_nm1_E_over_L_notmu_isotight" ) ;

      TH1F* h_hcal_nm1_E_over_L_mu_eloose = get_hist( "h_hcal_nm1_E_over_L_mu_eloose" ) ;
      TH1F* h_hcal_nm1_E_over_L_notmu_eloose = get_hist( "h_hcal_nm1_E_over_L_notmu_eloose" ) ;

      TH1F* h_hcal_nm1_E_over_L_mu_etight = get_hist( "h_hcal_nm1_E_over_L_mu_etight" ) ;
      TH1F* h_hcal_nm1_E_over_L_notmu_etight = get_hist( "h_hcal_nm1_E_over_L_notmu_etight" ) ;

      h_hcal_nm1_E_over_L_mu_isoloose -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_mu_isotight -> SetLineColor(4) ;
      h_hcal_nm1_E_over_L_notmu_isoloose -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_notmu_isotight -> SetLineColor(4) ;
      h_hcal_nm1_E_over_L_mu_eloose -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_mu_etight -> SetLineColor(4) ;
      h_hcal_nm1_E_over_L_notmu_eloose -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_notmu_etight -> SetLineColor(4) ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu -> Draw() ;
      h_hcal_nm1_E_over_L_mu_isoloose -> Draw("same") ;
      h_hcal_nm1_E_over_L_mu_isotight -> Draw("same") ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu -> Draw() ;
      h_hcal_nm1_E_over_L_notmu_isoloose -> Draw("same") ;
      h_hcal_nm1_E_over_L_notmu_isotight -> Draw("same") ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu -> Draw() ;
      h_hcal_nm1_E_over_L_mu_eloose -> Draw("same") ;
      h_hcal_nm1_E_over_L_mu_etight -> Draw("same") ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu -> Draw() ;
      h_hcal_nm1_E_over_L_notmu_eloose -> Draw("same") ;
      h_hcal_nm1_E_over_L_notmu_etight -> Draw("same") ;


    //------------------------------------



      nx = 2 ;
      ny = 2 ;

      gStyle -> SetOptStat("emr") ;

      canx0 += can_dxy ;
      cany0 += can_dxy ;

      TCanvas* can5 = get_canvas( "can5", "N-1, E/L plots, 2", canx0, cany0, nx*can_w_per_plot, ny*can_h_per_plot ) ;


      can = can5 ;

      can -> Clear() ;
      can -> cd() ;
      can -> Divide( nx, ny ) ;

      ci = 1 ;

      TH1F* h_hcal_nm1_E_over_L_mu_isoloose_eloose = get_hist( "h_hcal_nm1_E_over_L_mu_isoloose_eloose" ) ;
      TH1F* h_hcal_nm1_E_over_L_notmu_isoloose_eloose = get_hist( "h_hcal_nm1_E_over_L_notmu_isoloose_eloose" ) ;

      TH1F* h_hcal_nm1_E_over_L_mu_isotight_etight = get_hist( "h_hcal_nm1_E_over_L_mu_isotight_etight" ) ;
      TH1F* h_hcal_nm1_E_over_L_notmu_isotight_etight = get_hist( "h_hcal_nm1_E_over_L_notmu_isotight_etight" ) ;

      h_hcal_nm1_E_over_L_mu_isoloose_eloose -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_mu_isotight_etight -> SetLineColor(4) ;
      h_hcal_nm1_E_over_L_notmu_isoloose_eloose -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_notmu_isotight_etight -> SetLineColor(4) ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu -> Draw() ;
      h_hcal_nm1_E_over_L_mu_isoloose_eloose -> Draw("same") ;
      h_hcal_nm1_E_over_L_mu_isotight_etight -> Draw("same") ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu -> Draw() ;
      h_hcal_nm1_E_over_L_notmu_isoloose_eloose -> Draw("same") ;
      h_hcal_nm1_E_over_L_notmu_isotight_etight -> Draw("same") ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu_isotight_etight -> Draw() ;
      h_hcal_nm1_E_over_L_mu -> Draw("same") ;
      h_hcal_nm1_E_over_L_mu_isoloose_eloose -> Draw("same") ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu_isotight_etight -> Draw() ;
      h_hcal_nm1_E_over_L_notmu -> Draw("same") ;
      h_hcal_nm1_E_over_L_notmu_isoloose_eloose -> Draw("same") ;


    //---------

      TCanvas* can6 = get_canvas( "can6", "HCAL E vs L", 50, 50, 1400, 1000 ) ;

      can = can6 ;

      h_hcal_nm1_E_over_L_mu -> SetLineColor(4) ;

      float mumax = h_hcal_nm1_E_over_L_mu -> GetMaximum() ;
      h_hcal_nm1_E_over_L_mu_isotight -> SetMaximum( mumax ) ;
      h_hcal_nm1_E_over_L_mu_etight -> SetMaximum( mumax ) ;
      h_hcal_nm1_E_over_L_mu_isotight_etight -> SetMaximum( mumax ) ;

      float notmumax = h_hcal_nm1_E_over_L_notmu -> GetMaximum() ;
      h_hcal_nm1_E_over_L_notmu_isotight -> SetMaximum( notmumax ) ;
      h_hcal_nm1_E_over_L_notmu_etight -> SetMaximum( notmumax ) ;
      h_hcal_nm1_E_over_L_notmu_isotight_etight -> SetMaximum( notmumax ) ;

      h_hcal_nm1_E_over_L_notmu -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_notmu_isotight -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_notmu_etight -> SetLineColor(2) ;
      h_hcal_nm1_E_over_L_notmu_isotight_etight -> SetLineColor(2) ;

      h_hcal_nm1_E_over_L_mu -> SetFillColor( kBlue-10 ) ;
      h_hcal_nm1_E_over_L_mu_isotight -> SetFillColor( kBlue-10 ) ;
      h_hcal_nm1_E_over_L_mu_etight -> SetFillColor( kBlue-10 ) ;
      h_hcal_nm1_E_over_L_mu_isotight_etight -> SetFillColor( kBlue-10 ) ;

      h_hcal_nm1_E_over_L_notmu -> SetFillColor( kRed-10 ) ;
      h_hcal_nm1_E_over_L_notmu_isotight -> SetFillColor( kRed-10 ) ;
      h_hcal_nm1_E_over_L_notmu_etight -> SetFillColor( kRed-10 ) ;
      h_hcal_nm1_E_over_L_notmu_isotight_etight -> SetFillColor( kRed-10 ) ;



      can -> Clear() ;
      can -> Divide(4,2) ;

      ci = 1 ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu -> Draw() ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu_isotight -> Draw() ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu_etight -> Draw() ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_mu_isotight_etight -> Draw() ;


      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu -> Draw() ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu_isotight -> Draw() ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu_etight -> Draw() ;

      can -> cd( ci++ ) ;
      h_hcal_nm1_E_over_L_notmu_isotight_etight -> Draw() ;


    //---------

      nx = 1 ;
      ny = 2 ;

      gStyle -> SetOptStat(0) ;

      canx0 += can_dxy ;
      cany0 += can_dxy ;

      TCanvas* can7 = get_canvas( "can7", "chi2 energy vs iso", canx0, cany0, nx*can_w_per_plot, ny*can_h_per_plot ) ;


      can = can7 ;

      can -> Clear() ;
      can -> cd() ;
      can -> Divide( nx, ny ) ;

      ci = 1 ;

      TH2F* h_chi2_energy_vs_iso_mu = get_hist2d( "h_chi2_energy_vs_iso_mu" ) ;
      TH2F* h_chi2_energy_vs_iso_notmu = get_hist2d( "h_chi2_energy_vs_iso_notmu" ) ;
      h_chi2_energy_vs_iso_mu -> SetXTitle( "log10( chi2, iso )" ) ;
      h_chi2_energy_vs_iso_mu -> SetYTitle( "log10( chi2, energy )" ) ;
      h_chi2_energy_vs_iso_notmu -> SetXTitle( "log10( chi2, iso )" ) ;
      h_chi2_energy_vs_iso_notmu -> SetYTitle( "log10( chi2, energy )" ) ;

      can -> cd( ci++ ) ;
      h_chi2_energy_vs_iso_mu -> Draw("colz") ;
      grid_on() ;

      can -> cd( ci++ ) ;
      h_chi2_energy_vs_iso_notmu -> Draw("colz") ;
      grid_on() ;


   } // draw_plots1






