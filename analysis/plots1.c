#define plots1_cxx
#include "plots1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "histio.c"
#include "utils.c"
#include "draw_plots1.c"

void plots1::Loop()
{
   if (fChain == 0) return;

 //float iso_loose = 1e3 ;
 //float iso_tight = 1e2 ;
 //float energy_loose = 1e2 ;
 //float energy_tight = 1e1 ;

   float iso_loose = 1e3 ;
   float iso_tight = 1e1 ;
   float energy_loose = 1e2 ;
   float energy_tight = 1e1 ;
   float minpt = 0.5 ;

   Long64_t nentries = fChain->GetEntries();
   printf("\n\n Chain has %llu entries.\n\n", nentries ) ;

   gDirectory -> Delete( "h*" ) ;

   int bins = 60 ;
   float xlow = -3. ;
   float xhigh = 8. ;

   TH1F* h_chi2_energy_ecal_all = new TH1F( "h_chi2_energy_ecal_all", "log10 chi2 energy, ecal, N - 1, all", bins, xlow, xhigh ) ;
   TH1F* h_chi2_energy_hcal_all = new TH1F( "h_chi2_energy_hcal_all", "log10 chi2 energy, hcal, N - 1, all", bins, xlow, xhigh ) ;
   TH1F* h_chi2_energy_all = new TH1F( "h_chi2_energy_all", "log10 chi2 energy, ecal+hcal, N - 1, all", bins, xlow, xhigh ) ;

   TH1F* h_chi2_energy_ecal_mu = new TH1F( "h_chi2_energy_ecal_mu", "log10 chi2 energy, ecal, N - 1, muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_energy_hcal_mu = new TH1F( "h_chi2_energy_hcal_mu", "log10 chi2 energy, hcal, N - 1, muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_energy_mu = new TH1F( "h_chi2_energy_mu", "log10 chi2 energy, ecal+hcal, N - 1, muon", bins, xlow, xhigh ) ;

   TH1F* h_chi2_energy_ecal_notmu = new TH1F( "h_chi2_energy_ecal_notmu", "log10 chi2 energy, ecal, N - 1, not muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_energy_hcal_notmu = new TH1F( "h_chi2_energy_hcal_notmu", "log10 chi2 energy, hcal, N - 1, not muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_energy_notmu = new TH1F( "h_chi2_energy_notmu", "log10 chi2 energy, ecal+hcal, N - 1, not muon", bins, xlow, xhigh ) ;

   TH1F* h_chi2_iso_ecal_all = new TH1F( "h_chi2_iso_ecal_all", "log10 chi2 iso, ecal, N - 1, all", bins, xlow, xhigh ) ;
   TH1F* h_chi2_iso_hcal_all = new TH1F( "h_chi2_iso_hcal_all", "log10 chi2 iso, hcal, N - 1, all", bins, xlow, xhigh ) ;
   TH1F* h_chi2_iso_all = new TH1F( "h_chi2_iso_all", "log10 chi2 iso, ecal+hcal, N - 1, all", bins, xlow, xhigh ) ;

   TH1F* h_chi2_iso_ecal_mu = new TH1F( "h_chi2_iso_ecal_mu", "log10 chi2 iso, ecal, N - 1, muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_iso_hcal_mu = new TH1F( "h_chi2_iso_hcal_mu", "log10 chi2 iso, hcal, N - 1, muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_iso_mu = new TH1F( "h_chi2_iso_mu", "log10 chi2 iso, ecal+hcal, N - 1, muon", bins, xlow, xhigh ) ;

   TH1F* h_chi2_iso_ecal_notmu = new TH1F( "h_chi2_iso_ecal_notmu", "log10 chi2 iso, ecal, N - 1, not muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_iso_hcal_notmu = new TH1F( "h_chi2_iso_hcal_notmu", "log10 chi2 iso, hcal, N - 1, not muon", bins, xlow, xhigh ) ;
   TH1F* h_chi2_iso_notmu = new TH1F( "h_chi2_iso_notmu", "log10 chi2 iso, ecal+hcal, N - 1, not muon", bins, xlow, xhigh ) ;


   TH2F* h_chi2_energy_vs_iso_mu = new TH2F( "h_chi2_energy_vs_iso_mu", "log10 chi2, energy vs iso, muon", bins, xlow, xhigh, bins, xlow, xhigh ) ;
   TH2F* h_chi2_energy_vs_iso_notmu = new TH2F( "h_chi2_energy_vs_iso_notmu", "log10 chi2, energy vs iso, not muon", bins, xlow, xhigh, bins, xlow, xhigh ) ;


   xlow = -0.1 ;
   xhigh = 2. ;
   TH1F* h_energy_ecal_all = new TH1F( "h_energy_ecal_all", "Energy on track, ecal, all", bins, xlow, xhigh ) ;
   TH1F* h_energy_ecal_mu = new TH1F( "h_energy_ecal_mu", "Energy on track, ecal, muon", bins, xlow, xhigh ) ;
   TH1F* h_energy_ecal_notmu = new TH1F( "h_energy_ecal_notmu", "Energy on track, ecal, not muon", bins, xlow, xhigh ) ;

   TH1F* h_energy_ecal_all_isoloose = new TH1F( "h_energy_ecal_all_isoloose", "Energy on track, ecal, all, loose iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_ecal_mu_isoloose = new TH1F( "h_energy_ecal_mu_isoloose", "Energy on track, ecal, muon, loose iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_ecal_notmu_isoloose = new TH1F( "h_energy_ecal_notmu_isoloose", "Energy on track, ecal, not muon, loose iso", bins, xlow, xhigh ) ;

   TH1F* h_energy_ecal_all_isotight = new TH1F( "h_energy_ecal_all_isotight", "Energy on track, ecal, all, tight iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_ecal_mu_isotight = new TH1F( "h_energy_ecal_mu_isotight", "Energy on track, ecal, muon, tight iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_ecal_notmu_isotight = new TH1F( "h_energy_ecal_notmu_isotight", "Energy on track, ecal, not muon, tight iso", bins, xlow, xhigh ) ;

   xlow = -0.2 ;
   xhigh = 5. ;
   TH1F* h_energy_hcal_all = new TH1F( "h_energy_hcal_all", "Energy on track, hcal, all", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_mu = new TH1F( "h_energy_hcal_mu", "Energy on track, hcal, muon", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_notmu = new TH1F( "h_energy_hcal_notmu", "Energy on track, hcal, not muon", bins, xlow, xhigh ) ;

   TH1F* h_energy_hcal_all_isoloose = new TH1F( "h_energy_hcal_all_isoloose", "Energy on track, hcal, all, loose iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_mu_isoloose = new TH1F( "h_energy_hcal_mu_isoloose", "Energy on track, hcal, muon, loose iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_notmu_isoloose = new TH1F( "h_energy_hcal_notmu_isoloose", "Energy on track, hcal, not muon, loose iso", bins, xlow, xhigh ) ;

   TH1F* h_energy_hcal_all_isotight = new TH1F( "h_energy_hcal_all_isotight", "Energy on track, hcal, all, tight iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_mu_isotight = new TH1F( "h_energy_hcal_mu_isotight", "Energy on track, hcal, muon, tight iso", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_notmu_isotight = new TH1F( "h_energy_hcal_notmu_isotight", "Energy on track, hcal, not muon, tight iso", bins, xlow, xhigh ) ;

   TH1F* h_energy_hcal_ecalecut_all = new TH1F( "h_energy_hcal_ecalecut_all", "Energy on track, hcal, ecal E cut, all", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_mu = new TH1F( "h_energy_hcal_ecalecut_mu", "Energy on track, hcal, ecal E cut, muon", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_notmu = new TH1F( "h_energy_hcal_ecalecut_notmu", "Energy on track, hcal, ecal E cut, not muon", bins, xlow, xhigh ) ;

   TH1F* h_energy_hcal_ecalecut_ecalisocut_all = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_all", "Energy on track, hcal, ecal E, iso cuts, all", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_ecalisocut_mu = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_mu", "Energy on track, hcal, ecal E, iso cuts, muon", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_ecalisocut_notmu = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_notmu", "Energy on track, hcal, ecal E, iso cuts, not muon", bins, xlow, xhigh ) ;

   TH1F* h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_all = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_all", "Energy on track, hcal, ecal E, iso, chi2 cuts, all", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_mu = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_mu", "Energy on track, hcal, ecal E, iso, chi2 cuts, muon", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_notmu = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_notmu", "Energy on track, hcal, ecal E, iso, chi2 cuts, not muon", bins, xlow, xhigh ) ;

   TH1F* h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_all = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_all", "Energy on track, hcal, ecal E, iso, chi2 cuts, hcal iso cut, all", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_mu = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_mu", "Energy on track, hcal, ecal E, iso, chi2 cuts, hcal iso cut, muon", bins, xlow, xhigh ) ;
   TH1F* h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_notmu = new TH1F( "h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_notmu", "Energy on track, hcal, ecal E, iso, chi2 cuts, hcal iso cut, not muon", bins, xlow, xhigh ) ;

   TH1F* h_mcpid_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_notmu = new TH1F( "h_mcpid_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_notmu", "MC pid, ecal E, iso, chi2 cuts, hcal iso cut, not muon", 41, -0.5, 40.5 ) ;

   TH1F* h_rc_eta_all = new TH1F( "h_rc_eta_all", "track eta, all", bins, 2., 5. ) ;
   TH1F* h_rc_eta_mu = new TH1F( "h_rc_eta_mu", "track eta, muon", bins, 2., 5. ) ;
   TH1F* h_rc_eta_notmu = new TH1F( "h_rc_eta_notmu", "track eta, not muon", bins, 2., 5. ) ;

   TH1F* h_rc_pt_all = new TH1F( "h_rc_pt_all", "track pt, all", bins, -2., 15. ) ;
   TH1F* h_rc_pt_mu = new TH1F( "h_rc_pt_mu", "track pt, muon", bins, -2., 15. ) ;
   TH1F* h_rc_pt_notmu = new TH1F( "h_rc_pt_notmu", "track pt, not muon", bins, -2., 15. ) ;

   TH2F* h_nhits_ecal_vs_eta_all = new TH2F( "h_nhits_ecal_vs_eta_all", "Nhits vs eta, ecal, all", 40, 2., 5., 5, -0.5, 4.5 )  ;
   TH2F* h_nhits_hcal_vs_eta_all = new TH2F( "h_nhits_hcal_vs_eta_all", "Nhits vs eta, hcal, all", 40, 2., 5., 5, -0.5, 4.5 )  ;

   TH2F* h_nhits_ecal_vs_eta_mu = new TH2F( "h_nhits_ecal_vs_eta_mu", "Nhits vs eta, ecal, muon", 40, 2., 5., 5, -0.5, 4.5 )  ;
   TH2F* h_nhits_hcal_vs_eta_mu = new TH2F( "h_nhits_hcal_vs_eta_mu", "Nhits vs eta, hcal, muon", 40, 2., 5., 5, -0.5, 4.5 )  ;

   TH2F* h_nhits_ecal_vs_eta_notmu = new TH2F( "h_nhits_ecal_vs_eta_notmu", "Nhits vs eta, ecal, not muon", 40, 2., 5., 5, -0.5, 4.5 )  ;
   TH2F* h_nhits_hcal_vs_eta_notmu = new TH2F( "h_nhits_hcal_vs_eta_notmu", "Nhits vs eta, hcal, not muon", 40, 2., 5., 5, -0.5, 4.5 )  ;

   bins = 20 ;
   TH2F* h_hcal_nm1_E_vs_L_all = new TH2F( "h_hcal_nm1_E_vs_L_all", "hcal N-1, E vs L, all", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_mu = new TH2F( "h_hcal_nm1_E_vs_L_mu", "hcal N-1, E vs L, muon", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_notmu = new TH2F( "h_hcal_nm1_E_vs_L_notmu", "hcal N-1, E vs L, not muon", bins, 0., 100., bins, 0., 2. ) ;

   TH2F* h_hcal_nm1_E_vs_L_all_isoloose = new TH2F( "h_hcal_nm1_E_vs_L_all_isoloose", "hcal N-1, E vs L, all, loose iso", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_mu_isoloose = new TH2F( "h_hcal_nm1_E_vs_L_mu_isoloose", "hcal N-1, E vs L, muon, loose iso", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_notmu_isoloose = new TH2F( "h_hcal_nm1_E_vs_L_notmu_isoloose", "hcal N-1, E vs L, not muon, loose iso", bins, 0., 100., bins, 0., 2. ) ;

   TH2F* h_hcal_nm1_E_vs_L_all_isotight = new TH2F( "h_hcal_nm1_E_vs_L_all_isotight", "hcal N-1, E vs L, all, tight iso", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_mu_isotight = new TH2F( "h_hcal_nm1_E_vs_L_mu_isotight", "hcal N-1, E vs L, muon, tight iso", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_notmu_isotight = new TH2F( "h_hcal_nm1_E_vs_L_notmu_isotight", "hcal N-1, E vs L, not muon, tight iso", bins, 0., 100., bins, 0., 2. ) ;



   TH2F* h_hcal_nm1_E_vs_L_ecalecut_all = new TH2F( "h_hcal_nm1_E_vs_L_ecalecut_all", "hcal N-1, E vs L, ecal E cut, all", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_ecalecut_mu = new TH2F( "h_hcal_nm1_E_vs_L_ecalecut_mu", "hcal N-1, E vs L, ecal E cut, muon", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_ecalecut_notmu = new TH2F( "h_hcal_nm1_E_vs_L_ecalecut_notmu", "hcal N-1, E vs L, ecal E cut, not muon", bins, 0., 100., bins, 0., 2. ) ;

   TH2F* h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_all = new TH2F( "h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_all", "hcal N-1, E vs L, ecal E, iso cuts, all", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_mu = new TH2F( "h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_mu", "hcal N-1, E vs L, ecal E, iso cuts, muon", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_notmu = new TH2F( "h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_notmu", "hcal N-1, E vs L, ecal E, iso cuts, not muon", bins, 0., 100., bins, 0., 2. ) ;




   TH2F* h_hcal_nm1_E_vs_L_all_isoloose_eloose = new TH2F( "h_hcal_nm1_E_vs_L_all_isoloose_eloose", "hcal N-1, E vs L, all, loose iso, loose energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_mu_isoloose_eloose = new TH2F( "h_hcal_nm1_E_vs_L_mu_isoloose_eloose", "hcal N-1, E vs L, muon, loose iso, loose energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_notmu_isoloose_eloose = new TH2F( "h_hcal_nm1_E_vs_L_notmu_isoloose_eloose", "hcal N-1, E vs L, not muon, loose iso, loose energy chi2", bins, 0., 100., bins, 0., 2. ) ;

   TH2F* h_hcal_nm1_E_vs_L_all_isotight_eloose = new TH2F( "h_hcal_nm1_E_vs_L_all_isotight_eloose", "hcal N-1, E vs L, all, tight iso, loose energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_mu_isotight_eloose = new TH2F( "h_hcal_nm1_E_vs_L_mu_isotight_eloose", "hcal N-1, E vs L, muon, tight iso, loose energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_notmu_isotight_eloose = new TH2F( "h_hcal_nm1_E_vs_L_notmu_isotight_eloose", "hcal N-1, E vs L, not muon, tight iso, loose energy chi2", bins, 0., 100., bins, 0., 2. ) ;



   TH2F* h_hcal_nm1_E_vs_L_all_isoloose_etight = new TH2F( "h_hcal_nm1_E_vs_L_all_isoloose_etight", "hcal N-1, E vs L, all, loose iso, tight energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_mu_isoloose_etight = new TH2F( "h_hcal_nm1_E_vs_L_mu_isoloose_etight", "hcal N-1, E vs L, muon, loose iso, tight energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_notmu_isoloose_etight = new TH2F( "h_hcal_nm1_E_vs_L_notmu_isoloose_etight", "hcal N-1, E vs L, not muon, loose iso, tight energy chi2", bins, 0., 100., bins, 0., 2. ) ;

   TH2F* h_hcal_nm1_E_vs_L_all_isotight_etight = new TH2F( "h_hcal_nm1_E_vs_L_all_isotight_etight", "hcal N-1, E vs L, all, tight iso, tight energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_mu_isotight_etight = new TH2F( "h_hcal_nm1_E_vs_L_mu_isotight_etight", "hcal N-1, E vs L, muon, tight iso, tight energy chi2", bins, 0., 100., bins, 0., 2. ) ;
   TH2F* h_hcal_nm1_E_vs_L_notmu_isotight_etight = new TH2F( "h_hcal_nm1_E_vs_L_notmu_isotight_etight", "hcal N-1, E vs L, not muon, tight iso, tight energy chi2", bins, 0., 100., bins, 0., 2. ) ;




   bins = 80 ;

   TH1F* h_hcal_nm1_E_over_L_all = new TH1F( "h_hcal_nm1_E_over_L_all", "hcal N-1, E/L, all", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu = new TH1F( "h_hcal_nm1_E_over_L_mu", "hcal N-1, E/L, muon", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu = new TH1F( "h_hcal_nm1_E_over_L_notmu", "hcal N-1, E/L, not muon", bins, -0.01, 0.1 ) ;

   TH1F* h_hcal_nm1_E_over_L_all_isoloose = new TH1F( "h_hcal_nm1_E_over_L_all_isoloose", "hcal N-1, E/L, all, loose iso", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_isoloose = new TH1F( "h_hcal_nm1_E_over_L_mu_isoloose", "hcal N-1, E/L, muon, loose iso", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_isoloose = new TH1F( "h_hcal_nm1_E_over_L_notmu_isoloose", "hcal N-1, E/L, not muon, loose iso", bins, -0.01, 0.1 ) ;

   TH1F* h_hcal_nm1_E_over_L_all_isotight = new TH1F( "h_hcal_nm1_E_over_L_all_isotight", "hcal N-1, E/L, all, tight iso", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_isotight = new TH1F( "h_hcal_nm1_E_over_L_mu_isotight", "hcal N-1, E/L, muon, tight iso", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_isotight = new TH1F( "h_hcal_nm1_E_over_L_notmu_isotight", "hcal N-1, E/L, not muon, tight iso", bins, -0.01, 0.1 ) ;


   TH1F* h_hcal_nm1_E_over_L_all_eloose = new TH1F( "h_hcal_nm1_E_over_L_all_eloose", "hcal N-1, E/L, all, loose e", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_eloose = new TH1F( "h_hcal_nm1_E_over_L_mu_eloose", "hcal N-1, E/L, muon, loose e", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_eloose = new TH1F( "h_hcal_nm1_E_over_L_notmu_eloose", "hcal N-1, E/L, not muon, loose e", bins, -0.01, 0.1 ) ;

   TH1F* h_hcal_nm1_E_over_L_all_etight = new TH1F( "h_hcal_nm1_E_over_L_all_etight", "hcal N-1, E/L, all, tight e", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_etight = new TH1F( "h_hcal_nm1_E_over_L_mu_etight", "hcal N-1, E/L, muon, tight e", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_etight = new TH1F( "h_hcal_nm1_E_over_L_notmu_etight", "hcal N-1, E/L, not muon, tight e", bins, -0.01, 0.1 ) ;




   TH1F* h_hcal_nm1_E_over_L_all_isoloose_eloose = new TH1F( "h_hcal_nm1_E_over_L_all_isoloose_eloose", "hcal N-1, E/L, all, loose iso, loose energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_isoloose_eloose = new TH1F( "h_hcal_nm1_E_over_L_mu_isoloose_eloose", "hcal N-1, E/L, muon, loose iso, loose energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_isoloose_eloose = new TH1F( "h_hcal_nm1_E_over_L_notmu_isoloose_eloose", "hcal N-1, E/L, not muon, loose iso, loose energy chi2", bins, -0.01, 0.1 ) ;

   TH1F* h_hcal_nm1_E_over_L_all_isoloose_etight = new TH1F( "h_hcal_nm1_E_over_L_all_isoloose_etight", "hcal N-1, E/L, all, loose iso, tight energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_isoloose_etight = new TH1F( "h_hcal_nm1_E_over_L_mu_isoloose_etight", "hcal N-1, E/L, muon, loose iso, tight energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_isoloose_etight = new TH1F( "h_hcal_nm1_E_over_L_notmu_isoloose_etight", "hcal N-1, E/L, not muon, loose iso, tight energy chi2", bins, -0.01, 0.1 ) ;



   TH1F* h_hcal_nm1_E_over_L_all_isotight_eloose = new TH1F( "h_hcal_nm1_E_over_L_all_isotight_eloose", "hcal N-1, E/L, all, tight iso, loose energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_isotight_eloose = new TH1F( "h_hcal_nm1_E_over_L_mu_isotight_eloose", "hcal N-1, E/L, muon, tight iso, loose energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_isotight_eloose = new TH1F( "h_hcal_nm1_E_over_L_notmu_isotight_eloose", "hcal N-1, E/L, not muon, tight iso, loose energy chi2", bins, -0.01, 0.1 ) ;

   TH1F* h_hcal_nm1_E_over_L_all_isotight_etight = new TH1F( "h_hcal_nm1_E_over_L_all_isotight_etight", "hcal N-1, E/L, all, tight iso, tight energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_mu_isotight_etight = new TH1F( "h_hcal_nm1_E_over_L_mu_isotight_etight", "hcal N-1, E/L, muon, tight iso, tight energy chi2", bins, -0.01, 0.1 ) ;
   TH1F* h_hcal_nm1_E_over_L_notmu_isotight_etight = new TH1F( "h_hcal_nm1_E_over_L_notmu_isotight_etight", "hcal N-1, E/L, not muon, tight iso, tight energy chi2", bins, -0.01, 0.1 ) ;



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if ( jentry % 1000 == 0 ) printf("  %9llu / %9llu\n", jentry, nentries ) ;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      for ( int rci=0; rci<rc_nm1_chi2_energy_ecal->size(); rci++ ) {

         int other_rci = rc_rci->at(rci) ;

         float eta = rcEta -> at( other_rci ) ;
         float pt = rcPt -> at( other_rci ) ;
         int mcind = rcTrackId -> at(other_rci) - 1 ;
         int mcpid(0) ;
         if ( mcind >= 0 ) mcpid = mcPid->at(mcind) ;

         if ( pt < minpt ) continue ;


         if ( (*rc_nm1_chi2_energy_ecal)[rci] < 0 ) continue ;
         if ( (*rc_nm1_chi2_energy_hcal)[rci] < 0 ) continue ;

         if ( rc_nm1_hcal_L1->at(rci) < 10 ) continue ;
         if ( rc_nm1_hcal_E1->at(rci) < 0.20 ) continue ;

         h_rc_eta_all -> Fill( eta ) ;
         h_rc_pt_all -> Fill( pt ) ;

         h_nhits_ecal_vs_eta_all -> Fill( eta, rc_nm1_nhit_ecal->at(rci) ) ;
         h_nhits_hcal_vs_eta_all -> Fill( eta, rc_nm1_nhit_hcal->at(rci) ) ;

         float nm1_chi2_ecal_energy = (*rc_nm1_chi2_energy_ecal)[rci] ;
         float nm1_chi2_hcal_energy = (*rc_nm1_chi2_energy_hcal)[rci] ;
         float nm1_chi2_energy = nm1_chi2_ecal_energy + nm1_chi2_hcal_energy ;

         h_chi2_energy_ecal_all -> Fill( log10( nm1_chi2_ecal_energy ) ) ;
         h_chi2_energy_hcal_all -> Fill( log10( nm1_chi2_hcal_energy ) ) ;
         h_chi2_energy_all      -> Fill( log10( nm1_chi2_energy ) ) ;

         float nm1_chi2_ecal_iso = (*rc_nm1_chi2_iso1_ecal)[rci] ;
         float nm1_chi2_hcal_iso = (*rc_nm1_chi2_iso1_hcal)[rci] ;
         float nm1_chi2_iso = nm1_chi2_ecal_iso + nm1_chi2_hcal_iso ;

         h_chi2_iso_ecal_all -> Fill( log10( nm1_chi2_ecal_iso ) ) ;
         h_chi2_iso_hcal_all -> Fill( log10( nm1_chi2_hcal_iso ) ) ;
         h_chi2_iso_all      -> Fill( log10( nm1_chi2_iso ) ) ;

         h_energy_ecal_all -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
         h_energy_hcal_all -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;



         float nm1_E = rc_nm1_hcal_E1->at(rci) ;
         float nm1_L = rc_nm1_hcal_L1->at(rci) ;

         h_hcal_nm1_E_vs_L_all -> Fill( nm1_L, nm1_E ) ;
         h_hcal_nm1_E_over_L_all -> Fill( nm1_E / nm1_L ) ;

         if ( (*rc_nm1_energy_ecal)[rci] > 0.23 && (*rc_nm1_energy_ecal)[rci] < 0.34 ) {
            h_energy_hcal_ecalecut_all -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
            h_hcal_nm1_E_vs_L_ecalecut_all -> Fill( nm1_L, nm1_E ) ;
            if ( (*rc_nm1_chi2_iso1_ecal)[rci] < 1e0 ) {
               h_energy_hcal_ecalecut_ecalisocut_all -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
               h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_all -> Fill( nm1_L, nm1_E ) ;
               if ( nm1_chi2_ecal_energy < 1e2 ) {
                  h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_all -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                  if ( (*rc_nm1_chi2_iso1_hcal)[rci] < 1e0 ) {
                     h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_all -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                  }
               }
            }
         }





         bool ismu = (*rc_ismu)[rci] ;


         if ( ismu ) {

            h_chi2_energy_ecal_mu -> Fill( log10( (*rc_nm1_chi2_energy_ecal)[rci] ) ) ;
            h_chi2_energy_hcal_mu -> Fill( log10( (*rc_nm1_chi2_energy_hcal)[rci] ) ) ;
            h_chi2_energy_mu      -> Fill( log10( (*rc_nm1_chi2_energy_ecal)[rci] + (*rc_nm1_chi2_energy_hcal)[rci] ) ) ;

            h_chi2_iso_ecal_mu -> Fill( log10( nm1_chi2_ecal_iso ) ) ;
            h_chi2_iso_hcal_mu -> Fill( log10( nm1_chi2_hcal_iso ) ) ;
            h_chi2_iso_mu      -> Fill( log10( nm1_chi2_iso ) ) ;


            h_chi2_energy_vs_iso_mu -> Fill( log10( (*rc_nm1_chi2_iso1_ecal)[rci] + (*rc_nm1_chi2_iso1_hcal)[rci] ) ,
                                             log10( (*rc_nm1_chi2_energy_ecal)[rci] + (*rc_nm1_chi2_energy_hcal)[rci] ) ) ;

            h_energy_ecal_mu -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
            h_energy_hcal_mu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;

            if ( (*rc_nm1_energy_ecal)[rci] > 0.23 && (*rc_nm1_energy_ecal)[rci] < 0.34 ) {
               h_energy_hcal_ecalecut_mu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
               h_hcal_nm1_E_vs_L_ecalecut_mu -> Fill( nm1_L, nm1_E ) ;
               if ( (*rc_nm1_chi2_iso1_ecal)[rci] < 1e0 ) {
                  h_energy_hcal_ecalecut_ecalisocut_mu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                  h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_mu -> Fill( nm1_L, nm1_E ) ;
                  if ( nm1_chi2_ecal_energy < 1e2 ) {
                     h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_mu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                     if ( (*rc_nm1_chi2_iso1_hcal)[rci] < 1e0 ) {
                        h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_mu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                     }
                  }
               }
            }

            h_rc_eta_mu -> Fill( eta ) ;
            h_rc_pt_mu -> Fill( pt ) ;
            h_nhits_ecal_vs_eta_mu -> Fill( eta, rc_nm1_nhit_ecal->at(rci) ) ;
            h_nhits_hcal_vs_eta_mu -> Fill( eta, rc_nm1_nhit_hcal->at(rci) ) ;

            h_hcal_nm1_E_vs_L_mu -> Fill( nm1_L, nm1_E ) ;
            h_hcal_nm1_E_over_L_mu -> Fill( nm1_E / nm1_L ) ;

         } else {

            h_chi2_energy_ecal_notmu -> Fill( log10( (*rc_nm1_chi2_energy_ecal)[rci] ) ) ;
            h_chi2_energy_hcal_notmu -> Fill( log10( (*rc_nm1_chi2_energy_hcal)[rci] ) ) ;
            h_chi2_energy_notmu      -> Fill( log10( (*rc_nm1_chi2_energy_ecal)[rci] + (*rc_nm1_chi2_energy_hcal)[rci] ) ) ;

            h_chi2_iso_ecal_notmu -> Fill( log10( nm1_chi2_ecal_iso ) ) ;
            h_chi2_iso_hcal_notmu -> Fill( log10( nm1_chi2_hcal_iso ) ) ;
            h_chi2_iso_notmu      -> Fill( log10( nm1_chi2_iso ) ) ;


            h_chi2_energy_vs_iso_notmu -> Fill( log10( (*rc_nm1_chi2_iso1_ecal)[rci] + (*rc_nm1_chi2_iso1_hcal)[rci] ) ,
                                                log10( (*rc_nm1_chi2_energy_ecal)[rci] + (*rc_nm1_chi2_energy_hcal)[rci] ) ) ;

            h_energy_ecal_notmu -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
            h_energy_hcal_notmu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;

            if ( (*rc_nm1_energy_ecal)[rci] > 0.23 && (*rc_nm1_energy_ecal)[rci] < 0.34 ) {
               h_energy_hcal_ecalecut_notmu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
               h_hcal_nm1_E_vs_L_ecalecut_notmu -> Fill( nm1_L, nm1_E ) ;
               if ( (*rc_nm1_chi2_iso1_ecal)[rci] < 1e0 ) {
                  h_energy_hcal_ecalecut_ecalisocut_notmu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                  h_hcal_nm1_E_vs_L_ecalecut_ecalisocut_notmu -> Fill( nm1_L, nm1_E ) ;
                  if ( nm1_chi2_ecal_energy < 1e2 ) {
                     h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_notmu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                     if ( (*rc_nm1_chi2_iso1_hcal)[rci] < 1e0 ) {
                        h_energy_hcal_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_notmu -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
                        h_mcpid_ecalecut_ecalisocut_ecalechi2cut_hcalisocut_notmu -> Fill( mcpid ) ;

                     }
                  }
               }
            }

            h_rc_eta_notmu -> Fill( eta ) ;
            h_rc_pt_notmu -> Fill( pt ) ;
            h_nhits_ecal_vs_eta_notmu -> Fill( eta, rc_nm1_nhit_ecal->at(rci) ) ;
            h_nhits_hcal_vs_eta_notmu -> Fill( eta, rc_nm1_nhit_hcal->at(rci) ) ;

            h_hcal_nm1_E_vs_L_notmu -> Fill( nm1_L, nm1_E ) ;
            h_hcal_nm1_E_over_L_notmu -> Fill( nm1_E / nm1_L ) ;
         }

        //--- loose energy chi2 cut
         if ( nm1_chi2_energy < energy_loose ) {
            h_hcal_nm1_E_over_L_all_eloose -> Fill( nm1_E / nm1_L ) ;
            if ( ismu ) {
               h_hcal_nm1_E_over_L_mu_eloose -> Fill( nm1_E / nm1_L ) ;
            } else {
               h_hcal_nm1_E_over_L_notmu_eloose -> Fill( nm1_E / nm1_L ) ;
            }
         }

        //--- tight energy chi2 cut
         if ( nm1_chi2_energy < energy_tight ) {
            h_hcal_nm1_E_over_L_all_etight -> Fill( nm1_E / nm1_L ) ;
            if ( ismu ) {
               h_hcal_nm1_E_over_L_mu_etight -> Fill( nm1_E / nm1_L ) ;
            } else {
               h_hcal_nm1_E_over_L_notmu_etight -> Fill( nm1_E / nm1_L ) ;
            }
         }

     //--- Loose isolation cut

         if ( nm1_chi2_iso > iso_loose ) continue ;

         h_hcal_nm1_E_vs_L_all_isoloose -> Fill( nm1_L, nm1_E ) ;
         h_energy_hcal_all_isoloose -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
         h_energy_ecal_all_isoloose -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
         h_hcal_nm1_E_over_L_all_isoloose -> Fill( nm1_E / nm1_L ) ;
         if ( ismu ) {
            h_hcal_nm1_E_vs_L_mu_isoloose -> Fill( nm1_L, nm1_E ) ;
            h_energy_hcal_mu_isoloose -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
            h_energy_ecal_mu_isoloose -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
            h_hcal_nm1_E_over_L_mu_isoloose -> Fill( nm1_E / nm1_L ) ;
         } else {
            h_hcal_nm1_E_vs_L_notmu_isoloose -> Fill( nm1_L, nm1_E ) ;
            h_energy_hcal_notmu_isoloose -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
            h_energy_ecal_notmu_isoloose -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
            h_hcal_nm1_E_over_L_notmu_isoloose -> Fill( nm1_E / nm1_L ) ;
         }

        //--- add loose energy chi2 cut
         if ( nm1_chi2_energy < energy_loose ) {
            h_hcal_nm1_E_vs_L_all_isoloose_eloose -> Fill( nm1_L, nm1_E ) ;
            h_hcal_nm1_E_over_L_all_isoloose_eloose -> Fill( nm1_E / nm1_L ) ;
            if ( ismu ) {
               h_hcal_nm1_E_vs_L_mu_isoloose_eloose -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_mu_isoloose_eloose -> Fill( nm1_E / nm1_L ) ;
            } else {
               h_hcal_nm1_E_vs_L_notmu_isoloose_eloose -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_notmu_isoloose_eloose -> Fill( nm1_E / nm1_L ) ;
            }

         }

        //--- add tight energy chi2 cut
         if ( nm1_chi2_energy < energy_tight ) {
            h_hcal_nm1_E_vs_L_all_isoloose_etight -> Fill( nm1_L, nm1_E ) ;
            h_hcal_nm1_E_over_L_all_isoloose_etight -> Fill( nm1_E / nm1_L ) ;
            if ( ismu ) {
               h_hcal_nm1_E_vs_L_mu_isoloose_etight -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_mu_isoloose_etight -> Fill( nm1_E / nm1_L ) ;
            } else {
               h_hcal_nm1_E_vs_L_notmu_isoloose_etight -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_notmu_isoloose_etight -> Fill( nm1_E / nm1_L ) ;
            }

         }


     //--- tight isolation cut

         if ( nm1_chi2_iso > iso_tight ) continue ;

         h_hcal_nm1_E_vs_L_all_isotight -> Fill( nm1_L, nm1_E ) ;
         h_energy_hcal_all_isotight -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
         h_energy_ecal_all_isotight -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
         h_hcal_nm1_E_over_L_all_isotight -> Fill( nm1_E / nm1_L ) ;
         if ( ismu ) {
            h_hcal_nm1_E_vs_L_mu_isotight -> Fill( nm1_L, nm1_E ) ;
            h_energy_hcal_mu_isotight -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
            h_energy_ecal_mu_isotight -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
            h_hcal_nm1_E_over_L_mu_isotight -> Fill( nm1_E / nm1_L ) ;
         } else {
            h_hcal_nm1_E_vs_L_notmu_isotight -> Fill( nm1_L, nm1_E ) ;
            h_energy_hcal_notmu_isotight -> Fill( (*rc_nm1_energy_hcal)[rci] ) ;
            h_energy_ecal_notmu_isotight -> Fill( (*rc_nm1_energy_ecal)[rci] ) ;
            h_hcal_nm1_E_over_L_notmu_isotight -> Fill( nm1_E / nm1_L ) ;
         }

        //--- add loose energy chi2 cut
         if ( nm1_chi2_energy < energy_loose ) {
            h_hcal_nm1_E_vs_L_all_isotight_eloose -> Fill( nm1_L, nm1_E ) ;
            h_hcal_nm1_E_over_L_all_isotight_eloose -> Fill( nm1_E / nm1_L ) ;
            if ( ismu ) {
               h_hcal_nm1_E_vs_L_mu_isotight_eloose -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_mu_isotight_eloose -> Fill( nm1_E / nm1_L ) ;
            } else {
               h_hcal_nm1_E_vs_L_notmu_isotight_eloose -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_notmu_isotight_eloose -> Fill( nm1_E / nm1_L ) ;
            }

         }

        //--- add tight energy chi2 cut
         if ( nm1_chi2_energy < energy_tight ) {
            h_hcal_nm1_E_vs_L_all_isotight_etight -> Fill( nm1_L, nm1_E ) ;
            h_hcal_nm1_E_over_L_all_isotight_etight -> Fill( nm1_E / nm1_L ) ;
            if ( ismu ) {
               h_hcal_nm1_E_vs_L_mu_isotight_etight -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_mu_isotight_etight -> Fill( nm1_E / nm1_L ) ;
            } else {
               h_hcal_nm1_E_vs_L_notmu_isotight_etight -> Fill( nm1_L, nm1_E ) ;
               h_hcal_nm1_E_over_L_notmu_isotight_etight -> Fill( nm1_E / nm1_L ) ;
            }

         }


      } // rci



   } // jentry

   saveHist( "plots.root", "h*" ) ;

   draw_plots1( "plots.root" ) ;

} // Loop



