#define display_only_cxx
#include "display_only.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

float cos_angle = 0.999544 ;
float sin_angle = 0.030190 ;

float hcal_fcs_zmin = 782.63 - 710.16 ; // = 72.47
float hcal_fcs_zmax = hcal_fcs_zmin + 84.24 ; // or should this be 117.0?
float ecal_fcs_zmin = 13.9 ; // front space
float ecal_fcs_zmax = 71.21 - 18.37 ; // total - back space = 52.84

void display_only::Loop( int first_event )
{
   if (fChain == 0) return;

   TPolyMarker* tpm_fst = new TPolyMarker() ;
   tpm_fst -> SetMarkerStyle(20) ;
   tpm_fst -> SetMarkerColor(2) ;
   tpm_fst -> SetMarkerSize(0.7) ;

   TPolyMarker* tpm_ftt = new TPolyMarker() ;
   tpm_ftt -> SetMarkerStyle(20) ;
   tpm_ftt -> SetMarkerColor(4) ;
   tpm_ftt -> SetMarkerSize(0.7) ;

   TPolyMarker* tpm_fcsc = new TPolyMarker() ;
   tpm_fcsc -> SetMarkerStyle(20) ;
   tpm_fcsc -> SetMarkerColor(3) ;
   tpm_fcsc -> SetMarkerSize(0.7) ;

   TPolyMarker* tpm_fcs_mc_ecal = new TPolyMarker() ;
   tpm_fcs_mc_ecal -> SetMarkerStyle(25) ;
   tpm_fcs_mc_ecal -> SetMarkerColor(4) ;
   tpm_fcs_mc_ecal -> SetMarkerSize(2.0) ;

   TPolyMarker* tpm_fcs_mc_hcal = new TPolyMarker() ;
   tpm_fcs_mc_hcal -> SetMarkerStyle(25) ;
   tpm_fcs_mc_hcal -> SetMarkerColor(2) ;
   tpm_fcs_mc_hcal -> SetMarkerSize(4.0) ;

   TPolyMarker* tpm_tp = new TPolyMarker() ;
   tpm_tp -> SetMarkerStyle(24) ;
   tpm_tp -> SetMarkerColor(1) ;
   tpm_tp -> SetMarkerSize(1.0) ;

   TLine* tl_pseudotrack = new TLine() ;
   TLine* tl_mctrack = new TLine() ;
   tl_mctrack -> SetLineColor(4) ;
   tl_mctrack -> SetLineStyle(2) ;

   TPolyMarker* tpm_min = new TPolyMarker() ;
   tpm_min -> SetMarkerStyle(20) ;
   tpm_min -> SetMarkerSize(0.5) ;

   TBox* tb_ecal = new TBox() ;
   tb_ecal -> SetFillColor(4) ;
   tb_ecal -> SetFillStyle(1001) ;
   tb_ecal -> SetLineColor(4) ;
   TBox* tb_ecal_line = new TBox() ;
   tb_ecal_line -> SetFillStyle(0) ;
   tb_ecal_line -> SetLineColor(4) ;

   TBox* tb_hcal = new TBox() ;
   tb_hcal -> SetFillColor(2) ;
   tb_hcal -> SetFillStyle(1001) ;
   tb_hcal -> SetLineColor(2) ;
   TBox* tb_hcal_line = new TBox() ;
   tb_hcal_line -> SetFillStyle(0) ;
   tb_hcal_line -> SetLineColor(2) ;


   TCanvas* can_xy(0x0) ;
   TCanvas* can_rz(0x0) ;
      gStyle -> SetOptStat(0) ;
      can_xy = new TCanvas( "can_xy", "xy", 50, 50, 800, 800 ) ;
      can_xy -> Draw() ;
      can_rz = new TCanvas( "can_rz", "rz", 850, 50, 900, 700 ) ;
      can_rz -> Draw() ;
      gSystem -> ProcessEvents() ;

   TH2F* h_display_xy = new TH2F( "h_display_xy", "Display, xy", 800, -200., 200., 800, -200., 200. ) ;
   TH2F* h_display_rz = new TH2F( "h_display_rz", "Display, rz", 1000, 0., 1000., 800, 0., 100. ) ;



   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=first_event; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

///   for ( int tpi=0; tpi<tprojN; tpi++ ) {
///      if ( tprojIdD -> at(tpi) != 4 ) continue ;
///      float pt = sqrt( tprojPx->at(tpi) * tprojPx->at(tpi) + tprojPy->at(tpi) * tprojPy->at(tpi) ) ;
///      float theta = atan2( pt, tprojPz->at(tpi) ) ;
///      float phi = atan2( tprojPy->at(tpi), tprojPx->at(tpi) ) ;
///      float eta = -1 * log( tan( theta/2. ) ) ;
///      printf(" tproj:  %3d  idt = %5d,  x,y,z = (%7.2f, %7.2f, %7.2f)  Eta = %7.3f, phi = %7.3f, pt = %7.2f\n",
///         tpi, tprojIdT->at(tpi), tprojX->at(tpi), tprojY->at(tpi), tprojZ->at(tpi), eta, phi, pt ) ;
///   } // tpi

  //------------

         can_xy -> cd() ;
         h_display_xy -> SetXTitle( "x position (cm)" ) ;
         h_display_xy -> SetYTitle( "y position (cm)" ) ;
         h_display_xy -> Draw() ;

         double x[1] ;
         double y[1] ;

         for ( int i=0; i<fstX->size(); i++ ) {
            x[0] = fstX->at(i) ;
            y[0] = fstY->at(i) ;
            tpm_fst -> DrawPolyMarker( 1, x, y ) ;
         }
         for ( int i=0; i<fttX->size(); i++ ) {
            x[0] = fttX->at(i) ;
            y[0] = fttY->at(i) ;
            tpm_ftt -> DrawPolyMarker( 1, x, y ) ;
         }
         for ( int i=0; i<fcs_mc_ecalX->size(); i++ ) {
            float x, y ;
               x = fcs_mc_ecalX->at(i) ;
               y = fcs_mc_ecalY->at(i) ;
            float energy = fcs_mc_ecalE->at(i) ;
            float alpha = energy / 0.08 ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_ecal -> SetFillColorAlpha( kBlue, alpha ) ;
            tb_ecal -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
            tb_ecal_line -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
         }
         for ( int i=0; i<fcs_mc_hcalX->size(); i++ ) {
            float x, y ;
               x = fcs_mc_hcalX->at(i) ;
               y = fcs_mc_hcalY->at(i) ;
            float energy = fcs_mc_hcalE->at(i) ;
            float alpha = energy / 0.05 ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_hcal -> SetFillColorAlpha( kRed, alpha ) ;
            tb_hcal -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
            tb_hcal_line -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
         }

 //      x[0] = tpx_ecal ;
 //      y[0] = tpy_ecal ;
 //      tpm_tp -> DrawPolyMarker( 1, x, y ) ;

 //      x[0] = tpx_hcal ;
 //      y[0] = tpy_hcal ;
 //      tpm_tp -> DrawPolyMarker( 1, x, y ) ;

 //      /////////////tl_pseudotrack -> DrawLine( 0., 0., x[0], y[0] ) ;

 //      for ( int i=0; i<mcEta->size(); i++ ) {
 //         printf(" N = %5d  size = %lu :  %4d : pid = %8d\n", mcN, mcEta->size(), i, mcPid->at(i) ) ;
 //         if ( abs(mcPid->at(i)) != 5 ) continue ;
 //         float mc_mu_eta = mcEta->at(i) ;
 //         float mc_mu_pt = mcPt->at(i) ;
 //         float mc_mu_phi = mcPhi->at(i) ;
 //         float mc_mu_px = mc_mu_pt * cos( mc_mu_phi ) ;
 //         float mc_mu_py = mc_mu_pt * sin( mc_mu_phi ) ;
 //         float mc_mu_theta = 2. * atan( exp( -1 * mc_mu_eta ) ) ;
 //         float mc_mu_pz = mc_mu_pt / tan( mc_mu_theta ) ;
 //         float mx_mu_x = mc_mu_px * (783. / mc_mu_pz ) ;
 //         float mx_mu_y = mc_mu_py * (783. / mc_mu_pz ) ;
 //         printf(" %d  MC muon:  pt = %5.2f, eta = %7.3f, phi = %7.3f,  x,y = (%7.2f, %7.2f)\n", i, mc_mu_pt, mc_mu_eta, mc_mu_phi, mx_mu_x, mx_mu_y ) ;
 //         printf(" N = %5d  size = %lu :  %4d : pid = %8d\n", mcN, mcEta->size(), i, mcPid->at(i) ) ;
 //         tl_mctrack -> DrawLine( 0., 0., mx_mu_x, mx_mu_y ) ;
 //      } // i

         can_xy -> Update() ;
         can_xy -> Draw() ;
         gSystem -> ProcessEvents() ;




         can_rz -> cd() ;
         h_display_rz -> SetXTitle("z position (cm)") ;
         h_display_rz -> SetYTitle("rho position (cm)") ;
         h_display_rz -> Draw() ;

         double r[1] ;
         double z[1] ;

         for ( int i=0; i<fstX->size(); i++ ) {
            r[0] = sqrt( fstX->at(i) * fstX->at(i) + fstY->at(i) * fstY->at(i) ) ;
            z[0] = fstZ->at(i) ;
            tpm_fst -> DrawPolyMarker( 1, z, r ) ;
         }
         for ( int i=0; i<fttX->size(); i++ ) {
            r[0] = sqrt( fttX->at(i) * fttX->at(i) + fttY->at(i) * fttY->at(i) ) ;
            z[0] = fttZ->at(i) ;
            tpm_ftt -> DrawPolyMarker( 1, z, r ) ;
         }
         for ( int i=0; i<fcs_mc_ecalX->size(); i++ ) {
            float r1, z1, r2, z2 ;
            float x, y, r ;
            x = fcs_mc_ecalX->at(i)  ;
            y = fcs_mc_ecalY->at(i)  ;
            r = sqrt( x*x + y*y ) ;
            r1 = r - 5.572/2. ;
            r2 = r + 5.572/2. ;
            z1 = 710.16 + ecal_fcs_zmin ;
            z2 = 710.16 + (71.21 - 18.37) ;
            float energy = fcs_mc_ecalE->at(i) ;
            float alpha = energy / 0.08 ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_ecal -> SetFillColorAlpha( kBlue, alpha ) ;
            tb_ecal -> DrawBox( z1, r1, z2, r2 ) ;
            tb_ecal_line -> DrawBox( z1, r1, z2, r2 ) ;
         }
         for ( int i=0; i<fcs_mc_hcalX->size(); i++ ) {
            float r1, z1, r2, z2 ;
            float x, y, r ;
            x = fcs_mc_hcalX->at(i)  ;
            y = fcs_mc_hcalY->at(i)  ;
            r = sqrt( x*x + y*y ) ;
            r1 = r - 9.99/2. ;
            r2 = r + 9.99/2. ;
            z1 = 782.63 ;
            z2 = z1 + 84.24 ;
            float energy = fcs_mc_hcalE->at(i) ;
            float alpha = energy / 0.05 ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_hcal -> SetFillColorAlpha( kRed, alpha ) ;
            tb_hcal -> DrawBox( z1, r1, z2, r2 ) ;
            tb_hcal_line -> DrawBox( z1, r1, z2, r2 ) ;
         }

//       r[0] = sqrt( tpx_ecal*tpx_ecal + tpy_ecal*tpy_ecal ) ;
//       z[0] = tpz_ecal ;
//       tpm_tp -> DrawPolyMarker( 1, z, r ) ;

//       r[0] = sqrt( tpx_hcal*tpx_hcal + tpy_hcal*tpy_hcal ) ;
//       z[0] = tpz_hcal ;
//       tpm_tp -> DrawPolyMarker( 1, z, r ) ;

         ///////tl_pseudotrack -> DrawLine( 0., 0., z[0], r[0] ) ;


         can_rz -> Update() ;
         can_rz -> Draw() ;
         gSystem -> ProcessEvents() ;


      printf("\n\n ====== Event %5llu\n\n", jentry ) ;
      printf("  Reconstructed tracks:  %3d\n", rcN ) ;

      for ( int ti=0; ti<rcEta->size(); ti++ ) {
         printf(" %3d :  Eta = %7.3f, Phi = %7.3f, Pt = %7.3f,  q = %2d,  Id = %3d, IdT = %3d\n",
             ti, rcEta->at(ti), rcPhi->at(ti), rcPt->at(ti), rcCharge->at(ti), rcTrackId->at(ti), rcTrackIdT->at(ti) ) ;
         float mindr(1e9) ;
         int mcind(-1) ;
         for ( int mci=0; mci<mcN; mci++ ) {
            if ( abs(mcCharge->at(mci)) != 1 ) continue ;
            float deta = mcEta->at(mci) - rcEta->at(ti) ;
            float dphi = mcPhi->at(mci) - rcPhi->at(ti) ;
            if ( dphi > 3.14159265 ) dphi = dphi - 2*3.14159265 ;
            if ( dphi <-3.14159265 ) dphi = dphi + 2*3.14159265 ;
            float dr = sqrt( deta*deta + dphi*dphi ) ;
            if ( dr < mindr ) {
               mindr = dr ;
               mcind = mci ;
            }
         }
         if ( mcind < 0 ) continue ;
         printf("   MC : Eta = %7.3f, Phi = %7.3f, Pt = %7.3f,  q = %2d,  Id = %3d, PID = %8d, dr = %7.3f\n",
                mcEta->at(mcind), mcPhi->at(mcind), mcPt->at(mcind), mcCharge->at(mcind), mcind, mcPid->at(mcind), mindr ) ;
         int tpi = rcTrackIdT->at(ti) ;
         for ( int tpi=0; tpi<tprojIdT->size(); tpi++ ) {
            if ( tprojIdT->at(tpi) != rcTrackIdT->at(ti) ) continue ;
            if ( tprojIdD -> at(tpi) != 4 ) continue ;
            float pt = sqrt( tprojPx->at(tpi) * tprojPx->at(tpi) + tprojPy->at(tpi) * tprojPy->at(tpi) ) ;
            float theta = atan2( pt, tprojPz->at(tpi) ) ;
            float phi = atan2( tprojPy->at(tpi), tprojPx->at(tpi) ) ;
            float eta = -1 * log( tan( theta/2. ) ) ;
            printf(" tproj: Eta = %7.3f, Phi = %7.3f, Pt = %7.3f  idt = %5d,  x,y,z = (%7.2f, %7.2f, %7.2f)\n",
               eta, phi, pt, tprojIdT->at(tpi), tprojX->at(tpi), tprojY->at(tpi), tprojZ->at(tpi) ) ;
            if ( abs(mcPid->at(mcind)) == 5 ) {
               tl_pseudotrack -> SetLineWidth(3) ;
            } else {
               tl_pseudotrack -> SetLineWidth(1) ;
            }
            x[0] = tprojX->at(tpi) ;
            y[0] = tprojY->at(tpi) ;
            can_xy -> cd() ;
            tl_pseudotrack -> DrawLine( 0., 0., x[0], y[0] ) ;
            can_rz -> cd() ;
            x[0] = tprojZ->at(tpi) ;
            y[0] = sqrt( tprojX->at(tpi) * tprojX->at(tpi) + tprojY->at(tpi)*tprojY->at(tpi) ) ;
            tl_pseudotrack -> DrawLine( 0., 0., x[0], y[0] ) ;
         } // tpi
         printf("\n") ;
      } // ti

         can_xy -> Update() ;
         can_xy -> Draw() ;
         gSystem -> ProcessEvents() ;
         can_rz -> Update() ;
         can_rz -> Draw() ;
         gSystem -> ProcessEvents() ;

            char answ = getchar() ;
            if ( answ == 'q' ) return ;

   } // jentry
}

