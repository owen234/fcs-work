#define display_data_cxx
#include "display_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

float cos_angle = 0.999544 ;
float sin_angle = 0.030190 ;

float hcal_fcs_zmin = 782.63 - 710.16 ; // = 72.47
float hcal_fcs_zmax = hcal_fcs_zmin + 84.24 ; // or should this be 117.0?
float ecal_fcs_zmin = 13.9 ; // front space
float ecal_fcs_zmax = 71.21 - 18.37 ; // total - back space = 52.84

//---------------
void simple_projection( float px, float py, float pz, float z, float& projx, float& projy ) {
   projx = 0. ;
   projy = 0. ;
   if ( pz <= 0 ) return ;
   projx = (px/pz)*z ;
   projy = (py/pz)*z ;
}
//---------------
void eta_phi_pt_to_px_py_pz(  float eta, float phi, float pt, float& px, float& py, float& pz ) {
   px = 0. ;
   py = 0. ;
   pz = 0. ;
   px = pt*cos(phi) ;
   py = pt*sin(phi) ;
   float theta = 2 * atan( exp( -1*eta ) ) ;
   pz = pt / tan(theta) ;
}
//---------------
void px_py_pz_to_eta_phi_pt( float px, float py, float pz, float& eta, float& phi, float& pt ) {
   eta = 0. ;
   phi = 0. ;
   pt = 0. ;
   pt = sqrt( px*px + py*py ) ;
   phi = atan2( py, px ) ;
   float theta = atan2( pt, pz ) ;
   eta = -1 * log( tan(theta/2) ) ;
}
//---------------
void ecal_id_to_row_col( int id, int& row, int& col ) {
   row = id/22+1 ;
   col = id%22+1 ;
}
//---------------
void hcal_id_to_row_col( int id, int& row, int& col ) {
   row = id/13+1 ;
   col = id%13+1 ;
}
//---------------
void ecal_row_col_to_fcsxy( int row, int col,
      float& fcsx, float& fcsy ) {
    fcsx = (col-0.5)*5.572 ;
    fcsy = -1* ((row-0.5)*5.572 - 34.*5.572/2.) +0.31 - 5.572 ;
}
//---------------
void hcal_row_col_to_fcsxy( int row, int col,
      float& fcsx, float& fcsy ) {
    fcsx = (col-0.5)*9.99 + 1.54 ;
    fcsy = -1* ((row-0.5)*9.99 - 20.*9.99/2.) +1.70 ;
}
//---------------
void fcs_to_star_ecal( float fcsx, float fcsy, float fcsz, int ns,
                       float& starx, float& stary, float& starz ) {
    ////////////printf(" +++ fcs_to_star_ecal : input fcs  coords (%7.3f, %7.3f, %7.3f)\n", fcsx, fcsy, fcsz ) ;
    if ( ns == 1 ) {
       float tmpx = fcsx ;
       float tmpz = fcsz ;
       stary = fcsy ;
       starx =      cos_angle * tmpx + 1 * sin_angle * tmpz ;
       starz =      cos_angle * tmpz - 1 * sin_angle * tmpx ;
       ////////////printf(" +++ fcs_to_star_ecal : after  rotation, coords (%7.3f, %7.3f, %7.3f)\n", starx, stary, starz  ) ;
       starx = starx + 17.40 ;
       starz = starz + 710.16 ;
    } else {
       float tmpx = -1*fcsx ;
       float tmpz = fcsz ;
       stary = fcsy ;
       starx =      cos_angle * tmpx - 1 * sin_angle * tmpz ;
       starz =      cos_angle * tmpz + 1 * sin_angle * tmpx ;
       ////////////printf(" +++ fcs_to_star_ecal : after  rotation, coords (%7.3f, %7.3f, %7.3f)\n", starx, stary, starz  ) ;
       starx = starx - 17.40 ;
       starz = starz + 710.16 ;
    }

}
//---------------
void star_to_fcs_ecal( float starx, float stary, float starz,
                       float& fcsx, float& fcsy, float& fcsz ) {
    ////////////printf(" --- star_to_fcs_ecal : input star coords (%7.3f, %7.3f, %7.3f)\n", starx, stary, starz ) ;
    if ( starx > 0 ) {
       fcsx = starx - 17.40 ;
       fcsy = stary ;
       fcsz = starz - 710.16 ;
       float tmpx = fcsx ;
       float tmpz = fcsz ;
       ////////////printf(" --- star_to_fcs_ecal : before rotation, fcs coords (%7.3f, %7.3f, %7.3f)\n", fcsx, fcsy, fcsz ) ;
       fcsx =      cos_angle * tmpx - 1 * sin_angle * tmpz ;
       fcsz =      cos_angle * tmpz + 1 * sin_angle * tmpx ;
       ////////////printf(" --- star_to_fcs_ecal : after  rotation, fcs coords (%7.3f, %7.3f, %7.3f)\n", fcsx, fcsy, fcsz ) ;
    } else {
       fcsx = starx + 17.40 ;
       fcsy = stary  ;
       fcsz = starz - 710.16 ;
       float tmpx = fcsx ;
       float tmpz = fcsz ;
       ////////////printf(" --- star_to_fcs_ecal : before rotation, fcs coords (%7.3f, %7.3f, %7.3f)\n", fcsx, fcsy, fcsz ) ;
       fcsx =      cos_angle * tmpx + 1 * sin_angle * tmpz ;
       fcsz =      cos_angle * tmpz - 1 * sin_angle * tmpx ;
       ////////////printf(" --- star_to_fcs_ecal : after  rotation, fcs coords (%7.3f, %7.3f, %7.3f)\n", fcsx, fcsy, fcsz ) ;
       fcsx = -1 * fcsx ;
    }
}
//---------------
void display_data::Loop( int first_event )
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


         can_xy -> cd() ;
         h_display_xy -> SetXTitle( "x position (cm)" ) ;
         h_display_xy -> SetYTitle( "y position (cm)" ) ;
         h_display_xy -> Draw() ;

         double x[1] ;
         double y[1] ;

  /////  for ( int i=0; i<fstX->size(); i++ ) {
  /////     x[0] = fstX->at(i) ;
  /////     y[0] = fstY->at(i) ;
  /////     tpm_fst -> DrawPolyMarker( 1, x, y ) ;
  /////  }
  /////  for ( int i=0; i<fttX->size(); i++ ) {
  /////     x[0] = fttX->at(i) ;
  /////     y[0] = fttY->at(i) ;
  /////     tpm_ftt -> DrawPolyMarker( 1, x, y ) ;
  /////  }


         for ( int i=0; i<fcs_rec_ecalE->size(); i++ ) {

            float energy = fcs_rec_ecalE->at(i) ;
            if ( energy < 0.0001 ) continue ;

            float x, y ;
            ///////   x = fcs_rec_ecalX->at(i) ;
            ///////   y = fcs_rec_ecalY->at(i) ;
            int id = fcs_rec_ecalId -> at(i) ;
            int row, col ;
            ecal_id_to_row_col( id, row, col ) ;
            float fcsx, fcsy ;
            ecal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
            float starx(0.) ;
            float stary(0.) ;
            float starz(0.) ;
            int ns = fcs_rec_ecalDet->at(i) % 2 ;
            fcs_to_star_ecal( fcsx, fcsy, 33.3661, ns,   starx, stary, starz ) ;

            x = starx ;
            y = stary ;



            /////float alpha = energy / 0.08 ;
            float alpha = energy  ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_ecal -> SetFillColorAlpha( kBlue, alpha ) ;
            tb_ecal -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
            tb_ecal_line -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
            printf( "  ECAL hit %3d,  id=%4d, row = %2d, col = %2d, ns = %d,  x = %7.2f, y = %7.2f,  E = %7.3f\n",
                i, id, row, col, ns, x, y, energy ) ;
         }

         for ( int i=0; i<fcs_rec_hcalE->size(); i++ ) {

            float energy = fcs_rec_hcalE->at(i) ;
            if ( energy < 0.0001 ) continue ;

            float x, y ;
               //////////x = fcs_rec_hcalX->at(i) ;
               //////////y = fcs_rec_hcalY->at(i) ;

            int id = fcs_rec_hcalId -> at(i) ;
            int row, col ;
            hcal_id_to_row_col( id, row, col ) ;
            float fcsx, fcsy ;
            hcal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
            float starx(0.) ;
            float stary(0.) ;
            float starz(0.) ;
            int ns = fcs_rec_hcalDet->at(i) % 2 ;
            fcs_to_star_ecal( fcsx, fcsy, 114.6741, ns,   starx, stary, starz ) ;

            x = starx ;
            y = stary ;
            ////////float alpha = energy / 0.05 ;
            float alpha = energy / 2. ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_hcal -> SetFillColorAlpha( kRed, alpha ) ;
            tb_hcal -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
            tb_hcal_line -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
            printf( "  HCAL hit %3d,  id=%4d, row = %2d, col = %2d, ns = %d,  x = %7.2f, y = %7.2f,  E = %7.3f\n",
                i, id, row, col, ns, x, y, energy ) ;
         }


         can_xy -> Update() ;
         can_xy -> Draw() ;
         gSystem -> ProcessEvents() ;




         can_rz -> cd() ;
         h_display_rz -> SetXTitle("z position (cm)") ;
         h_display_rz -> SetYTitle("rho position (cm)") ;
         h_display_rz -> Draw() ;

         double r[1] ;
         double z[1] ;

  ////   for ( int i=0; i<fstX->size(); i++ ) {
  ////      r[0] = sqrt( fstX->at(i) * fstX->at(i) + fstY->at(i) * fstY->at(i) ) ;
  ////      z[0] = fstZ->at(i) ;
  ////      tpm_fst -> DrawPolyMarker( 1, z, r ) ;
  ////   }
  ////   for ( int i=0; i<fttX->size(); i++ ) {
  ////      r[0] = sqrt( fttX->at(i) * fttX->at(i) + fttY->at(i) * fttY->at(i) ) ;
  ////      z[0] = fttZ->at(i) ;
  ////      tpm_ftt -> DrawPolyMarker( 1, z, r ) ;
  ////   }




         for ( int i=0; i<fcs_rec_ecalX->size(); i++ ) {
            float r1, z1, r2, z2 ;
            float x, y, r ;
            /////////x = fcs_rec_ecalX->at(i)  ;
            /////////y = fcs_rec_ecalY->at(i)  ;
            int id = fcs_rec_ecalId -> at(i) ;
            int row, col ;
            ecal_id_to_row_col( id, row, col ) ;
            float fcsx, fcsy ;
            ecal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
            float starx(0.) ;
            float stary(0.) ;
            float starz(0.) ;
            int ns = fcs_rec_ecalDet->at(i) % 2 ;
            fcs_to_star_ecal( fcsx, fcsy, 33.3661, ns,   starx, stary, starz ) ;

            x = starx ;
            y = stary ;

            r = sqrt( x*x + y*y ) ;
            r1 = r - 5.572/2. ;
            r2 = r + 5.572/2. ;
            z1 = 710.16 + ecal_fcs_zmin ;
            z2 = 710.16 + (71.21 - 18.37) ;
            float energy = fcs_rec_ecalE->at(i) ;
            /////////float alpha = energy / 0.08 ;
            float alpha = energy  ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_ecal -> SetFillColorAlpha( kBlue, alpha ) ;
            tb_ecal -> DrawBox( z1, r1, z2, r2 ) ;
            tb_ecal_line -> DrawBox( z1, r1, z2, r2 ) ;
         }

         printf("\n\n Number of HCAL hits:  %lu\n", fcs_rec_hcalE->size() ) ;
         for ( int i=0; i<fcs_rec_hcalE->size(); i++ ) {

            float energy = fcs_rec_hcalE->at(i) ;
            if ( energy < 0.0001 ) continue ;

            float r1, z1, r2, z2 ;
            float x, y, r ;
            ////////////x = fcs_mc_hcalX->at(i)  ;
            ////////////y = fcs_mc_hcalY->at(i)  ;
            int id = fcs_rec_hcalId -> at(i) ;
            int row, col ;
            hcal_id_to_row_col( id, row, col ) ;
            float fcsx, fcsy ;
            hcal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
            float starx(0.) ;
            float stary(0.) ;
            float starz(0.) ;
            int ns = fcs_rec_hcalDet->at(i) % 2 ;
            fcs_to_star_ecal( fcsx, fcsy, 114.6741, ns,   starx, stary, starz ) ;

            x = starx ;
            y = stary ;
            r = sqrt( x*x + y*y ) ;
            r1 = r - 9.99/2. ;
            r2 = r + 9.99/2. ;
            z1 = 782.63 ;
            z2 = z1 + 84.24 ;
            /////////float alpha = energy / 0.05 ;
            float alpha = energy / 2 ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_hcal -> SetFillColorAlpha( kRed, alpha ) ;
            printf("  HCAL:  %4d, id=%4d, row=%2d, col=%2d, x=%7.2f, y=%7.2f\n", i, id, row, col, x, y ) ;
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
         printf(" %3d :  Eta = %7.3f, Phi = %7.3f, Pt = %7.3f,  nFST = %d, nFTT=%d\n",
             ti, rcEta->at(ti), rcPhi->at(ti), rcPt->at(ti), rcNumFST->at(ti), rcNumFTT->at(ti) ) ;
         if ( rcNumFST->at(ti) < 4 ) continue ;
         if ( rcEta->at(ti) < 2.5 ) continue ;
            printf("   ECAL proj:  (%8.3f, %8.3f, %8.3f)\n", rcProjEcalx->at(ti), rcProjEcaly->at(ti), rcProjEcalz->at(ti) ) ;
            printf("   HCAL proj:  (%8.3f, %8.3f, %8.3f)\n", rcProjHcalx->at(ti), rcProjHcaly->at(ti), rcProjHcalz->at(ti) ) ;
            printf("   proj px,py,pz:  (%8.3f, %8.3f, %8.3f)\n", rcProjEcalPx->at(ti), rcProjEcalPy->at(ti), rcProjEcalPz->at(ti) ) ;
            float tpx, tpy, tpz ;
            eta_phi_pt_to_px_py_pz( rcEta->at(ti), rcPhi->at(ti), rcPt->at(ti), tpx, tpy, tpz ) ;
            float peta, pphi, ppt ;
            px_py_pz_to_eta_phi_pt( rcProjEcalPx->at(ti), rcProjEcalPy->at(ti), rcProjEcalPz->at(ti), peta, pphi, ppt ) ;
            printf("   trk  px,py,pz:  (%8.3f, %8.3f, %8.3f)\n", tpx, tpy, tpz ) ;
            float deta = (peta-rcEta->at(ti)) ;
            float dphi = (pphi - rcPhi->at(ti)) ;
            if ( dphi >  3.14159265 ) dphi = dphi - 2*3.14159265 ;
            if ( dphi < -3.14159265 ) dphi = dphi + 2*3.14159265 ;
            printf("  proj  eta,phi,pt:  (%7.3f, %7.3f, %7.3f)    deta = %7.3f,  dphi = %7.3f\n", peta, pphi, ppt, deta, dphi ) ;

            float simple_projx, simple_projy ;
            simple_projection( tpx, tpy, tpz, rcProjEcalz->at(ti), simple_projx, simple_projy ) ;
            printf(" simple proj  (%8.3f, %8.3f)\n", simple_projx, simple_projy ) ;

            x[0] = rcProjEcalx->at(ti) ;
            y[0] = rcProjEcaly->at(ti) ;
            can_xy -> cd() ;
            tl_pseudotrack -> DrawLine( 0., 0., x[0], y[0] ) ;
            can_rz -> cd() ;
            x[0] = rcProjEcalz->at(ti) ;
            y[0] = sqrt( rcProjEcalx->at(ti) * rcProjEcalx->at(ti)  +  rcProjEcaly->at(ti) * rcProjEcaly->at(ti) ) ;
            tl_pseudotrack -> DrawLine( 0., 0., x[0], y[0] ) ;
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

