#define fill_hists1_cxx
#include "fill_hists1.h"
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
void fill_hists1::Loop( bool verb, int first_event )
{
   if (fChain == 0) return;

   gDirectory -> Delete( "h*" ) ;


   Long64_t nentries = fChain->GetEntriesFast();

   int bins = 60 ;

   TH1F* h_deta = new TH1F( "h_deta", "deta, track momentum vs proj. momentum", bins, -1., 1. ) ;
   TH1F* h_dphi = new TH1F( "h_dphi", "dphi, track momentum vs proj. momentum", bins, -1., 1. ) ;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=first_event; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;



      if ( verb ) printf("\n\n ====== Event %5llu\n\n", jentry ) ;
      if ( verb ) printf("  Reconstructed tracks:  %3d\n", rcN ) ;

      for ( int ti=0; ti<rcEta->size(); ti++ ) {
         if ( verb ) printf(" %3d :  Eta = %7.3f, Phi = %7.3f, Pt = %7.3f,  nFST = %d, nFTT=%d\n",
             ti, rcEta->at(ti), rcPhi->at(ti), rcPt->at(ti), rcNumFST->at(ti), rcNumFTT->at(ti) ) ;
         if ( rcNumFST->at(ti) < 4 ) continue ;
         if ( rcEta->at(ti) < 2.5 ) continue ;
            if ( verb ) printf("   ECAL proj:  (%8.3f, %8.3f, %8.3f)\n", rcProjEcalx->at(ti), rcProjEcaly->at(ti), rcProjEcalz->at(ti) ) ;
            if ( verb ) printf("   HCAL proj:  (%8.3f, %8.3f, %8.3f)\n", rcProjHcalx->at(ti), rcProjHcaly->at(ti), rcProjHcalz->at(ti) ) ;
            if ( verb ) printf("   proj px,py,pz:  (%8.3f, %8.3f, %8.3f)\n", rcProjEcalPx->at(ti), rcProjEcalPy->at(ti), rcProjEcalPz->at(ti) ) ;
            float tpx, tpy, tpz ;
            eta_phi_pt_to_px_py_pz( rcEta->at(ti), rcPhi->at(ti), rcPt->at(ti), tpx, tpy, tpz ) ;
            float peta, pphi, ppt ;
            px_py_pz_to_eta_phi_pt( rcProjEcalPx->at(ti), rcProjEcalPy->at(ti), rcProjEcalPz->at(ti), peta, pphi, ppt ) ;
            float deta = (peta-rcEta->at(ti)) ;
            float dphi = (pphi - rcPhi->at(ti)) ;
            if ( dphi > 3.14159265 ) dphi = dphi - 2*3.14159265 ;
            if ( dphi <-3.14159265 ) dphi = dphi + 2*3.14159265 ;
            if ( verb ) printf("   trk  px,py,pz:  (%8.3f, %8.3f, %8.3f)\n", tpx, tpy, tpz ) ;
            if ( verb ) printf("  proj  eta,phi,pt:  (%7.3f, %7.3f, %7.3f)    deta = %7.3f,  dphi = %7.3f\n", peta, pphi, ppt, deta, dphi ) ;

         h_deta -> Fill( deta ) ;
         h_dphi -> Fill( dphi ) ;

      } // ti


      if ( verb ) {
         char answ = getchar() ;
         if ( answ == 'q' ) return ;
      }


   } // jentry

   gDirectory -> ls() ;

}

