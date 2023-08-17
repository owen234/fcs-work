#define display_clusters1_cxx
#include "display_clusters1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

float cos_angle = 0.999544 ;
float sin_angle = 0.030190 ;

//==========================
//float hcal_fcs_zmin = 782.63 - 710.16 ; // = 72.47
//float hcal_fcs_zmax = hcal_fcs_zmin + 84.24 ; // or should this be 117.0?
//float ecal_fcs_zmin = 13.9 ; // front space
//float ecal_fcs_zmax = 71.21 - 18.37 ; // total - back space = 52.84
//==========================

float hcal_fcs_zmin = 782.63 - 710.16 ; // = 72.47
float hcal_fcs_zmax = hcal_fcs_zmin + 84.24 ; // or should this be 117.0?
float ecal_fcs_zmin = 13.9 ; // front space
float ecal_fcs_zmax = 71.21 - 18.37 ; // total - back space = 52.84

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
int ecal_row_col_to_id( int row, int col ) {
   ///return (row-1)*22 + col ;
   return (row-1)*22 + col -1 ;
}
//---------------
int hcal_row_col_to_id( int row, int col ) {
   ///return (row-1)*13 + col ;
   return (row-1)*13 + col -1 ;
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
void ecal_fcsxy_to_row_col( float fcsx, float fcsy,
       int& row, int& col ) {
   col = TMath::Nint( fcsx/5.572 + 0.5 ) ;
   row = TMath::Nint((-1 * (fcsy - 0.31 + 5.572) + 34.*5.572/2.) / 5.572 + 0.5) ;
}
//---------------
void hcal_fcsxy_to_row_col( float fcsx, float fcsy,
       int& row, int& col ) {
   col = TMath::Nint( fcsx/9.99 + 0.5 ) ;
   row = TMath::Nint((-1 * (fcsy - 1.70) + 20.*9.99/2.) / 9.99 + 0.5) ;
}
//---------------
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

void display_clusters1::Loop()
{

   gSystem -> Exec( "mkdir -p plots" ) ;

   gStyle -> SetPadRightMargin(0.05) ;
   gStyle -> SetPadLeftMargin(0.15) ;
   gStyle -> SetPadTopMargin(0.05) ;
   gStyle -> SetPadBottomMargin(0.15) ;
   gStyle -> SetOptTitle(0) ;
   if (fChain == 0) return;


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
   //tpm_tp -> SetMarkerStyle(24) ;
   tpm_tp -> SetMarkerStyle(20) ;
   tpm_tp -> SetMarkerColor(1) ;
   tpm_tp -> SetMarkerSize(1.0) ;

   TPolyMarker* tpm_cl = new TPolyMarker() ;
   tpm_cl -> SetMarkerStyle(22) ;
   tpm_cl -> SetMarkerColor(2) ;
   tpm_cl -> SetMarkerSize(1.5) ;


   TBox* tb_ecal = new TBox() ;
   tb_ecal -> SetFillColor(4) ;
   tb_ecal -> SetFillStyle(1001) ;
   tb_ecal -> SetLineColor(4) ;
   TBox* tb_ecal_line = new TBox() ;
   tb_ecal_line -> SetFillStyle(0) ;
   tb_ecal_line -> SetLineColor(4) ;

   TBox* tb_ecal_line2 = new TBox() ;
   tb_ecal_line2 -> SetFillStyle(0) ;
   tb_ecal_line2 -> SetLineColor(17) ;

   TBox* tb_hcal = new TBox() ;
   tb_hcal -> SetFillColor(2) ;
   tb_hcal -> SetFillStyle(1001) ;
   tb_hcal -> SetLineColor(2) ;
   TBox* tb_hcal_line = new TBox() ;
   tb_hcal_line -> SetFillStyle(0) ;
   tb_hcal_line -> SetLineColor(2) ;

   TEllipse* te_circle = new TEllipse() ;
   te_circle -> SetLineWidth(2) ;
   te_circle -> SetLineColor(2) ;
   te_circle -> SetFillStyle(0) ;


   TCanvas* can_xy(0x0) ;
   TCanvas* can_xy2(0x0) ;
      gStyle -> SetOptStat(0) ;
      can_xy = new TCanvas( "can_xy", "xy", 50, 50, 800, 800 ) ;
      can_xy -> Draw() ;
      can_xy2 = new TCanvas( "can_xy2", "xy2", 850, 50, 800, 800 ) ;
      can_xy2 -> Draw() ;
      gSystem -> ProcessEvents() ;

   TH2F* h_display_xy = new TH2F( "h_display_xy", "Display, xy", 800, -200., 200., 800, -200., 200. ) ;



   Long64_t nentries = fChain->GetEntries();
   printf("\n\n Chain has %llu entries\n\n", nentries ) ;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      printf("\n\n\n ============= Event %llu\n\n", jentry ) ;


         can_xy -> cd() ;
         h_display_xy -> SetXTitle( "x position (cm)" ) ;
         h_display_xy -> SetYTitle( "y position (cm)" ) ;
         h_display_xy -> SetTitleOffset( 1.4, "y") ;
         h_display_xy -> SetTitleOffset( 1.2, "x") ;
         h_display_xy -> Draw() ;

         double x[1] ;
         double y[1] ;

         for ( int i=0; i<fcs_rec_ecalE->size(); i++ ) {

            float energy = fcs_rec_ecalE->at(i) ;
            if ( energy < 0.0001 ) continue ;

            float x, y ;
            x = fcs_rec_ecalX->at(i) ;
            y = fcs_rec_ecalY->at(i) ;

            int id = fcs_rec_ecalId -> at(i) ;
            int row, col ;
            ecal_id_to_row_col( id, row, col ) ;
            int ns = fcs_rec_ecalDet->at(i) % 2 ;



            float alpha = energy  ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_ecal -> SetFillColorAlpha( kBlue, alpha ) ;
            tb_ecal -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
            tb_ecal_line -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
            printf( "  ECAL hit %3d,  id=%4d, row = %2d, col = %2d, ns = %d,  x = %7.2f  , y = %7.2f  ,  E = %7.3f\n",
                i, id, row, col, ns, x, y, energy ) ;
         }


         for ( int ti=0; ti<rcN; ti++ ) {
            float tpx = rcProjEcalx->at(ti) ;
            float tpy = rcProjEcaly->at(ti) ;
            x[0] = tpx ;
            y[0] = tpy ;
            tpm_tp -> DrawPolyMarker( 1, x, y ) ;
            float radius = 20. ;
            te_circle -> DrawEllipse( x[0], y[0], radius, radius, 0., 360., 0. ) ;
            printf(" --- Track %2d :  x = %7.2f   y = %7.2f\n", ti, tpx, tpy ) ;
            for ( int hi=0; hi<fcs_rec_ecalE->size(); hi++ ) {
               float energy = fcs_rec_ecalE->at(hi) ;
               if ( energy < 0.0001 ) continue ;
               float hx = fcs_rec_ecalX->at(hi) ;
               float hy = fcs_rec_ecalY->at(hi) ;
               float dx = hx - tpx ;
               float dy = hy - tpy ;
               float dr = sqrt( dx*dx + dy*dy ) ;
               if ( dr < radius ) {
                  printf( " ECAL hit %3d :  E = %7.3f  x = %7.2f   y = %7.2f     dx = %7.2f  dy = %7.2f  dr = %7.2f\n",
                     hi, energy, hx, hy, dx, dy, dr ) ;
               }
            } // hi
            for ( int hi=0; hi<fcs_rec_hcalE->size(); hi++ ) {
               float energy = fcs_rec_hcalE->at(hi) ;
               if ( energy < 0.0001 ) continue ;
               float hx = fcs_rec_hcalX->at(hi) ;
               float hy = fcs_rec_hcalY->at(hi) ;
               float dx = hx - tpx ;
               float dy = hy - tpy ;
               float dr = sqrt( dx*dx + dy*dy ) ;
               if ( dr < radius ) {
                  printf( " HCAL hit %3d :  E = %7.3f  x = %7.2f   y = %7.2f     dx = %7.2f  dy = %7.2f  dr = %7.2f\n",
                     hi, energy, hx, hy, dx, dy, dr ) ;
               }
            } // hi
         } // ti

         for ( int ci = 0; ci < fcs_cl_ecalN; ci++ ) {
            x[0] = fcs_cl_ecalX->at(ci) ;
            y[0] = fcs_cl_ecalY->at(ci) ;
            /////////tpm_cl -> DrawPolyMarker( 1, x, y ) ;
            printf("  ECAL cluster %2d :  x = %7.2f, y = %7.2f\n", ci, x[0], y[0] ) ;
         }
     /// for ( int ci = 0; ci < fcs_cl_hcalN; ci++ ) {
     ///    x[0] = fcs_cl_hcalX->at(ci) ;
     ///    y[0] = fcs_cl_hcalY->at(ci) ;
     ///    tpm_cl -> DrawPolyMarker( 1, x, y ) ;
     ///    printf("  HCAL cluster %2d :  x = %7.2f, y = %7.2f\n", ci, x[0], y[0] ) ;
     /// }

         can_xy -> Update() ;
         can_xy -> Draw() ;
         gSystem -> ProcessEvents() ;


  //------------------------------------------------------

            can_xy2 -> cd() ;
            h_display_xy -> Draw() ;

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

            float sx = fcs_rec_hcalX->at(i) ;
            float sy = fcs_rec_hcalY->at(i) ;
            x = starx ;
            y = stary ;
            ////////float alpha = energy / 0.05 ;
            float alpha = energy / 2. ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_hcal -> SetFillColorAlpha( kRed, alpha ) ;
            tb_hcal -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
            tb_hcal_line -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
        /// printf( "  HCAL hit %3d,  id=%4d, row = %2d, col = %2d, ns = %d,  x = %7.2f (%7.2f), y = %7.2f (%7.2f),  E = %7.3f\n",
        ///     i, id, row, col, ns, x, sx, y, sy, energy ) ;
         }

         for ( int ti=0; ti<rcN; ti++ ) {
            float tpx = rcProjEcalx->at(ti) ;
            float tpy = rcProjEcaly->at(ti) ;
            x[0] = tpx ;
            y[0] = tpy ;
            tpm_tp -> DrawPolyMarker( 1, x, y ) ;
            float radius = 20. ;
            te_circle -> DrawEllipse( x[0], y[0], radius, radius, 0., 360., 0. ) ;
            printf(" --- Track %2d :  x = %7.2f   y = %7.2f\n", ti, tpx, tpy ) ;
            for ( int hi=0; hi<fcs_rec_ecalE->size(); hi++ ) {
               float energy = fcs_rec_ecalE->at(hi) ;
               if ( energy < 0.0001 ) continue ;
               float hx = fcs_rec_ecalX->at(hi) ;
               float hy = fcs_rec_ecalY->at(hi) ;
               float dx = hx - tpx ;
               float dy = hy - tpy ;
               float dr = sqrt( dx*dx + dy*dy ) ;
               if ( dr < radius ) {
                  printf( " ECAL hit %3d :  E = %7.3f  x = %7.2f   y = %7.2f     dx = %7.2f  dy = %7.2f  dr = %7.2f\n",
                     hi, energy, hx, hy, dx, dy, dr ) ;
               }
            } // hi
            for ( int hi=0; hi<fcs_rec_hcalE->size(); hi++ ) {
               float energy = fcs_rec_hcalE->at(hi) ;
               if ( energy < 0.0001 ) continue ;
               float hx = fcs_rec_hcalX->at(hi) ;
               float hy = fcs_rec_hcalY->at(hi) ;
               float dx = hx - tpx ;
               float dy = hy - tpy ;
               float dr = sqrt( dx*dx + dy*dy ) ;
               if ( dr < radius ) {
       ////       printf( " HCAL hit %3d :  E = %7.3f  x = %7.2f   y = %7.2f     dx = %7.2f  dy = %7.2f  dr = %7.2f\n",
       ////          hi, energy, hx, hy, dx, dy, dr ) ;
               }
            } // hi
         } // ti


 /////   for ( int ci = 0; ci < fcs_cl_ecalN; ci++ ) {

 /////      printf("  ECAL cluster %d\n", ci ) ;
 /////      can_xy2 -> cd() ;
 /////      h_display_xy -> Draw() ;

 /////      for ( int hi = 0; hi < fcs_rec_ecalN; hi++ ) {


 /////         int id = fcs_rec_ecalId -> at(hi) ;
 /////         int row, col ;
 /////         ecal_id_to_row_col( id, row, col ) ;
 /////         float fcsx, fcsy ;
 /////         ecal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
 /////         float starx(0.) ;
 /////         float stary(0.) ;
 /////         float starz(0.) ;
 /////         int ns = fcs_rec_ecalDet->at(hi) % 2 ;
 /////         fcs_to_star_ecal( fcsx, fcsy, 33.3661, ns,   starx, stary, starz ) ;

 /////         float x, y ;
 /////         x = starx ;
 /////         y = stary ;

 /////         tb_ecal_line2 -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;

 /////         if ( fcs_rec_ecalClIndex->at(hi) != ci ) continue ;

 /////         tb_ecal_line -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;

 /////         float energy = fcs_rec_ecalE->at(hi)  ;

 /////         float alpha = energy ;
 /////         if ( alpha > 1 ) alpha = 1 ;
 /////         if ( alpha < 0 ) alpha = 0 ;
 /////         tb_ecal -> SetFillColorAlpha( kBlue, alpha ) ;
 /////         tb_ecal -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
 /////         printf( "  ECAL hit %3d,  id=%4d, row = %2d, col = %2d, ns = %d,  x = %7.2f, y = %7.2f,  E = %7.3f\n",
 /////             hi, id, row, col, ns, x, y, energy ) ;


 /////      } // hi


 /////      char answ = getchar() ;
 /////      if ( answ == 'q' ) return ;

 /////   } // ci

            can_xy2 -> Update() ;
            can_xy2 -> Draw() ;
            gSystem -> ProcessEvents() ;

            can_xy -> SaveAs("plots/display_xy.pdf" ) ;
            can_xy -> SaveAs("plots/display_xy.png" ) ;


            printf("\n\n\n ============= Event %llu\n\n", jentry ) ;

            char answ = getchar() ;
            if ( answ == 'q' ) return ;



   }

} // Loop




