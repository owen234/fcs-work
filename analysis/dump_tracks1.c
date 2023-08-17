#define dump_tracks1_cxx
#include "dump_tracks1.h"
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

void dump_tracks1::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   printf("\n\n Chain has %llu entries\n\n", nentries ) ;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      for ( int ti=0; ti<rcN; ti++ ) {

         if ( rcNumFST->at(ti) <= 0 ) continue ;

         printf("\n\n\n ====== Event %7llu :  track %3d\n", jentry, ti ) ;
         printf("  Pt = %7.3f   Eta = %7.3f   Phi = %7.3f\n", rcPt->at(ti), rcEta->at(ti), rcPhi->at(ti) ) ;
         printf("  N FST = %d  N FTT = %d\n", rcNumFST->at(ti), rcNumFTT->at(ti) ) ;

         int ecal_ci = rcEcalClIndex->at(ti) ;
         int hcal_ci = rcHcalClIndex->at(ti) ;
         if ( ecal_ci >= 0 ) {
            printf("  ECAL cluster  %2d :  E = %7.2f  Dx = %7.2f  Dy = %7.2f  Dr = %7.2f\n",
                ecal_ci, rcEcalClE->at(ti), rcEcalClDx->at(ti), rcEcalClDy->at(ti) , rcEcalClDr->at(ti) ) ;
         }
         if ( hcal_ci >= 0 ) {
            printf("  HCAL cluster  %2d :  E = %7.2f  Dx = %7.2f  Dy = %7.2f  Dr = %7.2f\n",
                hcal_ci, rcHcalClE->at(ti), rcHcalClDx->at(ti), rcHcalClDy->at(ti), rcHcalClDr->at(ti) ) ;
         }
         printf("  ECAL projection:  (x,y,z) = (%7.2f, %7.2f, %7.2f)\n", rcProjEcalx->at(ti), rcProjEcaly->at(ti), rcProjEcalz->at(ti) ) ;
         if ( ecal_ci >= 0 ) {
            printf("  ECAL cluster   :  (x,y,z) = (%7.2f, %7.2f, %7.2f)\n", fcs_cl_ecalX->at(ecal_ci), fcs_cl_ecalY->at(ecal_ci), fcs_cl_ecalZ->at(ecal_ci) ) ;
         }
         if ( hcal_ci >= 0 ) {
            printf("  HCAL cluster   :  (x,y,z) = (%7.2f, %7.2f, %7.2f)\n", fcs_cl_hcalX->at(hcal_ci), fcs_cl_hcalY->at(hcal_ci), fcs_cl_hcalZ->at(hcal_ci) ) ;
         }

         if ( ecal_ci >= 0 ) {

            printf("  ECAL cluster hits:  %2d\n", fcs_cl_ecalNhit->at(ecal_ci) ) ;

            float sum_ex(0.) ;
            float sum_ey(0.) ;
            float sum_e(0.) ;

            for ( int hi=0; hi<fcs_rec_ecalN; hi++ ) {
               if ( fcs_rec_ecalClIndex->at(hi) == ecal_ci ) {

                  int id = fcs_rec_ecalId -> at(hi) ;
                  int row, col ;
                  ecal_id_to_row_col( id, row, col ) ;
                  float fcsx, fcsy ;
                  ecal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
                  float starx(0.) ;
                  float stary(0.) ;
                  float starz(0.) ;
                  int ns = fcs_rec_ecalDet->at(hi) % 2 ;
                  fcs_to_star_ecal( fcsx, fcsy, 33.3661, ns,   starx, stary, starz ) ;

                  float e = fcs_rec_ecalE->at(hi) ;
                  float x = starx ;
                  float y = stary ;

                  sum_ex += x*e ;
                  sum_ey += y*e ;
                  sum_e += e ;

                  printf("   ECAL hit  %3d :  E = %7.3f   (x,y,z) = (%7.2f, %7.2f, %7.2f)\n",
                   hi, e, x, y, fcs_rec_ecalZ->at(hi) ) ;
               }
            } // hi
            float calc_cl_x = sum_ex / sum_e ;
            float calc_cl_y = sum_ey / sum_e ;
            float dx = calc_cl_x - fcs_cl_ecalX->at(ecal_ci) ;
            float dy = calc_cl_y - fcs_cl_ecalY->at(ecal_ci) ;
            printf("    ECAL calculated cluster position:  (x,y) = (%7.2f, %7.2f)   E = %7.3f   diffs  dx = %7.3f  dy = %7.3f\n",
               calc_cl_x, calc_cl_y, sum_e, dx, dy ) ;

         }


         if ( hcal_ci >= 0 ) {

            printf("  HCAL cluster hits:  %2d\n", fcs_cl_hcalNhit->at(hcal_ci) ) ;

            float sum_ex(0.) ;
            float sum_ey(0.) ;
            float sum_e(0.) ;

            for ( int hi=0; hi<fcs_rec_hcalN; hi++ ) {
               if ( fcs_rec_hcalClIndex->at(hi) == hcal_ci ) {

                  int id = fcs_rec_hcalId -> at(hi) ;
                  int row, col ;
                  hcal_id_to_row_col( id, row, col ) ;
                  float fcsx, fcsy ;
                  hcal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
                  float starx(0.) ;
                  float stary(0.) ;
                  float starz(0.) ;
                  int ns = fcs_rec_hcalDet->at(hi) % 2 ;
                  fcs_to_star_ecal( fcsx, fcsy, 33.3661, ns,   starx, stary, starz ) ;

                  float e = fcs_rec_hcalE->at(hi) ;
                  float x = starx ;
                  float y = stary ;

                  sum_ex += x*e ;
                  sum_ey += y*e ;
                  sum_e += e ;

                  printf("   HCAL hit  %3d :  E = %7.3f   (x,y,z) = (%7.2f, %7.2f, %7.2f)\n",
                   hi, e, x, y, fcs_rec_hcalZ->at(hi) ) ;
               }
            } // hi
            float calc_cl_x = sum_ex / sum_e ;
            float calc_cl_y = sum_ey / sum_e ;
            float dx = calc_cl_x - fcs_cl_hcalX->at(hcal_ci) ;
            float dy = calc_cl_y - fcs_cl_hcalY->at(hcal_ci) ;
            printf("    HCAL calculated cluster position:  (x,y) = (%7.2f, %7.2f)   E = %7.3f   diffs  dx = %7.3f  dy = %7.3f\n",
                 calc_cl_x, calc_cl_y, sum_e, dx, dy ) ;
         }




      } // ti

      char answ = getchar() ;
      if ( answ == 'q' ) return ;

   }

} // Loop
