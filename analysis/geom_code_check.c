#define geom_code_check_cxx
#include "geom_code_check.h"
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

float ecal_rec_over_mc_scale = 4.85 ;
float hcal_rec_over_mc_scale = 68.3 ;

float ecal_dedl = ecal_rec_over_mc_scale * 0.0580 / 39. ;
float hcal_dedl = hcal_rec_over_mc_scale * 0.0200 / 84. ;

float ecal_reso2 = ecal_rec_over_mc_scale * ecal_rec_over_mc_scale * 0.0023 * 0.0023 ; // = (0.0023 *  4.85)^2 = 0.011 ^2
float hcal_reso2 = hcal_rec_over_mc_scale * hcal_rec_over_mc_scale * 0.0005  * 0.0005  ; // = (0.0005 * 68.3 )^2 = 0.034 ^2


std::map<int,float> ecal_active_path_length ;
std::map<int,float> hcal_active_path_length ;

std::map<int,float> ecal_mc_hit_energies ;
std::map<int,float> hcal_mc_hit_energies ;

std::map<int,float> ecal_rec_hit_energies ;
std::map<int,float> hcal_rec_hit_energies ;

float ecal_mc_event_sum_energy ;
float hcal_mc_event_sum_energy ;

float ecal_rec_event_sum_energy ;
float hcal_rec_event_sum_energy ;

bool verbose ;



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

void geom_code_check::load_fcs_mc_hit_energy_maps() {
   ecal_mc_hit_energies.clear() ;
   hcal_mc_hit_energies.clear() ;
   if ( verbose ) printf("\n ==== load_fcs_mc_hit_energy_maps\n") ;
   for ( int i=0; i<fcs_mc_ecalN; i++ ) {
      int id = fcs_mc_ecalVid -> at(i) % 1000 - 1 ;
      int row, col ;
      ecal_id_to_row_col( id, row, col ) ;
      float fcsx, fcsy ;
      ecal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;

      ecal_mc_hit_energies[id] = fcs_mc_ecalE->at(i) ;
      ecal_mc_event_sum_energy += fcs_mc_ecalE->at(i) ;

      float test_starx(0.) ;
      float test_stary(0.) ;
      float test_starz(0.) ;
      int ns(0) ;
      if ( fcs_mc_ecalX->at(i) > 0 ) ns = 1 ;
      if ( fcs_mc_ecalX->at(i) < 0 ) ns = 0 ;
      fcs_to_star_ecal( fcsx, fcsy, 33.3661, ns,   test_starx, test_stary, test_starz ) ;

      if ( verbose ) printf("   mc ECAL id %5d  derived row,col = (%2d, %2d) fcsx,fcsy = (%7.2f, %7.2f)\n", id, row, col, fcsx, fcsy ) ;
      if ( verbose ) printf("   mc ECAL id %5d  ttree coords = (%7.2f, %7.2f, %7.2f)  derived coords = (%7.2f, %7.2f, %7.2f)   diff (%9.4f, %9.4f, %9.4f)\n",
         id, fcs_mc_ecalX->at(i), fcs_mc_ecalY->at(i), fcs_mc_ecalZ->at(i),
         test_starx, test_stary, test_starz,
         (test_starx-fcs_mc_ecalX->at(i)),
         (test_stary-fcs_mc_ecalY->at(i)),
         (test_starz-fcs_mc_ecalZ->at(i))
         ) ;
      if ( verbose ) printf("\n") ;
   }
   for ( int i=0; i<fcs_mc_hcalN; i++ ) {
      int id = fcs_mc_hcalVid -> at(i) % 1000 - 1 ;
      int row, col ;
      hcal_id_to_row_col( id, row, col ) ;
      float fcsx, fcsy ;
      hcal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;

      hcal_mc_hit_energies[id] = fcs_mc_hcalE->at(i) ;
      hcal_mc_event_sum_energy += fcs_mc_hcalE->at(i) ;

      float test_starx(0.) ;
      float test_stary(0.) ;
      float test_starz(0.) ;
      int ns(0) ;
      if ( fcs_mc_hcalX->at(i) > 0 ) ns = 1 ;
      if ( fcs_mc_hcalX->at(i) < 0 ) ns = 0 ;
      fcs_to_star_ecal( fcsx, fcsy, 114.6741, ns,   test_starx, test_stary, test_starz ) ;

      if ( verbose ) printf("   mc HCAL id %5d  derived row,col = (%2d, %2d) fcsx,fcsy = (%7.2f, %7.2f)\n", id, row, col, fcsx, fcsy ) ;
      if ( verbose ) printf("   mc HCAL id %5d  ttree coords = (%7.2f, %7.2f, %7.2f)  derived coords = (%7.2f, %7.2f, %7.2f)   diff (%9.4f, %9.4f, %9.4f)\n",
         id, fcs_mc_hcalX->at(i), fcs_mc_hcalY->at(i), fcs_mc_hcalZ->at(i),
         test_starx, test_stary, test_starz,
         (test_starx-fcs_mc_hcalX->at(i)),
         (test_stary-fcs_mc_hcalY->at(i)),
         (test_starz-fcs_mc_hcalZ->at(i))
         ) ;
      if ( verbose ) printf("\n") ;
   }
}
//---------------

void geom_code_check::load_fcs_rec_hit_energy_maps() {
   ecal_rec_hit_energies.clear() ;
   hcal_rec_hit_energies.clear() ;
   ecal_rec_event_sum_energy = 0. ;
   hcal_rec_event_sum_energy = 0. ;
   for ( int i=0; i<fcs_rec_ecalN; i++ ) {
      int id = fcs_rec_ecalId -> at(i) ;
      ecal_rec_hit_energies[id] = fcs_rec_ecalE->at(i) ;
      ecal_rec_event_sum_energy += fcs_rec_ecalE->at(i) ;
   }
   for ( int i=0; i<fcs_rec_hcalN; i++ ) {
      int id = fcs_rec_hcalId -> at(i) ;
      hcal_rec_hit_energies[id] = fcs_rec_hcalE->at(i) ;
      hcal_rec_event_sum_energy += fcs_rec_hcalE->at(i) ;
   }
   if ( verbose ) printf("\n ==== load_fcs_rec_hit_energy_maps\n") ;
   for ( auto it = ecal_rec_hit_energies.begin(); it != ecal_rec_hit_energies.end(); it++ ) {
      int id = it->first ;
      int row, col ;
      ecal_id_to_row_col( id, row, col ) ;
      float fcsx, fcsy ;
      ecal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
      if ( verbose ) printf("  rec ECAL id %5d  derived row,col = (%2d, %2d) fcsx,fcsy = (%7.2f, %7.2f), E = %6.4f\n", id, row, col, fcsx, fcsy, it->second ) ;
   }
   for ( auto it = hcal_rec_hit_energies.begin(); it != hcal_rec_hit_energies.end(); it++ ) {
      int id = it->first ;
      int row, col ;
      hcal_id_to_row_col( id, row, col ) ;
      float fcsx, fcsy ;
      hcal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
      if ( verbose ) printf("  rec HCAL id %5d  derived row,col = (%2d, %2d) fcsx,fcsy = (%7.2f, %7.2f), E = %6.4f\n", id, row, col, fcsx, fcsy, it->second ) ;
   }
   if ( verbose ) printf("  rec sum energies:  ECAL  %7.4f,  HCAL  %7.4f\n", ecal_rec_event_sum_energy, hcal_rec_event_sum_energy ) ;

}
//---------------


//==================================================================================================================================================





void geom_code_check::Loop( bool arg_verbose, int max_events )
{
   if (fChain == 0) return;







   verbose = arg_verbose ;

   Long64_t nentries = fChain->GetEntries();
   if ( max_events > 0 && max_events < nentries ) nentries = max_events ;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if ( verbose ) {
         printf("\n\n =========== Event %lld\n", jentry ) ;
      } else {
         printf(" =========== Event %lld out of %lld\n", jentry, nentries ) ;
      }

      load_fcs_mc_hit_energy_maps() ;
      load_fcs_rec_hit_energy_maps() ;



      if ( verbose ) {
         char answ = getchar() ;
         if ( answ == 'q' ) return ;
      }




   } // jentry


} // Loop




