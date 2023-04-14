#define analysis6_cxx
#include "analysis6.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "histio.c"

float cos_angle = 0.999544 ;
float sin_angle = 0.030190 ;

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
void get_gradients( TH1* hp, float fcsx, float fcsy, float& gradx, float& grady, float dist = 1.0 ) { // distance from min is 1cm.

   gradx = 0.1 ;
   grady = 0.1 ;

   if ( hp == 0x0 ) return ;

   TAxis* ax = hp -> GetXaxis() ;
   TAxis* ay = hp -> GetYaxis() ;

   int nbins = hp -> GetNbinsX() ;
   float xmin = ax -> GetBinLowEdge(1) ;
   float xmax = ax -> GetBinUpEdge(nbins) ;
   float bw = (xmax-xmin)/nbins ;

   int bin_min_x = ax -> FindBin( fcsx ) ;
   int bin_min_y = ay -> FindBin( fcsy ) ;

   if ( verbose ) printf(" get_gradients: bin for min point x,y (%7.2f, %7.2f) is %3d, %3d\n", fcsx, fcsy, bin_min_x, bin_min_y ) ;

   int bin_low_x  = ax -> FindBin( fcsx - dist ) ;
   int bin_high_x = ax -> FindBin( fcsx + dist ) ;

   int bin_low_y  = ay -> FindBin( fcsy - dist ) ;
   int bin_high_y = ay -> FindBin( fcsy + dist ) ;

   if ( bin_low_x < 1 ) bin_low_x = 1 ;
   if ( bin_high_x > nbins ) bin_high_x = nbins ;

   if ( bin_low_y < 1 ) bin_low_y = 1 ;
   if ( bin_high_y > nbins ) bin_high_y = nbins ;

   float val_min = hp -> GetBinContent( bin_min_x, bin_min_y ) ;

   float val_low_x = hp -> GetBinContent( bin_low_x, bin_min_y ) ;
   float val_high_x = hp -> GetBinContent( bin_high_x, bin_min_y ) ;

   float val_low_y = hp -> GetBinContent( bin_min_x, bin_low_y ) ;
   float val_high_y = hp -> GetBinContent( bin_min_x, bin_high_y ) ;

   if ( val_low_x  < 0 ) val_low_x  = val_min ;
   if ( val_high_x < 0 ) val_high_x = val_min ;

   if ( val_low_y  < 0 ) val_low_y  = val_min ;
   if ( val_high_y < 0 ) val_high_y = val_min ;


  //--- use the LOWER of the two differences (the weaker side).

   if ( fabs( val_low_x - val_min ) > fabs( val_high_x - val_min ) ) {
      gradx = val_high_x - val_min ;
   } else {
      gradx = val_low_x - val_min ;
   }

   if ( fabs( val_low_y - val_min ) > fabs( val_high_y - val_min ) ) {
      grady = val_high_y - val_min ;
   } else {
      grady = val_low_y - val_min ;
   }
   if ( gradx < 0. ) gradx = 0.1 ;
   if ( grady < 0. ) grady = 0.1 ;

   if ( verbose ) printf(" get_gradients:  min coords (%6.2f,%6.2f) , min val  %7.2f,  x: %7.2f, %7.2f (%7.2f, %10.1f),  y: %7.2f, %7.2f (%7.2f, %10.1f)\n",
       fcsx, fcsy, val_min, val_low_x, val_high_x, gradx, gradx/val_min, val_low_y, val_high_y, grady, grady/val_min ) ;

}


//---------------
void star_to_fcs_ecal( float starx, float stary, float starz,
                       float& fcsx, float& fcsy, float& fcsz ) {
    if ( starx > 0 ) {
       fcsx = starx - 17.40 ;
       fcsy = stary ;
       fcsz = starz - 710.16 ;
       float tmpx = fcsx ;
       float tmpz = fcsz ;
       fcsx =      cos_angle * tmpx - 1 * sin_angle * tmpz ;
       fcsz =      cos_angle * tmpz + 1 * sin_angle * tmpx ;
    } else {
       fcsx = starx + 17.40 ;
       fcsy = stary  ;
       fcsz = starz - 710.16 ;
       float tmpx = fcsx ;
       float tmpz = fcsz ;
       fcsx =      cos_angle * tmpx + 1 * sin_angle * tmpz ;
       fcsz =      cos_angle * tmpz - 1 * sin_angle * tmpx ;
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
void do_dz_steps( float fcsx, float fcsy, float fcsz, float dxdz, float dydz,
                  float& iso1_chi2_ecal,
                  float& iso1_chi2_hcal,
                  float& iso2_chi2_ecal,
                  float& iso2_chi2_hcal,
                  float& energy_chi2_ecal,
                  float& energy_chi2_hcal,
                  bool use_mc=true  ) {

   iso1_chi2_ecal = 1e9 ;
   iso1_chi2_hcal = 1e9 ;
   iso2_chi2_ecal = 1e9 ;
   iso2_chi2_hcal = 1e9 ;
   energy_chi2_ecal = 1e9 ;
   energy_chi2_hcal = 1e9 ;

   std::map<int,float>* ecal_hit_energies;
   std::map<int,float>* hcal_hit_energies;
   if ( use_mc ) {
      ecal_hit_energies = &ecal_mc_hit_energies ;
      hcal_hit_energies = &hcal_mc_hit_energies ;
   } else {
      ecal_hit_energies = &ecal_rec_hit_energies ;
      hcal_hit_energies = &hcal_rec_hit_energies ;
   }


   if ( verbose ) printf("\n\n ===== top of do_dz_steps\n") ;

   if ( verbose ) printf("  starting point:  fcs x,y,z = (%7.2f, %7.2f, %7.2f),  dxdz = %7.4f, dydz = %7.4f\n",
      fcsx, fcsy, fcsz, dxdz, dydz ) ;

   ecal_active_path_length.clear() ;
   hcal_active_path_length.clear() ;

   float dz = 1.0 ;
   float dl = dz * sqrt( 1. + dxdz*dxdz + dydz*dydz ) ;

   int nsteps = TMath::Nint( hcal_fcs_zmax / dz ) ;

   if ( verbose ) printf("  fcs zmax = %7.2f,  will do %d steps of %.3f\n", hcal_fcs_zmax, nsteps, dz ) ;

   int row, col ;
   ecal_fcsxy_to_row_col( fcsx, fcsy, row, col ) ;
   int id = ecal_row_col_to_id( row, col ) ;

   if ( verbose ) printf("  starting in ECAL cell row %d, col %d, id %d\n", row, col, id ) ;

   float x = fcsx ;
   float y = fcsy ;
   float z = fcsz ;

   for ( int si=0; si<nsteps; si++ ) {

      x = x + dz * dxdz ;
      y = y + dz * dydz ;
      z = z + dz ;
      //////printf("   dz step %5d :  x,y,z (%7.2f, %7.2f, %7.2f)\n", si, x, y, z ) ;

      if ( z > ecal_fcs_zmin && z < ecal_fcs_zmax ) {

         ecal_fcsxy_to_row_col( x, y, row, col ) ;
         id = ecal_row_col_to_id( row, col ) ;
         if ( ecal_active_path_length.find( id ) == ecal_active_path_length.end() ) {
            if ( verbose ) printf("  ECAL Adding id =  %5d to active_path_length at x,y,z = (%7.2f, %7.2f, %7.2f)\n", id, x, y, z ) ;
            ecal_active_path_length[id] = dl ;
         } else {
            ecal_active_path_length[id] = ecal_active_path_length[id] + dl ;
         }

      } else if ( z > hcal_fcs_zmin && z < hcal_fcs_zmax ) {

         hcal_fcsxy_to_row_col( x, y, row, col ) ;
         id = hcal_row_col_to_id( row, col ) ;
         if ( hcal_active_path_length.find( id ) == hcal_active_path_length.end() ) {
            if ( verbose ) printf("  HCAL Adding id =  %5d to active_path_length at x,y,z = (%7.2f, %7.2f, %7.2f)\n", id, x, y, z ) ;
            hcal_active_path_length[id] = dl ;
         } else {
            hcal_active_path_length[id] = hcal_active_path_length[id] + dl ;
         }

      }

   } // si

   if ( verbose ) printf("\n") ;
   for ( auto it = ecal_active_path_length.begin(); it != ecal_active_path_length.end(); it++ ) {
      if ( verbose ) printf("  ECAL id %5d path length = %7.2f\n", it->first, it->second ) ;
   }
   for ( auto it = hcal_active_path_length.begin(); it != hcal_active_path_length.end(); it++ ) {
      if ( verbose ) printf("  HCAL id %5d path length = %7.2f\n", it->first, it->second ) ;
   }




  //--- compute the isolation

   std::map<int,float> ecal_neighbor_energies1 ;
   std::map<int,float> hcal_neighbor_energies1 ;
   std::map<int,float> ecal_neighbor_energies2 ;
   std::map<int,float> hcal_neighbor_energies2 ;

   int n_cell_iso = 2 ;

   for ( auto it = ecal_active_path_length.begin(); it != ecal_active_path_length.end(); it++ ) {
      int hit_id = it->first ;
      int hit_row, hit_col ;
      ecal_id_to_row_col( hit_id, hit_row, hit_col ) ;
      for ( int neighbor_row = hit_row-n_cell_iso; neighbor_row <= hit_row+n_cell_iso; neighbor_row++ ) {
         for ( int neighbor_col = hit_col-n_cell_iso; neighbor_col <= hit_col+n_cell_iso; neighbor_col++ ) {
            if ( neighbor_row == hit_row && neighbor_col == hit_col ) continue ; // don't count the hit!
            int neighbor_id = ecal_row_col_to_id( neighbor_row, neighbor_col ) ;
            if ( ecal_active_path_length.find( neighbor_id ) != ecal_active_path_length.end() ) continue ; // don't include any other hits with length>0.
            if ( ecal_hit_energies->find( neighbor_id ) != ecal_hit_energies->end() ) {
               float energy = (*ecal_hit_energies)[ neighbor_id ] ;
               if ( use_mc ) energy = ecal_rec_over_mc_scale * energy ;
               ecal_neighbor_energies2[ neighbor_id ] = energy ;
               if (abs(neighbor_row-hit_row)<=1 ) ecal_neighbor_energies1[ neighbor_id ] = energy ;
               if ( verbose ) printf("   ECAL isolation:  hit id %5d (%2d,%2d)    neighbor %5d (%2d,%2d)    energy = %7.3f\n",
                   hit_id, hit_row, hit_col, neighbor_id, neighbor_row, neighbor_col, energy ) ;
            }
         }
      }
   }
   for ( auto it = hcal_active_path_length.begin(); it != hcal_active_path_length.end(); it++ ) {
      int hit_id = it->first ;
      int hit_row, hit_col ;
      hcal_id_to_row_col( hit_id, hit_row, hit_col ) ;
      for ( int neighbor_row = hit_row-n_cell_iso; neighbor_row <= hit_row+n_cell_iso; neighbor_row++ ) {
         for ( int neighbor_col = hit_col-n_cell_iso; neighbor_col <= hit_col+n_cell_iso; neighbor_col++ ) {
            if ( neighbor_row == hit_row && neighbor_col == hit_col ) continue ; // don't count the hit!
            int neighbor_id = hcal_row_col_to_id( neighbor_row, neighbor_col ) ;
            if ( hcal_active_path_length.find( neighbor_id ) != hcal_active_path_length.end() ) continue ; // don't include any other hits with length>0.
            if ( hcal_hit_energies->find( neighbor_id ) != hcal_hit_energies->end() ) {
               float energy = (*hcal_hit_energies)[ neighbor_id ] ;
               if ( use_mc ) energy = hcal_rec_over_mc_scale * energy ;
               hcal_neighbor_energies2[ neighbor_id ] = energy ;
               if (abs(neighbor_row-hit_row)<=1 ) hcal_neighbor_energies1[ neighbor_id ] = energy ;
               if ( verbose ) printf("   HCAL isolation:  hit id %5d (%2d,%2d)    neighbor %5d (%2d,%2d)    energy = %7.3f\n",
                   hit_id, hit_row, hit_col, neighbor_id, neighbor_row, neighbor_col, energy ) ;
            }
         }
      }
   }

   float ecal_iso1_e2_sum(0.) ;
   float hcal_iso1_e2_sum(0.) ;
   for ( auto it = ecal_neighbor_energies1.begin(); it != ecal_neighbor_energies1.end(); it++ ) {
      float energy = it->second ;
      if ( verbose ) printf("   calc ECAL iso1 sum:  id = %5d,  E = %7.4f\n", it->first, energy ) ;
      ecal_iso1_e2_sum += energy*energy ;
   }
   for ( auto it = hcal_neighbor_energies1.begin(); it != hcal_neighbor_energies1.end(); it++ ) {
      float energy = it->second ;
      if ( verbose ) printf("   calc HCAL iso1 sum:  id = %5d,  E = %7.4f\n", it->first, energy ) ;
      hcal_iso1_e2_sum += energy*energy ;
   }

   float ecal_iso2_e2_sum(0.) ;
   float hcal_iso2_e2_sum(0.) ;
   for ( auto it = ecal_neighbor_energies2.begin(); it != ecal_neighbor_energies2.end(); it++ ) {
      float energy = it->second ;
      if ( verbose ) printf("   calc ECAL iso2 sum:  id = %5d,  E = %7.4f\n", it->first, energy ) ;
      ecal_iso2_e2_sum += energy*energy ;
   }
   for ( auto it = hcal_neighbor_energies2.begin(); it != hcal_neighbor_energies2.end(); it++ ) {
      float energy = it->second ;
      if ( verbose ) printf("   calc HCAL iso2 sum:  id = %5d,  E = %7.4f\n", it->first, energy ) ;
      hcal_iso2_e2_sum += energy*energy ;
   }


   ////////iso_chi2 = ecal_iso_e2_sum / ecal_reso2 + hcal_iso_e2_sum / hcal_reso2 ;
   iso1_chi2_ecal = ecal_iso1_e2_sum / ecal_reso2 ;
   iso1_chi2_hcal = hcal_iso1_e2_sum / hcal_reso2 ;
   iso2_chi2_ecal = ecal_iso2_e2_sum / ecal_reso2 ;
   iso2_chi2_hcal = hcal_iso2_e2_sum / hcal_reso2 ;
   if ( verbose ) printf("   Isolation chi2 = %7.2f  (ecal %7.2f, hcal %7.2f)\n", iso1_chi2_ecal+iso1_chi2_hcal, iso1_chi2_ecal, iso1_chi2_hcal ) ;




   std::set<int> ecal_ids ;
   std::set<int> hcal_ids ;

   for ( auto it = ecal_active_path_length.begin(); it != ecal_active_path_length.end(); it++ ) { ecal_ids.insert( it->first ) ; }
   for ( auto it = hcal_active_path_length.begin(); it != hcal_active_path_length.end(); it++ ) { hcal_ids.insert( it->first ) ; }






   /////////energy_chi2 = 0. ;
   energy_chi2_ecal = 0. ;
   energy_chi2_hcal = 0. ;

   if ( verbose ) printf("\n computing energy chi2.\n" ) ;
   float ecal_on_track_esum(0.) ;
   for ( auto it = ecal_ids.begin(); it != ecal_ids.end(); it++ ) {

      float energy = 0. ;
      float length = 0. ;
      if ( ecal_hit_energies->find( *it ) != ecal_hit_energies->end() ) energy = (*ecal_hit_energies)[*it] ;
      if ( ecal_active_path_length.find( *it ) != ecal_active_path_length.end() ) length = ecal_active_path_length[*it] ;

      if ( use_mc ) energy = ecal_rec_over_mc_scale * energy ;

      float pred_e = length * ecal_dedl ;
      float diff2 = (energy - pred_e)*(energy - pred_e) ;
      //-- only include here if length is non-zero.
      if ( length > 0 ) {
         energy_chi2_ecal += diff2 / ecal_reso2 ;
         ecal_on_track_esum += energy ;
         if ( verbose ) printf("  ECAL :  id = %5d   L = %6.1f  E = %7.3f  pred E = %7.3f    diff2/reso2 = %7.3f\n", *it, length, energy, pred_e, diff2 / ecal_reso2 ) ;
      }

   }
   float hcal_on_track_esum(0.) ;
   for ( auto it = hcal_ids.begin(); it != hcal_ids.end(); it++ ) {

      float energy = 0. ;
      float length = 0. ;
      if ( hcal_hit_energies->find( *it ) != hcal_hit_energies->end() ) energy = (*hcal_hit_energies)[*it] ;
      if ( hcal_active_path_length.find( *it ) != hcal_active_path_length.end() ) length = hcal_active_path_length[*it] ;

      if ( use_mc ) energy = hcal_rec_over_mc_scale * energy ;

      float pred_e = length * hcal_dedl ;
      float diff2 = (energy - pred_e)*(energy - pred_e) ;
      //-- only include here if length is non-zero.
      if ( length > 0 ) {
         energy_chi2_hcal += diff2 / hcal_reso2 ;
         hcal_on_track_esum += energy ;
         if ( verbose ) printf("  HCAL :  id = %5d   L = %6.1f  E = %7.3f  pred E = %7.3f    diff2/reso2 = %7.3f\n", *it, length, energy, pred_e, diff2 / hcal_reso2 ) ;
      }

   }

   /////////float global_chi2 = iso_chi2 + energy_chi2 ;
   float global_chi2 = iso1_chi2_ecal + iso1_chi2_hcal + energy_chi2_ecal + energy_chi2_hcal ;

   //-- set it to a negative number of not enough energy to be valid.
   if ( ecal_on_track_esum < 0.05 ) { energy_chi2_ecal = -10000. ; iso1_chi2_ecal = -10000. ; iso2_chi2_ecal = -10000. ; }
   if ( hcal_on_track_esum < 0.10 ) { energy_chi2_hcal = -10000. ; iso1_chi2_hcal = -10000. ; iso2_chi2_hcal = -10000. ; }

   if ( verbose ) printf("  chi2 =  energy %.2f (ecal %.2f, hcal %.2f) + iso %.2f (ecal %.2f, hcal %.2f) = %.2f\n",
       energy_chi2_ecal+energy_chi2_hcal, energy_chi2_ecal, energy_chi2_hcal,
       iso1_chi2_ecal+iso1_chi2_hcal, iso1_chi2_ecal, iso1_chi2_hcal,
       global_chi2 ) ;
   if ( verbose ) printf("\n") ;



} // do_dz_steps
//---------------

void analysis6::load_fcs_mc_hit_energy_maps() {
   ecal_mc_hit_energies.clear() ;
   hcal_mc_hit_energies.clear() ;
   ecal_mc_event_sum_energy = 0. ;
   hcal_mc_event_sum_energy = 0. ;
   for ( int i=0; i<fcs_mc_ecalN; i++ ) {
      int id = fcs_mc_ecalVid -> at(i) % 1000 - 1 ;
      ecal_mc_hit_energies[id] = fcs_mc_ecalE->at(i) ;
      ecal_mc_event_sum_energy += fcs_mc_ecalE->at(i) ;
   }
   for ( int i=0; i<fcs_mc_hcalN; i++ ) {
      int id = fcs_mc_hcalVid -> at(i) % 1000 - 1 ;
      hcal_mc_hit_energies[id] = fcs_mc_hcalE->at(i) ;
      hcal_mc_event_sum_energy += fcs_mc_hcalE->at(i) ;
   }
   if ( verbose ) printf("\n ==== load_fcs_mc_hit_energy_maps\n") ;
   for ( auto it = ecal_mc_hit_energies.begin(); it != ecal_mc_hit_energies.end(); it++ ) {
      int id = it->first ;
      int row, col ;
      ecal_id_to_row_col( id, row, col ) ;
      float fcsx, fcsy ;
      ecal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
      if ( verbose ) printf("   mc ECAL id %5d  derived row,col = (%2d, %2d) fcsx,fcsy = (%7.2f, %7.2f), E = %6.4f\n", id, row, col, fcsx, fcsy, it->second ) ;
   }
   for ( auto it = hcal_mc_hit_energies.begin(); it != hcal_mc_hit_energies.end(); it++ ) {
      int id = it->first ;
      int row, col ;
      hcal_id_to_row_col( id, row, col ) ;
      float fcsx, fcsy ;
      hcal_row_col_to_fcsxy( row, col, fcsx, fcsy ) ;
      if ( verbose ) printf("   mc HCAL id %5d  derived row,col = (%2d, %2d) fcsx,fcsy = (%7.2f, %7.2f), E = %6.4f\n", id, row, col, fcsx, fcsy, it->second ) ;
   }
   if ( verbose ) printf("   mc sum energies:  ECAL  %7.4f,  HCAL  %7.4f\n", ecal_mc_event_sum_energy, hcal_mc_event_sum_energy ) ;
}
//---------------

void analysis6::load_fcs_rec_hit_energy_maps() {
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

float get_th2_min( TH2F* hp, float& minx, float& miny ) {

   minx = 99. ;
   miny = 99. ;
   if ( hp == 0x0 ) return 1e9 ;

  //-- give preference to the point closest to the center.

   int nbinsx = hp -> GetNbinsX() ;
   int nbinsy = hp -> GetNbinsY() ;
   int midbin = nbinsx/2+1 ;


   //float global_minval = hp -> GetBinContent( hp->GetMinimumBin() ) ;
   float global_minval = 1e12 ;
   for ( int xbi=1; xbi<=nbinsx; xbi++ ) {
      for ( int ybi=1; ybi<=nbinsy; ybi++ ) {
         float val = hp -> GetBinContent( xbi, ybi ) ;
         if ( val < 0 ) continue ;
         if ( val > 0 && val < global_minval ) { global_minval = val ; }
      }
   }
   if ( verbose ) printf(" get_th2_min global min val:  %7.2f\n", global_minval ) ;

   float mindist2 = 1e9 ;
   int minxbi = -1 ;
   int minybi = -1 ;
   float minval = 1e9 ;
   for ( int xbi=1; xbi<=nbinsx; xbi++ ) {
      for ( int ybi=1; ybi<=nbinsy; ybi++ ) {
         float val = hp -> GetBinContent( xbi, ybi ) ;
         if ( val < 0 ) continue ;
         if ( fabs(val - global_minval) < 3.9) {
            float dist2 = (xbi-midbin)*(xbi-midbin) + (ybi-midbin)*(ybi-midbin) ;
            if ( dist2 < mindist2 ) {
               minval = val ;
               minxbi = xbi ;
               minybi = ybi ;
               mindist2 = dist2 ;
               float bin_center_x = hp -> GetXaxis() -> GetBinCenter( xbi ) ;
               float bin_center_y = hp -> GetYaxis() -> GetBinCenter( ybi ) ;
               if ( verbose ) printf(" get_th2_min b:  %3d, %3d  (%7.2f, %7.2f),  val = %7.4f, dist2 = %7.1f\n", 
                   xbi, ybi, bin_center_x, bin_center_y, val, dist2 ) ;
            }
         }
      } // ybi
   } // xbi


   minx = hp -> GetXaxis() -> GetBinCenter( minxbi ) ;
   miny = hp -> GetYaxis() -> GetBinCenter( minybi ) ;

   if ( verbose ) printf(" get_th2_min final:  %3d, %3d,  val = %7.4f, dist2 = %7.1f\n", minxbi, minybi, minval, mindist2 ) ;

   return minval ;

} // get_th2_min

//==================================================================================================================================================





void analysis6::Loop( bool arg_verbose, int max_events )
{
   if (fChain == 0) return;

   float scan_dx = 0.2 ;
   //////////int   scan_npoints = 51 ;
   int   scan_npoints = 21 ;

   float scan_dy = scan_dx ;
   float scan_xmax = scan_dx * (scan_npoints-1)/2 ;
   float scan_ymax = scan_xmax ;

   bool use_mc = false ;

   TH2F* h_chi2 = new TH2F( "h_chi2", "chi2", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx, scan_npoints, -1*scan_ymax-0.5*scan_dx, scan_ymax+0.5*scan_dx ) ;
   TH2F* h_chi2_iso = new TH2F( "h_chi2_iso", "iso chi2", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx, scan_npoints, -1*scan_ymax-0.5*scan_dx, scan_ymax+0.5*scan_dx ) ;
   TH2F* h_chi2_energy = new TH2F( "h_chi2_energy", "energy chi2", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx, scan_npoints, -1*scan_ymax-0.5*scan_dx, scan_ymax+0.5*scan_dx ) ;

   TH1F* h_log10chi2_all = new TH1F( "h_log10chi2_all", "log10 chi2", 80, -3., 8. ) ;
   TH1F* h_log10chi2_mu = new TH1F( "h_log10chi2_mu", "log10 chi2, muons", 80, -3., 8. ) ;
   TH1F* h_log10chi2_notmu = new TH1F( "h_log10chi2_notmu", "log10 chi2, not muons", 80, -3., 8. ) ;

   TH1F* h_log10chi2_energy_all = new TH1F( "h_log10chi2_energy_all", "log10 chi2, energy", 80, -3., 8. ) ;
   TH1F* h_log10chi2_energy_mu = new TH1F( "h_log10chi2_energy_mu", "log10 chi2, energy, muons", 80, -3., 8. ) ;
   TH1F* h_log10chi2_energy_notmu = new TH1F( "h_log10chi2_energy_notmu", "log10 chi2, energy, not muons", 80, -3., 8. ) ;

   TH1F* h_log10chi2_iso_all = new TH1F( "h_log10chi2_iso_all", "log10 chi2, iso", 80, -3., 8. ) ;
   TH1F* h_log10chi2_iso_mu = new TH1F( "h_log10chi2_iso_mu", "log10 chi2, iso, muons", 80, -3., 8. ) ;
   TH1F* h_log10chi2_iso_notmu = new TH1F( "h_log10chi2_iso_notmu", "log10 chi2, iso, not muons", 80, -3., 8. ) ;

   TH2F* h_log10chi2_energy_vs_iso_all = new TH2F( "h_log10chi2_energy_vs_iso_all", "log10 chi2, energy vs iso",  80, -3., 8., 80, -3., 8. ) ;
   TH2F* h_log10chi2_energy_vs_iso_mu = new TH2F( "h_log10chi2_energy_vs_iso_mu", "log10 chi2, energy vs iso, muons",  80, -3., 8., 80, -3., 8. ) ;
   TH2F* h_log10chi2_energy_vs_iso_notmu = new TH2F( "h_log10chi2_energy_vs_iso_notmu", "log10 chi2, energy vs iso, not muons",  80, -3., 8., 80, -3., 8. ) ;

   TH1F* h_dist_from_trk_all = new TH1F( "h_dist_from_trk_all", "distance from track, all", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx  ) ;
   TH1F* h_dist_from_trk_mu = new TH1F( "h_dist_from_trk_mu", "distance from track, muons", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx  ) ;
   TH1F* h_dist_from_trk_notmu = new TH1F( "h_dist_from_trk_notmu", "distance from track, not muons", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx  ) ;


   h_log10chi2_mu -> SetLineColor(4) ;
   h_log10chi2_energy_mu -> SetLineColor(4) ;
   h_log10chi2_iso_mu -> SetLineColor(4) ;
   h_dist_from_trk_mu -> SetLineColor(4) ;

   h_log10chi2_notmu -> SetLineColor(2) ;
   h_log10chi2_energy_notmu -> SetLineColor(2) ;
   h_log10chi2_iso_notmu -> SetLineColor(2) ;
   h_dist_from_trk_notmu -> SetLineColor(2) ;


   TCanvas* can1 = new TCanvas( "can1", "", 50, 50, 1100, 1300 ) ;
   can1 -> Divide(3,3) ;


   TFile* tf_output = new TFile( "analysis-tree.root", "RECREATE" ) ;

   TTree* tt_out = fChain -> CloneTree( 0 ) ;


   vector<float> rc_chi2_energy_ecal ;
   vector<float> rc_chi2_energy_hcal ;
   vector<float> rc_chi2_energy ;
   vector<float> rc_chi2_iso1_ecal ;
   vector<float> rc_chi2_iso1_hcal ;
   vector<float> rc_chi2_iso1 ;
   vector<float> rc_chi2_iso2_ecal ;
   vector<float> rc_chi2_iso2_hcal ;
   vector<float> rc_chi2_iso2 ;
   vector<float> rc_dx ;
   vector<float> rc_dy ;
   vector<float> rc_dx_gradient ;
   vector<float> rc_dy_gradient ;
   vector<float> rc_energy_ecal ;
   vector<float> rc_energy_hcal ;
   vector<float> rc_Pt ;
   vector<float> rc_Eta ;
   vector<int>   rc_nhit_ecal ;
   vector<int>   rc_nhit_hcal ;
   vector<int>   rc_mcPid ;
   vector<int>   rc_rci ;
   vector<bool>  rc_ismu ;
   tt_out -> Branch( "rc_chi2_energy_ecal", &rc_chi2_energy_ecal ) ;
   tt_out -> Branch( "rc_chi2_energy_hcal", &rc_chi2_energy_hcal ) ;
   tt_out -> Branch( "rc_chi2_energy", &rc_chi2_energy ) ;
   tt_out -> Branch( "rc_chi2_iso1_ecal", &rc_chi2_iso1_ecal ) ;
   tt_out -> Branch( "rc_chi2_iso1_hcal", &rc_chi2_iso1_hcal ) ;
   tt_out -> Branch( "rc_chi2_iso1", &rc_chi2_iso1 ) ;
   tt_out -> Branch( "rc_chi2_iso2_ecal", &rc_chi2_iso2_ecal ) ;
   tt_out -> Branch( "rc_chi2_iso2_hcal", &rc_chi2_iso2_hcal ) ;
   tt_out -> Branch( "rc_chi2_iso2", &rc_chi2_iso2 ) ;
   tt_out -> Branch( "rc_dx", &rc_dx ) ;
   tt_out -> Branch( "rc_dy", &rc_dy ) ;
   tt_out -> Branch( "rc_dx_gradient", &rc_dx_gradient ) ;
   tt_out -> Branch( "rc_dy_gradient", &rc_dy_gradient ) ;
   tt_out -> Branch( "rc_energy_ecal", &rc_energy_ecal ) ;
   tt_out -> Branch( "rc_energy_hcal", &rc_energy_hcal ) ;
   tt_out -> Branch( "rc_Pt", &rc_Pt ) ;
   tt_out -> Branch( "rc_Eta", &rc_Eta ) ;
   tt_out -> Branch( "rc_nhit_ecal", &rc_nhit_ecal ) ;
   tt_out -> Branch( "rc_nhit_hcal", &rc_nhit_hcal ) ;
   tt_out -> Branch( "rc_mcPid", &rc_mcPid ) ;
   tt_out -> Branch( "rc_rci", &rc_rci ) ;
   tt_out -> Branch( "rc_ismu", &rc_ismu ) ;


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
         printf(" =========== Event %lld\n", jentry ) ;
      }

      load_fcs_mc_hit_energy_maps() ;
      load_fcs_rec_hit_energy_maps() ;


      for ( int rci=0; rci<rcEta->size(); rci++ ) {

         if ( rcPt->at(rci) < 0.2 ) continue ;

         int tpind(-1) ;
         for ( int i=0; i<tprojIdT->size(); i++ ) {
            if ( tprojIdT->at(i) != rcTrackIdT->at(rci) ) continue ;
            if ( tprojIdD -> at(i) != 4 ) continue ;
            tpind = i ;
            break ;
         } // i
         int mcind(-1) ;
         bool mc_muon_match(false) ;
         if ( rcTrackId -> at(rci) >= 1 ) {
            mcind = rcTrackId -> at(rci) - 1 ;
            if ( abs( mcPid->at(mcind) ) == 5 || abs( mcPid->at(mcind) ) == 6 ) mc_muon_match = true ;
         }



         if ( verbose ) {
            printf("\n\n\n +++++++++++++++++++++++++++++++++++++++++++\n") ;
            printf(" TRACK %3d:  Eta = %7.3f, Phi = %7.3f, Pt = %7.3f,  q = %2d,  Id = %3d, IdT = %3d\n",
                rci, rcEta->at(rci), rcPhi->at(rci), rcPt->at(rci), rcCharge->at(rci), rcTrackId->at(rci), rcTrackIdT->at(rci) ) ;
            if ( tpind >= 0 ) {
               float pt = sqrt( tprojPx->at(tpind) * tprojPx->at(tpind) + tprojPy->at(tpind) * tprojPy->at(tpind) ) ;
               float theta = atan2( pt, tprojPz->at(tpind) ) ;
               float phi = atan2( tprojPy->at(tpind), tprojPx->at(tpind) ) ;
               float eta = -1 * log( tan( theta/2. ) ) ;
               printf(" tproj %3d: Eta = %7.3f, Phi = %7.3f, Pt = %7.3f  idt = %5d,  x,y,z = (%7.2f, %7.2f, %7.2f)\n",
                  tpind, eta, phi, pt, tprojIdT->at(tpind), tprojX->at(tpind), tprojY->at(tpind), tprojZ->at(tpind) ) ;
            } else {
               printf("    *** can't find track projection?\n") ;
            }
            if ( mcind >= 0 ) {
               float deta = mcEta->at(mcind) - rcEta->at(rci) ;
               float dphi = mcPhi->at(mcind) - rcPhi->at(rci) ;
               if ( dphi > 3.14159265 ) dphi = dphi - 2*3.14159265 ;
               if ( dphi <-3.14159265 ) dphi = dphi + 2*3.14159265 ;
               float dr = sqrt( deta*deta + dphi*dphi ) ;
                printf("    MC %3d: Eta = %7.3f, Phi = %7.3f, Pt = %7.3f,  q = %2d,  Id = %3d, PID = %8d, dr = %7.3f",
                  mcind, mcEta->at(mcind), mcPhi->at(mcind), mcPt->at(mcind), mcCharge->at(mcind), mcind, mcPid->at(mcind), dr ) ;
               if ( mc_muon_match ) { printf(" ** muon") ; }
               printf("\n") ;
            } else {
               printf("    *** no MC match?\n") ;
            }
            printf("\n") ;
         }


         if ( mcind < 0 ) {
            printf("  %lld , trk %d : no MC ind.\n", jentry, rci ) ;
            continue ;
         }
         if ( tpind < 0 ) {
            printf("  %lld , trk %d : no projection ind.\n", jentry, rci ) ;
            continue ;
         }



         float ecal_trkproj_starx = tprojX->at(tpind) ;
         float ecal_trkproj_stary = tprojY->at(tpind) ;
         float ecal_trkproj_starz = tprojZ->at(tpind) ;

         float ecal_trkproj_fcsx ;
         float ecal_trkproj_fcsy ;
         float ecal_trkproj_fcsz ;
         star_to_fcs_ecal( ecal_trkproj_starx, ecal_trkproj_stary, ecal_trkproj_starz, ecal_trkproj_fcsx, ecal_trkproj_fcsy, ecal_trkproj_fcsz ) ;

         int trkproj_ecal_row ;
         int trkproj_ecal_col ;
         ecal_fcsxy_to_row_col( ecal_trkproj_fcsx, ecal_trkproj_fcsy, trkproj_ecal_row, trkproj_ecal_col ) ;
         if ( verbose ) printf("  Track projection, star x,y,z (%7.2f, %7.2f, %7.2f),  FCS x,y,z (%7.2f, %7.2f, %7.2f)   FCS ecal row = %2d, col = %2d\n",
              ecal_trkproj_starx, ecal_trkproj_stary, ecal_trkproj_starz, ecal_trkproj_fcsx, ecal_trkproj_fcsy, ecal_trkproj_fcsz,  trkproj_ecal_row, trkproj_ecal_col ) ;


         float trk_star_px = tprojPx->at(tpind) ;
         float trk_star_py = tprojPy->at(tpind) ;
         float trk_star_pz = tprojPz->at(tpind) ;

         float trk_fcs_py = trk_star_py ;

         float trk_fcs_px ;
         float trk_fcs_pz ;

         if ( ecal_trkproj_starx > 0 ) {
            trk_fcs_px = cos_angle * trk_star_px - 1 * sin_angle * trk_star_pz ;
            trk_fcs_pz = cos_angle * trk_star_pz + 1 * sin_angle * trk_star_px ;
         } else {
            trk_fcs_px = cos_angle * trk_star_px + 1 * sin_angle * trk_star_pz ;
            trk_fcs_pz = cos_angle * trk_star_pz - 1 * sin_angle * trk_star_px ;
         }
         if ( verbose ) printf("\n  track before rotation  px, pz = (%7.3f, %7.3f)   after rotation (%7.3f, %7.3f)\n",
            trk_star_px, trk_star_pz, trk_fcs_px, trk_fcs_pz ) ;

         float trk_dxdz = trk_fcs_px / trk_fcs_pz ;
         float trk_dydz = trk_fcs_py / trk_fcs_pz ;

       //--- Need to add x sign flip and rotation about y for dxdz.
         if ( ecal_trkproj_starx < 0 ) {
            if ( verbose ) printf("  Changing sign of dxdz from %.4f to %.4f\n", trk_dxdz, -1*trk_dxdz ) ;
            trk_dxdz = -1*trk_dxdz ;
         }


         float iso1_chi2 ;
         float iso1_chi2_ecal, iso1_chi2_hcal ;
         float iso2_chi2 ;
         float iso2_chi2_ecal, iso2_chi2_hcal ;
         float energy_chi2 ;
         float energy_chi2_ecal, energy_chi2_hcal ;



         verbose = false ;

         h_chi2 -> Reset() ;
         h_chi2_iso -> Reset() ;
         h_chi2_energy -> Reset() ;

         float best_chi2(99999999.) ;
         float min_x, min_y ;
         int necal(0), nhcal(0) ;
         int necal_enz(0), nhcal_enz(0) ;

         float best_chi2_dx(0.) ;
         float best_chi2_dy(0.) ;
         float worst_chi2(0.) ;
         for ( int sxi = 0; sxi<scan_npoints; sxi++ ) {

            float dx = -1 * scan_xmax + sxi*scan_dx ;
            float scan_x = ecal_trkproj_fcsx + dx ;

            for ( int syi = 0; syi<scan_npoints; syi++ ) {

               float dy = -1 * scan_ymax + syi*scan_dy ;
               float scan_y = ecal_trkproj_fcsy + dy ;

               do_dz_steps( scan_x, scan_y, ecal_trkproj_fcsz, trk_dxdz, trk_dydz,
                            iso1_chi2_ecal, iso1_chi2_hcal,
                            iso2_chi2_ecal, iso2_chi2_hcal,
                            energy_chi2_ecal, energy_chi2_hcal, use_mc ) ;

               iso1_chi2 = iso1_chi2_ecal + iso1_chi2_hcal ;
               iso2_chi2 = iso2_chi2_ecal + iso2_chi2_hcal ;
               energy_chi2 = energy_chi2_ecal + energy_chi2_hcal ;

             //-- check for insufficient energy in either and set it for both if found.
               if ( iso1_chi2_ecal < 0 ) iso1_chi2 = iso1_chi2_ecal ;
               if ( iso1_chi2_hcal < 0 ) iso1_chi2 = iso1_chi2_hcal ;
               if ( iso2_chi2_ecal < 0 ) iso2_chi2 = iso2_chi2_ecal ;
               if ( iso2_chi2_hcal < 0 ) iso2_chi2 = iso2_chi2_hcal ;
               if ( energy_chi2_ecal < 0 ) energy_chi2 = energy_chi2_ecal ;
               if ( energy_chi2_hcal < 0 ) energy_chi2 = energy_chi2_hcal ;


               float chi2 = iso1_chi2 + energy_chi2 ;
               h_chi2 -> Fill( dx, dy, chi2 ) ;
               h_chi2_iso -> Fill( dx, dy, iso1_chi2 ) ;
               h_chi2_energy -> Fill( dx, dy, energy_chi2 ) ;

               if ( chi2 > 0 && chi2 < best_chi2 ) {
                  best_chi2 = chi2 ;
                  best_chi2_dx = dx ;
                  best_chi2_dy = dy ;
               }
               if ( chi2 > 0 && chi2 > worst_chi2 ) {
                  worst_chi2 = chi2 ;
               }

            } // syi

         } // sxi

         if ( worst_chi2 <= 0 ) worst_chi2 = best_chi2 ;

       //--- allow for a closer point to be the best if the difference in chi2 is small.
         float min_global_chi2_at_closest_point = get_th2_min( h_chi2, min_x, min_y ) ;

         best_chi2_dx = min_x ;
         best_chi2_dy = min_y ;

         float best_chi2_dr = sqrt( best_chi2_dx*best_chi2_dx + best_chi2_dy*best_chi2_dy ) ;
         best_chi2 = min_global_chi2_at_closest_point ;

         nhcal = hcal_active_path_length.size() ;
         necal = ecal_active_path_length.size() ;
         int nfcs = nhcal + necal ;

         std::map<int,float>* ecal_hit_energies;
         std::map<int,float>* hcal_hit_energies;
         if ( use_mc ) {
            ecal_hit_energies = &ecal_mc_hit_energies ;
            hcal_hit_energies = &hcal_mc_hit_energies ;
         } else {
            ecal_hit_energies = &ecal_rec_hit_energies ;
            hcal_hit_energies = &hcal_rec_hit_energies ;
         }

         float energy_on_track_ecal(0.) ;
         float energy_on_track_hcal(0.) ;
         for ( auto it = ecal_active_path_length.begin(); it != ecal_active_path_length.end(); it++ ) {
            int hit_id = it->first ;
            if ( ecal_hit_energies->find( hit_id ) != ecal_hit_energies->end() ) energy_on_track_ecal += (*ecal_hit_energies)[ hit_id ] ;
         }
         for ( auto it = hcal_active_path_length.begin(); it != hcal_active_path_length.end(); it++ ) {
            int hit_id = it->first ;
            if ( hcal_hit_energies->find( hit_id ) != hcal_hit_energies->end() ) energy_on_track_hcal += (*hcal_hit_energies)[ hit_id ] ;
         }


         if ( arg_verbose ) verbose = true ;

         if ( verbose ) printf("\n\n Rerunning at best chi2 point.   dx = %5.2f  dy = %5.2f\n", best_chi2_dx, best_chi2_dy ) ;
         do_dz_steps( ecal_trkproj_fcsx + best_chi2_dx, ecal_trkproj_fcsy + best_chi2_dy, ecal_trkproj_fcsz, trk_dxdz, trk_dydz,
                      iso1_chi2_ecal, iso1_chi2_hcal,
                      iso2_chi2_ecal, iso2_chi2_hcal,
                      energy_chi2_ecal, energy_chi2_hcal, use_mc ) ;

         iso1_chi2 = iso1_chi2_ecal + iso1_chi2_hcal ;
         iso2_chi2 = iso2_chi2_ecal + iso2_chi2_hcal ;
         energy_chi2 = energy_chi2_ecal + energy_chi2_hcal ;

       //-- check for insufficient energy in either and set it for both if found.
         if ( iso1_chi2_ecal < 0 ) iso1_chi2 = iso1_chi2_ecal ;
         if ( iso1_chi2_hcal < 0 ) iso1_chi2 = iso1_chi2_hcal ;
         if ( iso2_chi2_ecal < 0 ) iso2_chi2 = iso2_chi2_ecal ;
         if ( iso2_chi2_hcal < 0 ) iso2_chi2 = iso2_chi2_hcal ;
         if ( energy_chi2_ecal < 0 ) energy_chi2 = energy_chi2_ecal ;
         if ( energy_chi2_hcal < 0 ) energy_chi2 = energy_chi2_hcal ;

         if ( iso1_chi2_ecal == 0 && energy_chi2_ecal > 0 && energy_chi2_ecal < 1e7 ) iso1_chi2_ecal = 0.01 ;
         if ( iso1_chi2_hcal == 0 && energy_chi2_hcal > 0 && energy_chi2_hcal < 1e7 ) iso1_chi2_hcal = 0.01 ;
         if ( iso2_chi2_ecal == 0 && energy_chi2_ecal > 0 && energy_chi2_ecal < 1e7 ) iso2_chi2_ecal = 0.01 ;
         if ( iso2_chi2_hcal == 0 && energy_chi2_hcal > 0 && energy_chi2_hcal < 1e7 ) iso2_chi2_hcal = 0.01 ;

         if ( iso1_chi2 == 0 && energy_chi2 > 0 && energy_chi2 < 1e7 ) iso1_chi2 = 0.01 ; // for log10 plots
         if ( iso2_chi2 == 0 && energy_chi2 > 0 && energy_chi2 < 1e7 ) iso2_chi2 = 0.01 ; // for log10 plots

         float gradx, grady ;
         get_gradients( h_chi2, best_chi2_dx, best_chi2_dy, gradx, grady ) ;
         float best_dist(99.) ;

         if ( iso1_chi2 > 0 && iso1_chi2 < 1e6 && energy_chi2 > 0 && energy_chi2 < 1e6 ) {
            if ( gradx > grady ) {
               best_dist = best_chi2_dx ;
               if ( gradx > 50 ) {
                  h_dist_from_trk_all -> Fill( best_chi2_dx ) ;
                  if ( mc_muon_match ) {
                     h_dist_from_trk_mu -> Fill( best_chi2_dx ) ;
                  } else {
                     h_dist_from_trk_notmu -> Fill( best_chi2_dx ) ;
                  }
               }
            } else {
               best_dist = best_chi2_dy ;
               if ( grady > 50 ) {
                  h_dist_from_trk_all -> Fill( best_chi2_dy ) ;
                  if ( mc_muon_match ) {
                     h_dist_from_trk_mu -> Fill( best_chi2_dy ) ;
                  } else {
                     h_dist_from_trk_notmu -> Fill( best_chi2_dy ) ;
                  }
               }
            }
         }


         if (  energy_chi2 > 0 ) {
            printf(" event %6llu track %3d : necal = %d, nhcal = %d  Best chi2 = %9.2f ( %9.2f, %9.2f),  worst chi2 = %9.2f,   best/worst = %.5f,  dist=%7.3f",
               jentry, rci, necal, nhcal, best_chi2, energy_chi2, iso1_chi2, worst_chi2, best_chi2/worst_chi2, best_dist ) ;
            if ( mc_muon_match ) { printf(" ** muon") ; }
            printf("\n") ;
         }


         if ( verbose ) printf("\n\n\n") ;

         h_log10chi2_all -> Fill( log10(iso1_chi2 + energy_chi2) ) ;
         h_log10chi2_energy_all -> Fill( log10(energy_chi2) ) ;
         h_log10chi2_iso_all -> Fill( log10(iso1_chi2) ) ;
         h_log10chi2_energy_vs_iso_all -> Fill( log10(iso1_chi2), log10(energy_chi2) ) ;

         if ( mc_muon_match ) {
            h_log10chi2_mu -> Fill( log10(iso1_chi2 + energy_chi2) ) ;
            h_log10chi2_energy_mu -> Fill( log10(energy_chi2) ) ;
            h_log10chi2_iso_mu -> Fill( log10(iso1_chi2) ) ;
            h_log10chi2_energy_vs_iso_mu -> Fill( log10(iso1_chi2), log10(energy_chi2) ) ;
         } else {
            h_log10chi2_notmu -> Fill( log10(iso1_chi2 + energy_chi2) ) ;
            h_log10chi2_energy_notmu -> Fill( log10(energy_chi2) ) ;
            h_log10chi2_iso_notmu -> Fill( log10(iso1_chi2) ) ;
            h_log10chi2_energy_vs_iso_notmu -> Fill( log10(iso1_chi2), log10(energy_chi2) ) ;
         }

       //--- fill output branch values.
         rc_chi2_energy_ecal.push_back( energy_chi2_ecal ) ;
         rc_chi2_energy_hcal.push_back( energy_chi2_hcal ) ;
         rc_chi2_energy.push_back( energy_chi2 ) ;
         rc_chi2_iso1_ecal.push_back( iso1_chi2_ecal ) ;
         rc_chi2_iso1_hcal.push_back( iso1_chi2_hcal ) ;
         rc_chi2_iso1.push_back( iso1_chi2 ) ;
         rc_chi2_iso2_ecal.push_back( iso2_chi2_ecal ) ;
         rc_chi2_iso2_hcal.push_back( iso2_chi2_hcal ) ;
         rc_chi2_iso2.push_back( iso2_chi2 ) ;
         rc_nhit_ecal.push_back( necal ) ;
         rc_nhit_hcal.push_back( nhcal ) ;
         rc_dx.push_back( best_chi2_dx ) ;
         rc_dy.push_back( best_chi2_dy ) ;
         rc_dx_gradient.push_back( gradx ) ;
         rc_dy_gradient.push_back( grady ) ;
         rc_energy_ecal.push_back( energy_on_track_ecal ) ;
         rc_energy_hcal.push_back( energy_on_track_hcal ) ;
         int mcPidval(0) ;
         if ( mcind >= 0 ) mcPidval = mcPid->at(mcind) ;
         rc_mcPid.push_back( mcPidval ) ;
         rc_rci.push_back( rci ) ;
         bool ismu(false) ;
         if ( mcPidval==5 || mcPidval==6 ) ismu = true ;
         rc_ismu.push_back( ismu ) ;
         rc_Pt.push_back( rcPt->at(rci) ) ;
         rc_Eta.push_back( rcEta->at(rci) ) ;

         if ( jentry%10 == 0 ) {

            int can_ind(1) ;

            can1 -> cd(can_ind++) ;
            h_log10chi2_energy_mu -> Draw() ;
            h_log10chi2_energy_notmu -> Draw("same") ;
            h_log10chi2_energy_all -> Draw("same") ;
            gPad -> SetGridx(1) ;

            can1 -> cd(can_ind++) ;
            h_log10chi2_iso_mu -> Draw() ;
            h_log10chi2_iso_notmu -> Draw("same") ;
            h_log10chi2_iso_all -> Draw("same") ;
            gPad -> SetGridx(1) ;


            can1 -> cd(can_ind++) ;
            h_dist_from_trk_mu -> Draw() ;
            gPad -> SetGridx(1) ;




            can1 -> cd(can_ind++) ;
            h_log10chi2_energy_all -> Draw() ;
            h_log10chi2_energy_mu -> Draw("same") ;
            h_log10chi2_energy_notmu -> Draw("same") ;
            gPad -> SetGridx(1) ;

            can1 -> cd(can_ind++) ;
            h_log10chi2_iso_all -> Draw() ;
            h_log10chi2_iso_mu -> Draw("same") ;
            h_log10chi2_iso_notmu -> Draw("same") ;
            gPad -> SetGridx(1) ;

            can1 -> cd(can_ind++) ;
            h_dist_from_trk_notmu -> Draw() ;
            gPad -> SetGridx(1) ;





            can1 -> cd(can_ind++) ;
            h_log10chi2_energy_vs_iso_mu -> Draw("colz") ;
            can1 -> cd(can_ind++) ;
            h_log10chi2_energy_vs_iso_notmu -> Draw("colz") ;

            can1 -> Draw() ;
            can1 -> Update() ;
            gSystem -> ProcessEvents() ;

         }



      } // rci

      if ( verbose ) {
         char answ = getchar() ;
         if ( answ == 'q' ) return ;
      }

      tt_out -> Fill() ;

      rc_chi2_energy_ecal.clear() ;
      rc_chi2_energy_hcal.clear() ;
      rc_chi2_energy.clear() ;
      rc_chi2_iso1_ecal.clear() ;
      rc_chi2_iso1_hcal.clear() ;
      rc_chi2_iso1.clear() ;
      rc_chi2_iso2_ecal.clear() ;
      rc_chi2_iso2_hcal.clear() ;
      rc_chi2_iso2.clear() ;
      rc_dx.clear() ;
      rc_dy.clear() ;
      rc_dx_gradient.clear() ;
      rc_dy_gradient.clear() ;
      rc_energy_ecal.clear() ;
      rc_energy_hcal.clear() ;
      rc_Pt.clear() ;
      rc_Eta.clear() ;
      rc_nhit_ecal.clear() ;
      rc_nhit_hcal.clear() ;
      rc_mcPid.clear() ;
      rc_rci.clear() ;
      rc_ismu.clear() ;


   } // jentry

   tt_out -> AutoSave() ;

   tf_output -> Close() ;

   saveHist("hists.root","h*" ) ;

} // Loop




