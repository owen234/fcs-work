#define analysis3_cxx
#include "analysis3.h"
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
//float hcal_reso2 = hcal_rec_over_mc_scale * hcal_rec_over_mc_scale * 0.002  * 0.002  ; // = (0.0020 * 68.3 )^2 = 0.137 ^2
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
                  float& iso_chi2, float& energy_chi2, bool use_mc=true  ) {

   iso_chi2 = 1e9 ;
   energy_chi2 = 1e9 ;

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

   std::map<int,float> ecal_neighbor_energies ;
   std::map<int,float> hcal_neighbor_energies ;

   for ( auto it = ecal_active_path_length.begin(); it != ecal_active_path_length.end(); it++ ) {
      int hit_id = it->first ;
      int hit_row, hit_col ;
      ecal_id_to_row_col( hit_id, hit_row, hit_col ) ;
      for ( int neighbor_row = hit_row-1; neighbor_row <= hit_row+1; neighbor_row++ ) {
         for ( int neighbor_col = hit_col-1; neighbor_col <= hit_col+1; neighbor_col++ ) {
            if ( neighbor_row == hit_row && neighbor_col == hit_col ) continue ; // don't count the hit!
            int neighbor_id = ecal_row_col_to_id( neighbor_row, neighbor_col ) ;
            if ( ecal_active_path_length.find( neighbor_id ) != ecal_active_path_length.end() ) continue ; // don't include any other hits with length>0.
            if ( ecal_hit_energies->find( neighbor_id ) != ecal_hit_energies->end() ) {
               float energy = (*ecal_hit_energies)[ neighbor_id ] ;
               if ( use_mc ) energy = ecal_rec_over_mc_scale * energy ;
               ecal_neighbor_energies[ neighbor_id ] = energy ;
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
      for ( int neighbor_row = hit_row-1; neighbor_row <= hit_row+1; neighbor_row++ ) {
         for ( int neighbor_col = hit_col-1; neighbor_col <= hit_col+1; neighbor_col++ ) {
            if ( neighbor_row == hit_row && neighbor_col == hit_col ) continue ; // don't count the hit!
            int neighbor_id = hcal_row_col_to_id( neighbor_row, neighbor_col ) ;
            if ( hcal_active_path_length.find( neighbor_id ) != hcal_active_path_length.end() ) continue ; // don't include any other hits with length>0.
            if ( hcal_hit_energies->find( neighbor_id ) != hcal_hit_energies->end() ) {
               float energy = (*hcal_hit_energies)[ neighbor_id ] ;
               if ( use_mc ) energy = hcal_rec_over_mc_scale * energy ;
               hcal_neighbor_energies[ neighbor_id ] = energy ;
               if ( verbose ) printf("   HCAL isolation:  hit id %5d (%2d,%2d)    neighbor %5d (%2d,%2d)    energy = %7.3f\n",
                   hit_id, hit_row, hit_col, neighbor_id, neighbor_row, neighbor_col, energy ) ;
            }
         }
      }
   }

   float ecal_iso_e2_sum(0.) ;
   float hcal_iso_e2_sum(0.) ;
   for ( auto it = ecal_neighbor_energies.begin(); it != ecal_neighbor_energies.end(); it++ ) {
      float energy = it->second ;
      if ( verbose ) printf("   calc ECAL iso sum:  id = %5d,  E = %7.4f\n", it->first, energy ) ;
      ecal_iso_e2_sum += energy*energy ;
   }
   for ( auto it = hcal_neighbor_energies.begin(); it != hcal_neighbor_energies.end(); it++ ) {
      float energy = it->second ;
      if ( verbose ) printf("   calc HCAL iso sum:  id = %5d,  E = %7.4f\n", it->first, energy ) ;
      hcal_iso_e2_sum += energy*energy ;
   }


   iso_chi2 = ecal_iso_e2_sum / ecal_reso2 + hcal_iso_e2_sum / hcal_reso2 ;
   if ( verbose ) printf("   Isolation chi2 = %7.2f\n", iso_chi2 ) ;




   std::set<int> ecal_ids ;
   std::set<int> hcal_ids ;

   for ( auto it = ecal_active_path_length.begin(); it != ecal_active_path_length.end(); it++ ) { ecal_ids.insert( it->first ) ; }
   for ( auto it = hcal_active_path_length.begin(); it != hcal_active_path_length.end(); it++ ) { hcal_ids.insert( it->first ) ; }






   energy_chi2 = 0. ;

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
         energy_chi2 += diff2 / ecal_reso2 ;
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
         energy_chi2 += diff2 / hcal_reso2 ;
         hcal_on_track_esum += energy ;
         if ( verbose ) printf("  HCAL :  id = %5d   L = %6.1f  E = %7.3f  pred E = %7.3f    diff2/reso2 = %7.3f\n", *it, length, energy, pred_e, diff2 / hcal_reso2 ) ;
      }

   }

   float global_chi2 = iso_chi2 + energy_chi2 ;

   //-- set it to a negative number of not enough energy to be valid.
   if ( ecal_on_track_esum < 0.05 ) { energy_chi2 = -10000. ; iso_chi2 = -10000. ; }
   if ( hcal_on_track_esum < 0.10 ) { energy_chi2 = -10000. ; iso_chi2 = -10000. ; }

   if ( verbose ) printf("  chi2 =  energy %.2f + iso %.2f = %.2f\n", energy_chi2, iso_chi2, global_chi2 ) ;
   if ( verbose ) printf("\n") ;



} // do_dz_steps
//---------------

void analysis3::load_fcs_mc_hit_energy_maps() {
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

void analysis3::load_fcs_rec_hit_energy_maps() {
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
void analysis3::Loop( bool use_mc, bool quiet )
{

   gStyle -> SetPadRightMargin(0.11) ;

   bool dump = true ;

   bool do_display = true ;
   bool wait = true ;
   verbose = true ;

   if ( quiet ) {
      do_display = false ;
      wait = false ;
      verbose = false ;
      dump = false ;
   }

   float scan_dx = 0.2 ;
   int   scan_npoints = 51 ;

   //float scan_dx = 0.4 ;
   //int   scan_npoints = 61 ;

   float scan_dy = scan_dx ;
   float scan_xmax = scan_dx * (scan_npoints-1)/2 ;
   float scan_ymax = scan_xmax ;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


   TH1F* h_theta = new TH1F( "h_theta", "theta", 80, 0., 0.01 ) ;

   TH3F* h_dr = new TH3F( "h_dr", "dr", 50, -100., 100., 50, -100., 100., 50, -30., 30. ) ;
   TH3F* h_drdir = new TH3F( "h_drdir", "dr", 50, -100., 100., 50, -100., 100., 50, -3.2, 3.2 ) ;


  //--- display TH2s

   TH2F* h_display_xy = new TH2F( "h_display_xy", "Display, xy", 800, -100., 100., 800, -100., 100. ) ;
   TH2F* h_display_rz = new TH2F( "h_display_rz", "Display, rz", 1000, 0., 1000., 800, 0., 100. ) ;

   TH2F* h_chi2 = new TH2F( "h_chi2", "chi2", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx, scan_npoints, -1*scan_ymax-0.5*scan_dx, scan_ymax+0.5*scan_dx ) ;
   TH2F* h_chi2_iso = new TH2F( "h_chi2_iso", "iso chi2", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx, scan_npoints, -1*scan_ymax-0.5*scan_dx, scan_ymax+0.5*scan_dx ) ;
   TH2F* h_chi2_energy = new TH2F( "h_chi2_energy", "energy chi2", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx, scan_npoints, -1*scan_ymax-0.5*scan_dx, scan_ymax+0.5*scan_dx ) ;

   TH2F* h_fcs_x_pred_vs_calc = new TH2F("h_fcs_x_pred_vs_calc", "FCS x, pred vs calc", 200, 0., 150, 200, 0., 150. ) ;
   TH2F* h_fcs_y_pred_vs_calc = new TH2F("h_fcs_y_pred_vs_calc", "FCS y, pred vs calc", 200, -150, 150, 200, -150, 150. ) ;

   TH1F* h_hits_energy_ecal_obs = new TH1F( "h_hits_energy_ecal_obs", "", 12, 0.5, 12.5 ) ;
   TH1F* h_hits_energy_ecal_pred = new TH1F( "h_hits_energy_ecal_pred", "", 12, 0.5, 12.5 ) ;

   TH1F* h_hits_energy_hcal_obs = new TH1F( "h_hits_energy_hcal_obs", "", 12, 0.5, 12.5 ) ;
   TH1F* h_hits_energy_hcal_pred = new TH1F( "h_hits_energy_hcal_pred", "", 12, 0.5, 12.5 ) ;

   h_hits_energy_ecal_obs -> SetLineWidth(3) ;
   h_hits_energy_ecal_pred -> SetLineWidth(3) ;
   h_hits_energy_hcal_obs -> SetLineWidth(3) ;
   h_hits_energy_hcal_pred -> SetLineWidth(3) ;

   h_hits_energy_ecal_obs -> SetLineColor(4) ;
   h_hits_energy_ecal_pred -> SetLineColor(4) ;
   h_hits_energy_hcal_obs -> SetLineColor(2) ;
   h_hits_energy_hcal_pred -> SetLineColor(2) ;

   h_hits_energy_ecal_pred -> SetLineStyle(2) ;
   h_hits_energy_hcal_pred -> SetLineStyle(2) ;

   h_hits_energy_ecal_obs -> SetMarkerStyle(20) ;
   h_hits_energy_hcal_obs -> SetMarkerStyle(20) ;

   h_hits_energy_ecal_obs -> SetMarkerColor(4) ;
   h_hits_energy_hcal_obs -> SetMarkerColor(2) ;




  //--- analysis hists

   TH1F* h_best_log10chi2 = new TH1F( "h_best_log10chi2", "best log10 chi2", 80, -3., 8. ) ;
   TH1F* h_best_log10chi2_iso_comp = new TH1F( "h_best_log10chi2_iso_comp", "best log10 chi2, isolation component", 80, -3., 8. ) ;
   TH1F* h_best_log10chi2_energy_comp = new TH1F( "h_best_log10chi2_energy_comp", "best log10 chi2, energy component", 80, -3., 8. ) ;
   TH1F* h_best_chi2_dist = new TH1F( "h_best_chi2_dist", "dist. from trk proj at best chi2", 60, 0., 15. ) ;
   TH2F* h_best_chi2_xy = new TH2F( "h_best_chi2_xy", "relative position from trk proj at best chi2", scan_npoints, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx, scan_npoints, -1*scan_ymax-0.5*scan_dx, scan_ymax+0.5*scan_dx ) ;
   TH2F* h_best_vs_worst_log10chi2 = new TH2F( "h_best_vs_worst_log10chi2", "best vs worst log10 chi2", 80, -3., 8., 80, -3., 8. ) ;
   TH1F* h_best_over_worst_chi2 = new TH1F( "h_best_over_worst_chi2", "best/worst chi2", 60, 0., 1. ) ;
   TH1F* h_log10_best_over_worst_chi2 = new TH1F( "h_log10_best_over_worst_chi2", "log10 best/worst chi2", 60, -7., 1. ) ;

   TH1F* h_ncell_hcal = new TH1F( "h_ncell_hcal", "Number of HCAL cells", 8, -0.5, 7.5 ) ;
   TH1F* h_ncell_ecal = new TH1F( "h_ncell_ecal", "Number of ECAL cells", 8, -0.5, 7.5 ) ;
   TH1F* h_ncell_fcs  = new TH1F( "h_ncell_fcs", "Number of ECAL+HCAL cells", 8, -0.5, 7.5 ) ;
   TH2F* h_ncell_hcal_vs_ecal = new TH2F( "h_ncell_hcal_vs_ecal", "Number of HCAL vs ECAL cells", 8, -0.5, 7.5, 8, -0.5, 7.5 ) ;

   TH2F* h_ntracks_in_star_xy = new TH2F( "h_ntracks_in_star_xy", "Number of tracks in star xy plane", 60, -150., 150., 60, -150., 150. ) ;

   TH3F* h_ncell_ecal_in_star_xy = new TH3F( "h_ncell_ecal_in_star_xy", "Number of ECAL cells in star xy plane", 60, -150., 150., 60, -150., 150., 5, -0.5, 4.5 ) ;
   TH3F* h_ncell_hcal_in_star_xy = new TH3F( "h_ncell_hcal_in_star_xy", "Number of HCAL cells in star xy plane", 60, -150., 150., 60, -150., 150., 5, -0.5, 4.5 ) ;
   TH3F* h_ncell_fcs_in_star_xy = new TH3F( "h_ncell_fcs_in_star_xy", "Number of ECAL+HCAL cells in star xy plane", 60, -150., 150., 60, -150., 150., 8, -0.5, 7.5 ) ;

   TH1F* h_cutflow = new TH1F( "h_cutflow","Cutflow", 10, 0.5, 10.5 ) ;
   TH1F* h_cutflow2 = new TH1F( "h_cutflow2","Cutflow2", 10, 0.5, 10.5 ) ;


   TH1F* h_best_log10chi2_2cell = new TH1F( "h_best_log10chi2_2cell", "best log10 chi2, 2 FCS cells", 80, -3., 8. ) ;
   TH1F* h_best_log10chi2_3cell = new TH1F( "h_best_log10chi2_3cell", "best log10 chi2, 3 FCS cells", 80, -3., 8. ) ;
   TH1F* h_best_log10chi2_ge4cell = new TH1F( "h_best_log10chi2_ge4cell", "best log10 chi2, >=4 FCS cells", 80, -3., 8. ) ;




   TH1F* h_ecal_mc_track_sum = new TH1F( "h_ecal_mc_track_sum", "ECAL hits on track sum", 100, 0., 0.15 ) ;
   TH1F* h_hcal_mc_track_sum = new TH1F( "h_hcal_mc_track_sum", "HCAL hits on track sum", 100, 0., 0.06 ) ;

   TH1F* h_ecal_mc_track_sum_chi2cut = new TH1F( "h_ecal_mc_track_sum_chi2cut", "ECAL hits on track sum, chi2 cut", 100, 0., 0.15 ) ;
   TH1F* h_hcal_mc_track_sum_chi2cut = new TH1F( "h_hcal_mc_track_sum_chi2cut", "HCAL hits on track sum, chi2 cut", 100, 0., 0.06 ) ;

   TH1F* h_ecal_mc_event_sum = new TH1F( "h_ecal_event_sum", "ECAL event sum",         100, 0., 2. ) ;
   TH1F* h_hcal_mc_event_sum = new TH1F( "h_hcal_event_sum", "HCAL event sum",         100, 0., 0.5 ) ;

   TH1F* h_ecal_mc_event_sum_chi2cut = new TH1F( "h_ecal_mc_event_sum_chi2cut", "ECAL event sum, chi2 cut",         100, 0., 2. ) ;
   TH1F* h_hcal_mc_event_sum_chi2cut = new TH1F( "h_hcal_mc_event_sum_chi2cut", "HCAL event sum, chi2 cut",         100, 0., 0.5 ) ;




   TH1F* h_ecal_rec_track_sum = new TH1F( "h_ecal_rec_track_sum", "ECAL hits on track sum", 100, 0., 0.15*ecal_rec_over_mc_scale ) ;
   TH1F* h_hcal_rec_track_sum = new TH1F( "h_hcal_rec_track_sum", "HCAL hits on track sum", 100, 0., 0.06*hcal_rec_over_mc_scale ) ;

   TH1F* h_ecal_rec_track_sum_chi2cut = new TH1F( "h_ecal_rec_track_sum_chi2cut", "ECAL hits on track sum, chi2 cut", 100, 0., 0.15*ecal_rec_over_mc_scale ) ;
   TH1F* h_hcal_rec_track_sum_chi2cut = new TH1F( "h_hcal_rec_track_sum_chi2cut", "HCAL hits on track sum, chi2 cut", 100, 0., 0.06*hcal_rec_over_mc_scale ) ;

   TH1F* h_ecal_rec_event_sum = new TH1F( "h_ecal_event_sum", "ECAL event sum",         100, 0., 2.*ecal_rec_over_mc_scale ) ;
   TH1F* h_hcal_rec_event_sum = new TH1F( "h_hcal_event_sum", "HCAL event sum",         100, 0., 0.5*hcal_rec_over_mc_scale ) ;

   TH1F* h_ecal_rec_event_sum_chi2cut = new TH1F( "h_ecal_rec_event_sum_chi2cut", "ECAL event sum, chi2 cut",         100, 0., 2.*ecal_rec_over_mc_scale ) ;
   TH1F* h_hcal_rec_event_sum_chi2cut = new TH1F( "h_hcal_rec_event_sum_chi2cut", "HCAL event sum, chi2 cut",         100, 0., 0.5*hcal_rec_over_mc_scale ) ;






   TH1F* h_ecal_track_nhits = new TH1F( "h_ecal_track_nhits", "ECAL N hits on track", 7, -0.5, 6.5 ) ;
   TH1F* h_ecal_event_nhits = new TH1F( "h_ecal_event_nhits", "ECAL N hits in event", 201, -0.5, 200.5 ) ;

   TH1F* h_hcal_track_nhits = new TH1F( "h_hcal_track_nhits", "HCAL N hits on track", 7, -0.5, 6.5 ) ;
   TH1F* h_hcal_event_nhits = new TH1F( "h_hcal_event_nhits", "HCAL N hits in event", 201, -0.5, 200.5 ) ;

   TH1F* h_ecal_track_nhits_chi2cut = new TH1F( "h_ecal_track_nhits_chi2cut", "ECAL N hits on track, chi2 cut", 7, -0.5, 6.5 ) ;
   TH1F* h_ecal_event_nhits_chi2cut = new TH1F( "h_ecal_event_nhits_chi2cut", "ECAL N hits in event, chi2 cut", 201, -0.5, 200.5 ) ;

   TH1F* h_hcal_track_nhits_chi2cut = new TH1F( "h_hcal_track_nhits_chi2cut", "HCAL N hits on track, chi2 cut", 7, -0.5, 6.5 ) ;
   TH1F* h_hcal_event_nhits_chi2cut = new TH1F( "h_hcal_event_nhits_chi2cut", "HCAL N hits in event, chi2 cut", 201, -0.5, 200.5 ) ;


   TH2F* h_ecal_event_sum_rec_vs_mc = new TH2F( "h_ecal_event_sum_rec_vs_mc", "ECAL event sum, rec vs mc",         100, 0., 2.  , 100, 0., 10. ) ;
   TH2F* h_hcal_event_sum_rec_vs_mc = new TH2F( "h_hcal_event_sum_rec_vs_mc", "HCAL event sum, rec vs mc",         100, 0., 0.5 , 100, 0., 40. ) ;

   TH1F* h_ecal_energy_rec_over_mc = new TH1F( "h_ecal_energy_rec_over_mc", "ECAL, rec energy over mc energy", 60, 4., 5.5 ) ;
   TH1F* h_hcal_energy_rec_over_mc = new TH1F( "h_hcal_energy_rec_over_mc", "HCAL, rec energy over mc energy", 60, 55., 75. ) ;


   TH2F* h_track_nhcal_vs_necal = new TH2F( "h_track_nhcal_vs_necal", "N hits on track, HCAL vs ECAL", 7, -0.5, 6.5, 7, -0.5, 6.5 ) ;
   TH2F* h_track_nhcal_vs_necal_chi2cut = new TH2F( "h_track_nhcal_vs_necal_chi2cut", "N hits on track, HCAL vs ECAL, chi2 cut", 7, -0.5, 6.5, 7, -0.5, 6.5 ) ;



   TH1F* h_chi2_grad_val = new TH1F( "h_chi2_grad_val", "chi2 gradient val", 60, 0., 600. ) ;
   TH1F* h_chi2_log10grad_val = new TH1F( "h_chi2_log10grad_val", "chi2 gradient val", 60, -1.5, 3. ) ;
   TH1F* h_chi2_grad_dist_from_trk_gradcut = new TH1F( "h_chi2_grad_dist_from_trk_gradcut", "chi2 gradient distance from track, grad cut", 60, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx  ) ;
   TH1F* h_chi2_grad_dist_from_trk_gradcut_x = new TH1F( "h_chi2_grad_dist_from_trk_gradcut_x", "chi2 gradient distance from track, grad cut, x", 60, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx  ) ;
   TH1F* h_chi2_grad_dist_from_trk_gradcut_y = new TH1F( "h_chi2_grad_dist_from_trk_gradcut_y", "chi2 gradient distance from track, grad cut, y", 60, -1*scan_xmax-0.5*scan_dx, scan_xmax+0.5*scan_dx  ) ;
   TH1F* h_chi2_grad_necal = new TH1F( "h_chi2_grad_necal", "chi2 gradient, n ECAL", 7, -0.5, 6.5 ) ;
   TH1F* h_chi2_grad_nhcal = new TH1F( "h_chi2_grad_nhcal", "chi2 gradient, n HCAL", 7, -0.5, 6.5 ) ;
   TH2F* h_chi2_grad_nhcal_vs_necal = new TH2F( "h_chi2_grad_nhcal_vs_necal", "chi2 gradient, n HCAL vs n ECAL", 7, -0.5, 6.5, 7, -0.5, 6.5 ) ;



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
   TCanvas* can_chi2(0x0) ;
   //if ( do_display ) {
      gStyle -> SetOptStat(0) ;
      can_xy = new TCanvas( "can_xy", "xy", 50, 50, 800, 800 ) ;
      can_xy -> Draw() ;
      can_rz = new TCanvas( "can_rz", "rz", 850, 50, 900, 700 ) ;
      can_rz -> Draw() ;
      can_chi2 = new TCanvas( "can_chi2", "chi2", 1750, 50, 800, 800 ) ;
      can_chi2 -> Draw() ;
      gSystem -> ProcessEvents() ;
   //}


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (fttN<=0) continue ;
      if ( mcN<=0) continue ;
      if ( rcN<=0) continue ;
      if ( tprojN<=0) continue ;
      if ( fcsN<=0 ) continue ;

      float mcpx(0.) ;
      float mcpy(0.) ;
      float mcpz(0.) ;
      float mcz(0.) ;
      for ( int i=0; i<fttN; i++ ) {
         if ( fttVolumeId->at(i) == 3 && fttZ->at(i) > 330. && fttMcpz->at(i) > 0.5 ) {
            mcpx = fttMcpx->at(i) ;
            mcpy = fttMcpy->at(i) ;
            mcpz = fttMcpz->at(i) ;
            mcz = fttZ->at(i) ;
         }
      }
      float tpx_ecal(0.) ;
      float tpy_ecal(0.) ;
      float tpz_ecal(0.) ;
      float tpx_hcal(0.) ;
      float tpy_hcal(0.) ;
      float tpz_hcal(0.) ;
      int tpi(-1) ;
      float mindz(99999.) ;
      for ( int i=0; i<tprojN; i++ ) {
         float dz = tprojZ->at(i) - mcz ;
         if ( fabs(dz) < fabs(mindz) ) {
            mindz = dz ;
            tpi = i ;
         }
         if ( tprojIdD->at(i) == 4 ) {
            tpx_ecal = tprojX->at(i) ;
            tpy_ecal = tprojY->at(i) ;
            tpz_ecal = tprojZ->at(i) ;
         }
         if ( tprojIdD->at(i) == 5 ) {
            tpx_hcal = tprojX->at(i) ;
            tpy_hcal = tprojY->at(i) ;
            tpz_hcal = tprojZ->at(i) ;
         }
      }
      float tppx(0.) ;
      float tppy(0.) ;
      float tppz(0.) ;
      if ( tpi >= 0 ) {
         tppx = tprojPx->at(tpi) ;
         tppy = tprojPy->at(tpi) ;
         tppz = tprojPz->at(tpi) ;
      }

      float fcsx(0.) ;
      float fcsy(0.) ;
      float fcsz(0.) ;
      float fcse(0.) ;
      if ( fcsN == 1 ) {
         fcsx = fcsX->at(0) ;
         fcsy = fcsY->at(0) ;
         fcsz = fcsZ->at(0) ;
         fcse = fcsE->at(0) ;
      }

      float mcP = sqrt( mcpx*mcpx + mcpy*mcpy + mcpz*mcpz ) ;
      float  tP = sqrt( tppx*tppx + tppy*tppy + tppz*tppz ) ;
      if ( mcP <=0 ) continue ;
      if (  tP <=0 ) continue ;
      float dot = mcpx*tppx + mcpy*tppy + mcpz*tppz ;
      float costheta = dot/(mcP*tP) ;
      float theta = acos(costheta) ;

      float dx = tpx_ecal - fcsx ;
      float dy = tpy_ecal - fcsy ;
      float dz = tpz_ecal - fcsz ;
      float dr = sqrt( dx*dx + dy*dy + dz*dz ) ;
      float drdir = atan2( dy, dx ) ;

      h_dr -> Fill( fcsx, fcsy, dr ) ;
      h_drdir -> Fill( fcsx, fcsy, drdir ) ;

      h_theta -> Fill( theta ) ;

      printf(" === Event %3llu\n", jentry ) ;

      if ( dump ) {
         printf("  fttN = %d, mcN = %d, rcN = %d, tprojN = %d, fcsN = %d, fcs_mc_ecalN = %d, fcs_mc_hcalN = %d\n",
          fttN, mcN, rcN, tprojN, fcsN, fcs_mc_ecalN, fcs_mc_hcalN ) ;
         for ( int i=0; i<fttN; i++ ) {
            printf("  %2d: ftt hit  z = %7.1f, VID = %d,  MC px,py,pz = %7.2f, %7.2f, %7.2f\n",
               i, fttZ->at(i), fttVolumeId->at(i), fttMcpx->at(i), fttMcpy->at(i), fttMcpz->at(i) ) ;
         }
         for ( int i=0; i<tprojN; i++ ) {
            printf("  %2d: proj  IdD = %d, IdT = %d, x,y,z = %7.2f, %7.2f, %7.2f,  px,py,pz = %7.2f, %7.2f, %7.2f\n",
              i, tprojIdD->at(i), tprojIdT->at(i), tprojX->at(i), tprojY->at(i), tprojZ->at(i),
              tprojPx->at(i), tprojPy->at(i), tprojPz->at(i) ) ;
         }
         for ( int i=0; i<fcsN; i++ ) {
            printf("  %2d: FCS cluster  x,y,z = %7.2f, %7.2f, %7.2f   E = %7.4f\n",
               i, fcsX->at(i), fcsY->at(i), fcsZ->at(i), fcsE->at(i) ) ;
         }
         for ( int i=0; i<fcs_mc_ecalN; i++ ) {
            printf("  %2d: FCS ecal  cluster  Vid = %5d  x,y,z = %7.2f, %7.2f, %7.2f   E = %7.4f\n",
               i, fcs_mc_ecalVid->at(i), fcs_mc_ecalX->at(i), fcs_mc_ecalY->at(i), fcs_mc_ecalZ->at(i), fcs_mc_ecalE->at(i) ) ;
         }
         for ( int i=0; i<fcs_mc_hcalN; i++ ) {
            printf("  %2d: FCS hcal cluster  Vid = %5d  x,y,z = %7.2f, %7.2f, %7.2f   E = %7.4f\n",
               i, fcs_mc_hcalVid->at(i), fcs_mc_hcalX->at(i), fcs_mc_hcalY->at(i), fcs_mc_hcalZ->at(i), fcs_mc_hcalE->at(i) ) ;
         }

         printf("  dz , FTT hit, proj %7.4f\n", mindz ) ;
         printf("  MC    px, py, pz = %7.2f, %7.2f, %7.2f\n", mcpx, mcpy, mcpz ) ;
         printf("  proj  px, py, pz = %7.2f, %7.2f, %7.2f\n", tppx, tppy, tppz ) ;
         printf("  mcP = %7.3f, tP = %7.3f,  dot = %7.5f,  cos(theta) = %7.5f, theta = %10.6f\n", mcP, tP, dot, costheta, theta ) ;
         printf("  proj|fcs :  x  %7.2f | %7.2f     y  %7.2f | %7.2f    z  %7.2f | %7.2f\n", 
            tpx_ecal, fcsx,  tpy_ecal, fcsy,  tpz_ecal, fcsz ) ;
         printf("  proj-fcs :  dx, dy, dz, dr, drdir   %7.2f, %7.2f, %7.2f, %7.2f   %7.3f\n", dx, dy, dz, dr, drdir ) ;
      }


      if ( do_display ) {

         can_xy -> cd() ;
         h_display_xy -> SetXTitle( "x position (cm)" ) ;
         h_display_xy -> SetYTitle( "y position (cm)" ) ;
         h_display_xy -> Draw() ;

         double x[1] ;
         double y[1] ;

         for ( int i=0; i<fstN; i++ ) {
            x[0] = fstX->at(i) ;
            y[0] = fstY->at(i) ;
            tpm_fst -> DrawPolyMarker( 1, x, y ) ;
         }
         for ( int i=0; i<fttN; i++ ) {
            x[0] = fttX->at(i) ;
            y[0] = fttY->at(i) ;
            tpm_ftt -> DrawPolyMarker( 1, x, y ) ;
         }
         for ( int i=0; i<fcsN; i++ ) {
            x[0] = fcsX->at(i) ;
            y[0] = fcsY->at(i) ;
            /////////tpm_fcsc -> DrawPolyMarker( 1, x, y ) ;
         }
         for ( int i=0; i<fcs_mc_ecalN; i++ ) {
            //x[0] = fcs_mc_ecalX->at(i) ;
            //y[0] = fcs_mc_ecalY->at(i) ;
            //tpm_fcs_mc_ecal -> DrawPolyMarker( 1, x, y ) ;
            float x, y ;
          //if ( use_mc ) {
               x = fcs_mc_ecalX->at(i) ;
               y = fcs_mc_ecalY->at(i) ;
          //} else {
          //   x = fcs_rec_ecalX->at(i) ;
          //   y = fcs_rec_ecalY->at(i) ;
          //}
            float energy = fcs_mc_ecalE->at(i) ;
            float alpha = energy / 0.08 ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_ecal -> SetFillColorAlpha( kBlue, alpha ) ;
            tb_ecal -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
            tb_ecal_line -> DrawBox( x - 5.572/2.,  y - 5.572/2.,  x + 5.572/2, y + 5.572/2 ) ;
         }
         for ( int i=0; i<fcs_mc_hcalN; i++ ) {
            //x[0] = fcs_mc_hcalX->at(i) ;
            //y[0] = fcs_mc_hcalY->at(i) ;
            //tpm_fcs_mc_hcal -> DrawPolyMarker( 1, x, y ) ;
            float x, y ;
          //if ( use_mc ) {
               x = fcs_mc_hcalX->at(i) ;
               y = fcs_mc_hcalY->at(i) ;
          //} else {
          //   x = fcs_rec_ecalX->at(i) ;
          //   y = fcs_rec_ecalY->at(i) ;
          //}
            float energy = fcs_mc_hcalE->at(i) ;
            float alpha = energy / 0.05 ;
            if ( alpha > 1 ) alpha = 1 ;
            if ( alpha < 0 ) alpha = 0 ;
            tb_hcal -> SetFillColorAlpha( kRed, alpha ) ;
            tb_hcal -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
            tb_hcal_line -> DrawBox( x - 9.99/2.,  y - 9.99/2.,  x + 9.99/2, y + 9.99/2 ) ;
         }

         x[0] = tpx_ecal ;
         y[0] = tpy_ecal ;
         tpm_tp -> DrawPolyMarker( 1, x, y ) ;

         x[0] = tpx_hcal ;
         y[0] = tpy_hcal ;
         tpm_tp -> DrawPolyMarker( 1, x, y ) ;

         tl_pseudotrack -> DrawLine( 0., 0., x[0], y[0] ) ;

         can_xy -> Update() ;
         can_xy -> Draw() ;
         gSystem -> ProcessEvents() ;




         can_rz -> cd() ;
         h_display_rz -> SetXTitle("z position (cm)") ;
         h_display_rz -> SetYTitle("rho position (cm)") ;
         h_display_rz -> Draw() ;

         double r[1] ;
         double z[1] ;

         for ( int i=0; i<fstN; i++ ) {
            r[0] = sqrt( fstX->at(i) * fstX->at(i) + fstY->at(i) * fstY->at(i) ) ;
            z[0] = fstZ->at(i) ;
            tpm_fst -> DrawPolyMarker( 1, z, r ) ;
         }
         for ( int i=0; i<fttN; i++ ) {
            r[0] = sqrt( fttX->at(i) * fttX->at(i) + fttY->at(i) * fttY->at(i) ) ;
            z[0] = fttZ->at(i) ;
            tpm_ftt -> DrawPolyMarker( 1, z, r ) ;
         }
         for ( int i=0; i<fcsN; i++ ) {
            r[0] = sqrt( fcsX->at(i) * fcsX->at(i) + fcsY->at(i) * fcsY->at(i) ) ;
            z[0] = fcsZ->at(i) ;
            /////////tpm_fcsc -> DrawPolyMarker( 1, z, r ) ;
         }
         for ( int i=0; i<fcs_mc_ecalN; i++ ) {
            //r[0] = sqrt( fcs_mc_ecalX->at(i) * fcs_mc_ecalX->at(i) + fcs_mc_ecalY->at(i) * fcs_mc_ecalY->at(i) ) ;
            //z[0] = fcs_mc_ecalZ->at(i) ;
            //tpm_fcs_mc_ecal -> DrawPolyMarker( 1, z, r ) ;
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
         for ( int i=0; i<fcs_mc_hcalN; i++ ) {
            //r[0] = sqrt( fcs_mc_hcalX->at(i) * fcs_mc_hcalX->at(i) + fcs_mc_hcalY->at(i) * fcs_mc_hcalY->at(i) ) ;
            //z[0] = fcs_mc_hcalZ->at(i) ;
            //tpm_fcs_mc_hcal -> DrawPolyMarker( 1, z, r ) ;
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

         r[0] = sqrt( tpx_ecal*tpx_ecal + tpy_ecal*tpy_ecal ) ;
         z[0] = tpz_ecal ;
         tpm_tp -> DrawPolyMarker( 1, z, r ) ;

         r[0] = sqrt( tpx_hcal*tpx_hcal + tpy_hcal*tpy_hcal ) ;
         z[0] = tpz_hcal ;
         tpm_tp -> DrawPolyMarker( 1, z, r ) ;

         tl_pseudotrack -> DrawLine( 0., 0., z[0], r[0] ) ;

         can_rz -> Update() ;
         can_rz -> Draw() ;
         gSystem -> ProcessEvents() ;

      } // do_display ?

     //----------

      if ( verbose ) {
         printf("\n\n --- FCS MC ECAL hits\n") ;
         for ( int i=0; i<fcs_mc_ecalN; i++ ) {
             float fcsx, fcsy, fcsz ;
             float starx, stary, starz ;
             starx = fcs_mc_ecalX->at(i) ;
             stary = fcs_mc_ecalY->at(i) ;
             starz = fcs_mc_ecalZ->at(i) ;
             int id = fcs_mc_ecalVid -> at(i) % 1000 - 1 ;
             int row, col ;
             ecal_id_to_row_col( id, row, col ) ;
             int id_check = ecal_row_col_to_id( row, col ) ;
             float pred_fcsx ;
             float pred_fcsy ;
             ecal_row_col_to_fcsxy( row, col, pred_fcsx, pred_fcsy) ;
             int check_col ;
             int check_row ;
             ecal_fcsxy_to_row_col( pred_fcsx, pred_fcsy, check_row, check_col ) ;

             star_to_fcs_ecal( starx, stary, starz, fcsx, fcsy, fcsz ) ;
             printf("  %3d : volume_id = %6d check (%6d) row=%2d, col=%2d check (%2d, %2d):   Star x,y,z = (%7.2f, %7.2f, %7.2f),   FCS x,y,z = (%7.2f, %7.2f, %7.2f)  pred FCS x,y = (%7.2f, %7.2f)\n",
                i, fcs_mc_ecalVid -> at(i), id_check, row, col, check_row, check_col, starx, stary, starz, fcsx, fcsy, fcsz, pred_fcsx, pred_fcsy ) ;
             h_fcs_x_pred_vs_calc -> Fill( fcsx, pred_fcsx ) ;
             h_fcs_y_pred_vs_calc -> Fill( fcsy, pred_fcsy ) ;

         }
         printf("\n\n --- FCS rec ECAL hits\n") ;
         for ( int i=0; i<fcs_rec_ecalN; i++ ) {
             printf(" id %4d :  Star x,y,z = (%7.2f, %7.2f, %7.2f),   FCS x,y = (%7.2f, %7.2f)  E = %7.3f\n",
               fcs_rec_ecalId->at(i), fcs_rec_ecalX->at(i), fcs_rec_ecalY->at(i), fcs_rec_ecalZ->at(i), fcs_rec_ecalLX->at(i), fcs_rec_ecalLY->at(i), fcs_rec_ecalE->at(i) ) ;
         }

         printf("\n\n --- FCS MC HCAL hits\n") ;
         for ( int i=0; i<fcs_mc_hcalN; i++ ) {
             float fcsx, fcsy, fcsz ;
             float starx, stary, starz ;
             starx = fcs_mc_hcalX->at(i) ;
             stary = fcs_mc_hcalY->at(i) ;
             starz = fcs_mc_hcalZ->at(i) ;
             int id = fcs_mc_hcalVid -> at(i) % 1000 - 1 ;
             int row, col ;
             hcal_id_to_row_col( id, row, col ) ;
             int id_check = hcal_row_col_to_id( row, col ) ;
             float pred_fcsx(0.) ;
             float pred_fcsy(0.) ;
             hcal_row_col_to_fcsxy( row, col, pred_fcsx, pred_fcsy) ;
             int check_col(0) ;
             int check_row(0) ;
             hcal_fcsxy_to_row_col( pred_fcsx, pred_fcsy, check_row, check_col ) ;

             star_to_fcs_ecal( starx, stary, starz, fcsx, fcsy, fcsz ) ;
             printf("  %3d : volume_id = %6d check (%6d) row=%2d, col=%2d check (%2d, %2d):   Star x,y,z = (%7.2f, %7.2f, %7.2f),   FCS x,y,z = (%7.2f, %7.2f, %7.2f)  pred FCS x,y = (%7.2f, %7.2f)\n",
                i, fcs_mc_hcalVid -> at(i), id_check, row, col, check_row, check_col, starx, stary, starz, fcsx, fcsy, fcsz, pred_fcsx, pred_fcsy ) ;

         }
         printf("\n\n --- FCS rec HCAL hits\n") ;
         for ( int i=0; i<fcs_rec_hcalN; i++ ) {
             printf(" id %4d :  Star x,y,z = (%7.2f, %7.2f, %7.2f),   FCS x,y = (%7.2f, %7.2f)  E = %7.3f\n",
               fcs_rec_hcalId->at(i), fcs_rec_hcalX->at(i), fcs_rec_hcalY->at(i), fcs_rec_hcalZ->at(i), fcs_rec_hcalLX->at(i), fcs_rec_hcalLY->at(i), fcs_rec_hcalE->at(i) ) ;
         }
      } // verbose ?

      float best_chi2(99999999.) ;
      float min_x, min_y ;
      int necal(0), nhcal(0) ;
      int necal_enz(0), nhcal_enz(0) ;

      if ( verbose ) printf("\n") ;
      if ( tprojN == 11) {

         load_fcs_mc_hit_energy_maps() ;
         load_fcs_rec_hit_energy_maps() ;

         h_ecal_event_sum_rec_vs_mc -> Fill( ecal_mc_event_sum_energy, ecal_rec_event_sum_energy ) ;
         h_hcal_event_sum_rec_vs_mc -> Fill( hcal_mc_event_sum_energy, hcal_rec_event_sum_energy ) ;

         if ( ecal_mc_event_sum_energy > 0 ) h_ecal_energy_rec_over_mc -> Fill( ecal_rec_event_sum_energy/ecal_mc_event_sum_energy ) ;
         if ( hcal_mc_event_sum_energy > 0 ) h_hcal_energy_rec_over_mc -> Fill( hcal_rec_event_sum_energy/hcal_mc_event_sum_energy ) ;

         float ecal_trkproj_starx = tprojX->at(9) ;
         float ecal_trkproj_stary = tprojY->at(9) ;
         float ecal_trkproj_starz = tprojZ->at(9) ;

         float ecal_trkproj_fcsx ;
         float ecal_trkproj_fcsy ;
         float ecal_trkproj_fcsz ;
         star_to_fcs_ecal( ecal_trkproj_starx, ecal_trkproj_stary, ecal_trkproj_starz, ecal_trkproj_fcsx, ecal_trkproj_fcsy, ecal_trkproj_fcsz ) ;

         int trkproj_ecal_row ;
         int trkproj_ecal_col ;
         ecal_fcsxy_to_row_col( ecal_trkproj_fcsx, ecal_trkproj_fcsy, trkproj_ecal_row, trkproj_ecal_col ) ;
         if ( verbose ) printf("  Track projection, star x,y,z (%7.2f, %7.2f, %7.2f),  FCS x,y,z (%7.2f, %7.2f, %7.2f)   FCS ecal row = %2d, col = %2d\n",
              ecal_trkproj_starx, ecal_trkproj_stary, ecal_trkproj_starz, ecal_trkproj_fcsx, ecal_trkproj_fcsy, ecal_trkproj_fcsz,  trkproj_ecal_row, trkproj_ecal_col ) ;


         float trk_star_px = tprojPx->at(9) ;
         float trk_star_py = tprojPy->at(9) ;
         float trk_star_pz = tprojPz->at(9) ;

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


         float iso_chi2, energy_chi2 ;
         do_dz_steps( ecal_trkproj_fcsx, ecal_trkproj_fcsy, ecal_trkproj_fcsz, trk_dxdz, trk_dydz, iso_chi2, energy_chi2, use_mc ) ;

         h_chi2 -> Reset() ;
         h_chi2_iso -> Reset() ;
         h_chi2_energy -> Reset() ;

         float best_chi2_dx(0.) ;
         float best_chi2_dy(0.) ;
         float worst_chi2(0.) ;
         for ( int sxi = 0; sxi<scan_npoints; sxi++ ) {

            float dx = -1 * scan_xmax + sxi*scan_dx ;
            float scan_x = ecal_trkproj_fcsx + dx ;

            for ( int syi = 0; syi<scan_npoints; syi++ ) {

               float dy = -1 * scan_ymax + syi*scan_dy ;
               float scan_y = ecal_trkproj_fcsy + dy ;

               do_dz_steps( scan_x, scan_y, ecal_trkproj_fcsz, trk_dxdz, trk_dydz, iso_chi2, energy_chi2, use_mc ) ;
               float chi2 = iso_chi2 + energy_chi2 ;
               h_chi2 -> Fill( dx, dy, chi2 ) ;
               h_chi2_iso -> Fill( dx, dy, iso_chi2 ) ;
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

       //--- allow for a closer point to be the best if the difference in chi2 is small.
         float min_global_chi2_at_closest_point = get_th2_min( h_chi2, min_x, min_y ) ;

         best_chi2_dx = min_x ;
         best_chi2_dy = min_y ;

         float best_chi2_dr = sqrt( best_chi2_dx*best_chi2_dx + best_chi2_dy*best_chi2_dy ) ;
         best_chi2 = min_global_chi2_at_closest_point ;

         nhcal = hcal_active_path_length.size() ;
         necal = ecal_active_path_length.size() ;
         int nfcs = nhcal + necal ;



         if ( verbose ) printf("\n\n Rerunning at best chi2 point.   dx = %5.2f  dy = %5.2f\n", best_chi2_dx, best_chi2_dy ) ;
         load_fcs_mc_hit_energy_maps() ;
         load_fcs_rec_hit_energy_maps() ;
         do_dz_steps( ecal_trkproj_fcsx + best_chi2_dx, ecal_trkproj_fcsy + best_chi2_dy, ecal_trkproj_fcsz, trk_dxdz, trk_dydz, iso_chi2, energy_chi2, use_mc ) ;

         if ( verbose ) printf("   Best chi2 = %.2f,  worst chi2 = %.2f,   best/worst = %.5f\n\n", best_chi2, worst_chi2, best_chi2/worst_chi2 ) ;




         float gradx, grady ;
         get_gradients( h_chi2, best_chi2_dx, best_chi2_dy, gradx, grady ) ;

         if ( gradx > grady ) {
            h_chi2_grad_val -> Fill( gradx ) ;
            h_chi2_log10grad_val -> Fill( log10(gradx) ) ;
            if ( gradx > 50 ) {
               h_chi2_grad_dist_from_trk_gradcut -> Fill( best_chi2_dx ) ;
               h_chi2_grad_dist_from_trk_gradcut_x -> Fill( best_chi2_dx ) ;
               h_chi2_grad_necal -> Fill( necal ) ;
               h_chi2_grad_nhcal -> Fill( nhcal ) ;
               h_chi2_grad_nhcal_vs_necal -> Fill( necal, nhcal ) ;
            }
         } else {
            h_chi2_grad_val -> Fill( grady ) ;
            h_chi2_log10grad_val -> Fill( log10(grady) ) ;
            if ( grady > 50 ) {
               h_chi2_grad_dist_from_trk_gradcut -> Fill( best_chi2_dy ) ;
               h_chi2_grad_dist_from_trk_gradcut_y -> Fill( best_chi2_dy ) ;
               h_chi2_grad_necal -> Fill( necal ) ;
               h_chi2_grad_nhcal -> Fill( nhcal ) ;
               h_chi2_grad_nhcal_vs_necal -> Fill( necal, nhcal ) ;
            }
         }




         h_hits_energy_ecal_obs -> Reset() ;
         h_hits_energy_ecal_pred -> Reset() ;

         h_hits_energy_hcal_obs -> Reset() ;
         h_hits_energy_hcal_pred -> Reset() ;

         int hit_index(2) ;
         float ecal_mc_energy_on_track(0.) ;
         float ecal_rec_energy_on_track(0.) ;
         for ( auto it = ecal_active_path_length.begin(); it != ecal_active_path_length.end(); it++ ) {
            int id = it->first ;
            if ( ecal_mc_hit_energies.find( id )  != ecal_mc_hit_energies.end()  ) ecal_mc_energy_on_track  += ecal_mc_hit_energies[id] ;
            if ( ecal_rec_hit_energies.find( id ) != ecal_rec_hit_energies.end() ) ecal_rec_energy_on_track += ecal_rec_hit_energies[id] ;
            if ( use_mc ) {
               if ( ecal_mc_hit_energies.find( id )  != ecal_mc_hit_energies.end()  ) { if ( ecal_mc_hit_energies[id] > 0.10 ) necal_enz ++ ; }
               h_hits_energy_ecal_obs -> SetBinContent( hit_index, ecal_mc_hit_energies[id] ) ;
               float pred_energy = ecal_active_path_length[id] * ecal_rec_over_mc_scale * ecal_dedl ;
               h_hits_energy_ecal_pred -> SetBinContent( hit_index, pred_energy ) ;
               hit_index ++ ;
            } else {
               if ( ecal_rec_hit_energies.find( id )  != ecal_rec_hit_energies.end()  ) { if ( ecal_rec_hit_energies[id] > 0.10 ) necal_enz ++ ; }
               h_hits_energy_ecal_obs -> SetBinContent( hit_index, ecal_rec_hit_energies[id] ) ;
               float pred_energy = ecal_active_path_length[id] * ecal_dedl ;
               h_hits_energy_ecal_pred -> SetBinContent( hit_index, pred_energy ) ;
               hit_index ++ ;
            }
         }
         float hcal_mc_energy_on_track(0.) ;
         float hcal_rec_energy_on_track(0.) ;
         hit_index = 7 ;
         for ( auto it = hcal_active_path_length.begin(); it != hcal_active_path_length.end(); it++ ) {
            int id = it->first ;
            if ( hcal_mc_hit_energies.find( id )  != hcal_mc_hit_energies.end()  ) hcal_mc_energy_on_track += hcal_mc_hit_energies[id] ;
            if ( hcal_rec_hit_energies.find( id ) != hcal_rec_hit_energies.end() ) hcal_rec_energy_on_track += hcal_rec_hit_energies[id] ;
            if ( use_mc ) {
               if ( hcal_mc_hit_energies.find( id )  != hcal_mc_hit_energies.end()  ) { if ( hcal_mc_hit_energies[id] > 0.10 ) nhcal_enz ++ ; }
               h_hits_energy_hcal_obs -> SetBinContent( hit_index, hcal_mc_hit_energies[id] ) ;
               float pred_energy = hcal_active_path_length[id] * hcal_rec_over_mc_scale * hcal_dedl ;
               h_hits_energy_hcal_pred -> SetBinContent( hit_index, pred_energy ) ;
               hit_index ++ ;
            } else {
               if ( hcal_rec_hit_energies.find( id ) != hcal_rec_hit_energies.end() ) { if ( hcal_rec_hit_energies[id] > 0.10 ) nhcal_enz ++ ; }
               h_hits_energy_hcal_obs -> SetBinContent( hit_index, hcal_rec_hit_energies[id] ) ;
               float pred_energy = hcal_active_path_length[id] * hcal_dedl ;
               h_hits_energy_hcal_pred -> SetBinContent( hit_index, pred_energy ) ;
               hit_index ++ ;
            }
         }

      // printf("  Sum energy on track:  ECAL %7.4f,   HCAL %7.4f\n", ecal_energy_on_track, hcal_energy_on_track ) ;
      // printf("  Sum event energy:     ECAL %7.4f,   HCAL %7.4f\n", ecal_mc_event_sum_energy, hcal_mc_event_sum_energy ) ;

      // for ( int i=0; i<fcsN; i++ ) {
      //    printf("  FCS cluster energy : %7.4f\n", fcsE->at(i) ) ;
      // }


         h_ecal_mc_track_sum -> Fill( ecal_mc_energy_on_track ) ;
         h_hcal_mc_track_sum -> Fill( hcal_mc_energy_on_track ) ;

         h_ecal_mc_event_sum -> Fill( ecal_mc_event_sum_energy ) ;
         h_hcal_mc_event_sum -> Fill( hcal_mc_event_sum_energy ) ;


         h_ecal_rec_track_sum -> Fill( ecal_rec_energy_on_track ) ;
         h_hcal_rec_track_sum -> Fill( hcal_rec_energy_on_track ) ;

         h_ecal_rec_event_sum -> Fill( ecal_rec_event_sum_energy ) ;
         h_hcal_rec_event_sum -> Fill( hcal_rec_event_sum_energy ) ;




         h_ecal_track_nhits -> Fill( necal ) ;
         h_hcal_track_nhits -> Fill( nhcal ) ;

         h_ecal_event_nhits -> Fill( fcs_mc_ecalN ) ;
         h_hcal_event_nhits -> Fill( fcs_mc_hcalN ) ;

         h_track_nhcal_vs_necal -> Fill( necal, nhcal ) ;

         if ( best_chi2 < 100 ) {

            h_ecal_mc_track_sum_chi2cut -> Fill( ecal_mc_energy_on_track ) ;
            h_hcal_mc_track_sum_chi2cut -> Fill( hcal_mc_energy_on_track ) ;

            h_ecal_mc_event_sum_chi2cut -> Fill( ecal_mc_event_sum_energy ) ;
            h_hcal_mc_event_sum_chi2cut -> Fill( hcal_mc_event_sum_energy ) ;


            h_ecal_rec_track_sum_chi2cut -> Fill( ecal_rec_energy_on_track ) ;
            h_hcal_rec_track_sum_chi2cut -> Fill( hcal_rec_energy_on_track ) ;

            h_ecal_rec_event_sum_chi2cut -> Fill( ecal_rec_event_sum_energy ) ;
            h_hcal_rec_event_sum_chi2cut -> Fill( hcal_rec_event_sum_energy ) ;




            h_ecal_track_nhits_chi2cut -> Fill( necal ) ;
            h_hcal_track_nhits_chi2cut -> Fill( nhcal ) ;

            h_ecal_event_nhits_chi2cut -> Fill( fcs_mc_ecalN ) ;
            h_hcal_event_nhits_chi2cut -> Fill( fcs_mc_hcalN ) ;

            h_track_nhcal_vs_necal_chi2cut -> Fill( necal, nhcal ) ;

         }


         h_best_log10chi2 -> Fill( log10( best_chi2 ) ) ;
         h_best_log10chi2_iso_comp -> Fill( log10( iso_chi2 ) ) ;
         h_best_log10chi2_energy_comp -> Fill( log10( energy_chi2 ) ) ;

         if ( nfcs == 2 ) h_best_log10chi2_2cell -> Fill( log10( best_chi2 ) ) ;
         if ( nfcs == 3 ) h_best_log10chi2_3cell -> Fill( log10( best_chi2 ) ) ;
         if ( nfcs >= 4 ) h_best_log10chi2_ge4cell -> Fill( log10( best_chi2 ) ) ;

         h_best_chi2_dist -> Fill( best_chi2_dr ) ;
         h_best_chi2_xy -> Fill( best_chi2_dx, best_chi2_dy ) ;
         h_best_vs_worst_log10chi2 -> Fill( log10(worst_chi2), log10(best_chi2) ) ;
         h_best_over_worst_chi2 -> Fill( best_chi2 / worst_chi2 ) ;
         h_log10_best_over_worst_chi2 -> Fill( log10(best_chi2 / worst_chi2) ) ;

         h_ncell_hcal -> Fill( nhcal ) ;
         h_ncell_ecal -> Fill( necal ) ;
         h_ncell_fcs -> Fill( necal+nhcal ) ;
         h_ncell_hcal_vs_ecal -> Fill( necal, nhcal ) ;

         h_ncell_ecal_in_star_xy -> Fill( ecal_trkproj_starx, ecal_trkproj_stary, necal ) ;
         h_ncell_hcal_in_star_xy -> Fill( ecal_trkproj_starx, ecal_trkproj_stary, nhcal ) ;
         h_ncell_fcs_in_star_xy -> Fill( ecal_trkproj_starx, ecal_trkproj_stary, necal+nhcal ) ;

         h_ntracks_in_star_xy -> Fill( ecal_trkproj_starx, ecal_trkproj_stary ) ;

         h_cutflow -> Fill( 1. ) ;
         if ( best_chi2 < 100 ) {
            h_cutflow -> Fill( 2. ) ;
            if ( best_chi2_dr < 5. ) {
               h_cutflow -> Fill( 3. ) ;
               if ( best_chi2 / worst_chi2 < 0.05 ) {
                  h_cutflow -> Fill( 4. ) ;
               }
            }
         }

         h_cutflow2 -> Fill( 1. ) ;
         if ( best_chi2 < 10 ) {
            h_cutflow2 -> Fill( 2. ) ;
            if ( best_chi2_dr < 3. ) {
               h_cutflow2 -> Fill( 3. ) ;
               if ( best_chi2 / worst_chi2 < 0.01 ) {
                  h_cutflow2 -> Fill( 4. ) ;
               }
            }
         }




         if ( do_display ) {

            can_chi2 -> cd() ;
            can_chi2 -> Clear() ;
            can_chi2 -> Divide(2,2) ;

            can_chi2 -> cd(1) ;
            h_chi2_energy -> SetMinimum(0.01) ;
            h_chi2_energy -> SetMaximum(5000.) ;
            h_chi2_energy -> Draw("colz") ;
            gPad -> SetLogz(1) ;
            gPad -> SetGridx(1) ;
            gPad -> SetGridy(1) ;

            can_chi2 -> cd(2) ;
            h_chi2_iso -> SetMinimum(0.01) ;
            h_chi2_iso -> SetMaximum(5000.) ;
            h_chi2_iso -> Draw("colz") ;
            gPad -> SetLogz(1) ;
            gPad -> SetGridx(1) ;
            gPad -> SetGridy(1) ;

            can_chi2 -> cd(3) ;
            h_chi2 -> SetXTitle("Position rel. to track proj. x (cm)" ) ;
            h_chi2 -> SetYTitle("Position rel. to track proj. y (cm)" ) ;
            h_chi2 -> SetMinimum(0.01) ;
            h_chi2 -> SetMaximum(5000.) ;
            h_chi2 -> Draw("colz") ;
            gPad -> SetLogz(1) ;
            gPad -> SetGridx(1) ;
            gPad -> SetGridy(1) ;

            double tpmx[1] ;
            double tpmy[1] ;
            tpmx[0] = min_x ;
            tpmy[0] = min_y ;
            tpm_min -> DrawPolyMarker( 1, tpmx, tpmy ) ;

            can_chi2 -> cd(4) ;
            h_hits_energy_ecal_obs -> SetYTitle("Hit Energy (GeV)") ;
            h_hits_energy_ecal_obs -> SetLabelOffset( 0.99, "x" ) ;
            h_hits_energy_ecal_obs -> SetMaximum(1.7) ;
            h_hits_energy_ecal_obs -> Draw() ;
            h_hits_energy_ecal_pred -> Draw("same") ;
            h_hits_energy_hcal_obs -> Draw("same") ;
            h_hits_energy_hcal_pred -> Draw("same") ;
            h_hits_energy_ecal_obs-> Draw("psame") ;
            h_hits_energy_hcal_obs -> Draw("psame") ;


            can_chi2 -> Update() ;
            can_chi2 -> Draw() ;
            gSystem -> ProcessEvents() ;

         } // do_display ?

      }


      if ( wait ) {
       //if (best_chi2 < 400 ) {
       //if ( necal_enz > 1 ) {
            char answ = getchar() ;
            if ( answ == 'q' ) return ;
       //}
       //}
      }

      if ( quiet && jentry%10 == 0  ) {
         can_chi2 -> Clear() ;
         can_chi2 -> cd() ;
         h_best_log10chi2 -> Draw() ;
         can_chi2 -> Update() ;
         can_chi2 -> Draw() ;
         gSystem -> ProcessEvents() ;
      }




   } // jentry


   saveHist("hists.root","h*" ) ;

}
