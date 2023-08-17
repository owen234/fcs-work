#define mip_search1_cxx
#include "mip_search1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "histio.c"

//---------------

void  mip_search1::doTrackLoop( std::map<const char*,TH1F*> hist_map, std::map<const char*,TH2F*> hist_map2d, bool use_previous_event_tracks = false ) {

      vector<float>* this_rcProjEcalx ;
      vector<float>* this_rcProjEcaly ;
      vector<float>* this_rcProjHcalx ;
      vector<float>* this_rcProjHcaly ;
      vector<int>* this_rcNumFST ;

      if ( use_previous_event_tracks ) {
         this_rcProjEcalx = previous_event_rcProjEcalx ;
         this_rcProjEcaly = previous_event_rcProjEcaly ;
         this_rcProjHcalx = previous_event_rcProjHcalx ;
         this_rcProjHcaly = previous_event_rcProjHcaly ;
         this_rcNumFST = previous_event_rcNumFST ;
      } else {
         this_rcProjEcalx = rcProjEcalx ;
         this_rcProjEcaly = rcProjEcaly ;
         this_rcProjHcalx = rcProjHcalx ;
         this_rcProjHcaly = rcProjHcaly ;
         this_rcNumFST = rcNumFST ;
      }

      if ( verbose ) printf("\n\n At top of doTrackLoop, use_previous_event_tracks = %s   N tracks = %lu\n",
          (use_previous_event_tracks?"true":"false"), this_rcProjEcalx->size() ) ;

   //--- loop over tracks

      for ( int ti=0; ti<this_rcProjEcalx->size(); ti++ ) {


      //-- track selection cuts here.

   ////  int nfst = this_rcNumFST->at(ti) ;
   ////  if ( nfst <= 0 ) continue ;






         float tpex = this_rcProjEcalx->at(ti) ;
         float tpey = this_rcProjEcaly->at(ti) ;

         float tphx = this_rcProjHcalx->at(ti) ;
         float tphy = this_rcProjHcaly->at(ti) ;

         if ( verbose ) printf( "  %s  track %2d :   (x,y) = (%7.2f, %7.2f)\n",
           (use_previous_event_tracks?"prev":"same"), ti, tpex, tpey ) ;





   //--- loop over ECAL hits

         int ecal_nhit(0) ;
         float ecal_sum_e(0.) ;
         float ecal_sum_ex(0.) ;
         float ecal_sum_ey(0.) ;

         for ( int hi=0; hi<fcs_rec_ecalE->size(); hi++ ) {

            float energy = fcs_rec_ecalE->at(hi) ;

            if ( energy < 0.01 ) continue ;

            float hx = fcs_rec_ecalX->at(hi) ;
            float hy = fcs_rec_ecalY->at(hi) ;

            float dx = hx - tpex ;
            float dy = hy - tpey ;
            float dr = sqrt( dx*dx + dy*dy ) ;

            if ( dr < radius ) {

               ecal_sum_e += energy ;
               ecal_sum_ex += hx * energy ;
               ecal_sum_ey += hy * energy ;

               ecal_nhit++ ;

               if ( verbose ) printf("      ECAL hit %4d :  E = %7.3f,  (x,y) = (%7.2f, %7.2f)   (dx,dy) = (%7.2f, %7.2f)   dr = %7.2f\n", hi, energy, hx, hy, dx, dy, dr ) ;

            } // inside radius?


         } // hi

         float ecal_ave_x = -999. ;
         float ecal_ave_y = -999. ;
         if ( ecal_nhit > 0 && ecal_sum_e > 0 ) {
            ecal_ave_x = ecal_sum_ex / ecal_sum_e ;
            ecal_ave_y = ecal_sum_ey / ecal_sum_e ;
         }
         float ecal_ave_dx = ecal_ave_x - tpex ;
         float ecal_ave_dy = ecal_ave_y - tpey ;
         float ecal_ave_dr = sqrt( ecal_ave_dx*ecal_ave_dx + ecal_ave_dy*ecal_ave_dy ) ;

         if ( verbose ) printf("    track %2d :  ECAL Esum = %7.3f, Nhit = %2d,  (avex,avey) = (%7.2f, %7.2f)  (dx,dy) = (%7.2f, %7.2f)   dr = %7.2f\n",
            ti, ecal_sum_e, ecal_nhit, ecal_ave_x, ecal_ave_y, ecal_ave_dx, ecal_ave_dy, ecal_ave_dr ) ;


         hist_map["ecal_sum_e"] -> Fill( ecal_sum_e ) ;
         hist_map["ecal_nhit"] -> Fill( ecal_nhit ) ;
         hist_map["ecal_dx"] -> Fill( ecal_ave_dx ) ;
         hist_map["ecal_dy"] -> Fill( ecal_ave_dy ) ;
         hist_map["ecal_dr"] -> Fill( ecal_ave_dr ) ;

         hist_map2d["ecal_nhit_vs_sum_e"] -> Fill( ecal_sum_e, ecal_nhit ) ;

         hist_map2d["ecal_nhit_vs_dx"] -> Fill( ecal_ave_dx, ecal_nhit ) ;




   //--- loop over HCAL hits

         int hcal_nhit(0) ;
         float hcal_sum_e(0.) ;
         float hcal_sum_ex(0.) ;
         float hcal_sum_ey(0.) ;

         for ( int hi=0; hi<fcs_rec_hcalE->size(); hi++ ) {

            float energy = fcs_rec_hcalE->at(hi) ;

            if ( energy < 0.01 ) continue ;

            float hx = fcs_rec_hcalX->at(hi) ;
            float hy = fcs_rec_hcalY->at(hi) ;

            float dx = hx - tphx ;
            float dy = hy - tphy ;
            float dr = sqrt( dx*dx + dy*dy ) ;

            if ( dr < radius ) {

               hcal_sum_e += energy ;
               hcal_sum_ex += hx * energy ;
               hcal_sum_ey += hy * energy ;

               hcal_nhit++ ;

               if ( verbose ) printf("      hcal hit %4d :  E = %7.3f,  (x,y) = (%7.2f, %7.2f)   (dx,dy) = (%7.2f, %7.2f)   dr = %7.2f\n", hi, energy, hx, hy, dx, dy, dr ) ;

            } // inside radius?


         } // hi

         float hcal_ave_x = -999. ;
         float hcal_ave_y = -999. ;
         if ( hcal_nhit > 0 && hcal_sum_e > 0 ) {
            hcal_ave_x = hcal_sum_ex / hcal_sum_e ;
            hcal_ave_y = hcal_sum_ey / hcal_sum_e ;
         }
         float hcal_ave_dx = hcal_ave_x - tphx ;
         float hcal_ave_dy = hcal_ave_y - tphy ;
         float hcal_ave_dr = sqrt( hcal_ave_dx*hcal_ave_dx + hcal_ave_dy*hcal_ave_dy ) ;

         if ( verbose ) printf("    track %2d :  hcal Esum = %7.3f, Nhit = %2d,  (avex,avey) = (%7.2f, %7.2f)  (dx,dy) = (%7.2f, %7.2f)   dr = %7.2f\n",
            ti, hcal_sum_e, hcal_nhit, hcal_ave_x, hcal_ave_y, hcal_ave_dx, hcal_ave_dy, hcal_ave_dr ) ;


         hist_map["hcal_sum_e"] -> Fill( hcal_sum_e ) ;
         hist_map["hcal_nhit"] -> Fill( hcal_nhit ) ;
         hist_map["hcal_dx"] -> Fill( hcal_ave_dx ) ;
         hist_map["hcal_dy"] -> Fill( hcal_ave_dy ) ;
         hist_map["hcal_dr"] -> Fill( hcal_ave_dr ) ;

         hist_map2d["hcal_nhit_vs_sum_e"] -> Fill( hcal_sum_e, hcal_nhit ) ;

         if ( tphx > 0 ) {
            hist_map["hcal_dx_xgt0"] -> Fill( hcal_ave_dx ) ;
         } else {
            hist_map["hcal_dx_xlt0"] -> Fill( hcal_ave_dx ) ;
         }


         if ( ecal_sum_e > 0.150 && ecal_sum_e < 0.350 && ecal_nhit > 0 && ecal_nhit < 3 ) {

            hist_map["hcal_ecalmip_sum_e"] -> Fill( hcal_sum_e ) ;
            hist_map["hcal_ecalmip_nhit"] -> Fill( hcal_nhit ) ;
            hist_map["hcal_ecalmip_dx"] -> Fill( hcal_ave_dx ) ;
            hist_map["hcal_ecalmip_dy"] -> Fill( hcal_ave_dy ) ;
            hist_map["hcal_ecalmip_dr"] -> Fill( hcal_ave_dr ) ;

            hist_map2d["hcal_ecalmip_nhit_vs_sum_e"] -> Fill( hcal_sum_e, hcal_nhit ) ;

         }


         if ( hcal_nhit > 0 ) {

            hist_map["ecal_hcalhits_sum_e"] -> Fill( ecal_sum_e ) ;
            hist_map["ecal_hcalhits_nhit"] -> Fill( ecal_nhit ) ;
            hist_map["ecal_hcalhits_dx"] -> Fill( ecal_ave_dx ) ;
            hist_map["ecal_hcalhits_dy"] -> Fill( ecal_ave_dy ) ;
            hist_map["ecal_hcalhits_dr"] -> Fill( ecal_ave_dr ) ;

            hist_map2d["ecal_hcalhits_nhit_vs_sum_e"] -> Fill( ecal_sum_e, ecal_nhit ) ;

         }



      } // ti

      if ( verbose ) printf(" End of mip_search1::doTrackLoop\n\n") ;

}

//---------------
void mip_search1::Loop()
{

   verbose = true ;
   verbose = false ;

   if (fChain == 0) return;

   radius = 20. ;

   gDirectory -> Delete( "h*" ) ;

   TH1F* h_ecal_sum_e = new TH1F( "h_ecal_sum_e", "h_ecal_sum_e", 100, 0., 2. ) ;
   TH1F* h_ecal_nhit  = new TH1F( "h_ecal_nhit", "h_ecal_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_ecal_dx = new TH1F( "h_ecal_dx", "h_ecal_dx", 60, -1*radius, radius ) ;
   TH1F* h_ecal_dy = new TH1F( "h_ecal_dy", "h_ecal_dy", 60, -1*radius, radius ) ;
   TH1F* h_ecal_dr = new TH1F( "h_ecal_dr", "h_ecal_dr", 60, 0., radius ) ;

   TH1F* h_ecal_hcalhits_sum_e = new TH1F( "h_ecal_hcalhits_sum_e", "h_ecal_hcalhits_sum_e", 100, 0., 2. ) ;
   TH1F* h_ecal_hcalhits_nhit  = new TH1F( "h_ecal_hcalhits_nhit", "h_ecal_hcalhits_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_ecal_hcalhits_dx = new TH1F( "h_ecal_hcalhits_dx", "h_ecal_hcalhits_dx", 60, -1*radius, radius ) ;
   TH1F* h_ecal_hcalhits_dy = new TH1F( "h_ecal_hcalhits_dy", "h_ecal_hcalhits_dy", 60, -1*radius, radius ) ;
   TH1F* h_ecal_hcalhits_dr = new TH1F( "h_ecal_hcalhits_dr", "h_ecal_hcalhits_dr", 60, 0., radius ) ;

   TH1F* h_hcal_sum_e = new TH1F( "h_hcal_sum_e", "h_hcal_sum_e", 100, 0., 5. ) ;
   TH1F* h_hcal_nhit  = new TH1F( "h_hcal_nhit", "h_hcal_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_hcal_dx = new TH1F( "h_hcal_dx", "h_hcal_dx", 60, -1*radius, radius ) ;
   TH1F* h_hcal_dy = new TH1F( "h_hcal_dy", "h_hcal_dy", 60, -1*radius, radius ) ;
   TH1F* h_hcal_dr = new TH1F( "h_hcal_dr", "h_hcal_dr", 60, 0., radius ) ;
   TH1F* h_hcal_dx_xgt0 = new TH1F( "h_hcal_dx_xgt0", "h_hcal_dx_xgt0", 60, -1*radius, radius ) ;
   TH1F* h_hcal_dx_xlt0 = new TH1F( "h_hcal_dx_xlt0", "h_hcal_dx_xlt0", 60, -1*radius, radius ) ;


   TH1F* h_hcal_ecalmip_sum_e = new TH1F( "h_hcal_ecalmip_sum_e", "h_hcal_ecalmip_sum_e", 100, 0., 5. ) ;
   TH1F* h_hcal_ecalmip_nhit  = new TH1F( "h_hcal_ecalmip_nhit", "h_hcal_ecalmip_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_hcal_ecalmip_dx = new TH1F( "h_hcal_ecalmip_dx", "h_hcal_ecalmip_dx", 60, -1*radius, radius ) ;
   TH1F* h_hcal_ecalmip_dy = new TH1F( "h_hcal_ecalmip_dy", "h_hcal_ecalmip_dy", 60, -1*radius, radius ) ;
   TH1F* h_hcal_ecalmip_dr = new TH1F( "h_hcal_ecalmip_dr", "h_hcal_ecalmip_dr", 60, 0., radius ) ;

   TH2F* h_ecal_nhit_vs_sum_e = new TH2F( "h_ecal_nhit_vs_sum_e", "h_ecal_nhit_vs_sum_e", 100, 0., 2., 21, -0.5, 20.5 ) ;
   TH2F* h_hcal_nhit_vs_sum_e = new TH2F( "h_hcal_nhit_vs_sum_e", "h_hcal_nhit_vs_sum_e", 100, 0., 5., 21, -0.5, 20.5 ) ;
   TH2F* h_hcal_ecalmip_nhit_vs_sum_e = new TH2F( "h_hcal_ecalmip_nhit_vs_sum_e", "h_hcal_ecalmip_nhit_vs_sum_e", 100, 0., 5., 21, -0.5, 20.5 ) ;
   TH2F* h_ecal_hcalhits_nhit_vs_sum_e = new TH2F( "h_ecal_hcalhits_nhit_vs_sum_e", "h_ecal_hcalhits_nhit_vs_sum_e", 100, 0., 2., 21, -0.5, 20.5 ) ;

   TH2F* h_ecal_nhit_vs_dx = new TH2F( "h_ecal_nhit_vs_dx", "h_ecal_nhit_vs_dx", 60, -1*radius, radius, 21, -0.5, 20.5 ) ;


   TH1F* h_mix_ecal_sum_e = new TH1F( "h_mix_ecal_sum_e", "h_mix_ecal_sum_e", 100, 0., 2. ) ;
   TH1F* h_mix_ecal_nhit  = new TH1F( "h_mix_ecal_nhit", "h_mix_ecal_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_mix_ecal_dx = new TH1F( "h_mix_ecal_dx", "h_mix_ecal_dx", 60, -1*radius, radius ) ;
   TH1F* h_mix_ecal_dy = new TH1F( "h_mix_ecal_dy", "h_mix_ecal_dy", 60, -1*radius, radius ) ;
   TH1F* h_mix_ecal_dr = new TH1F( "h_mix_ecal_dr", "h_mix_ecal_dr", 60, 0., radius ) ;

   TH1F* h_mix_ecal_hcalhits_sum_e = new TH1F( "h_mix_ecal_hcalhits_sum_e", "h_mix_ecal_hcalhits_sum_e", 100, 0., 2. ) ;
   TH1F* h_mix_ecal_hcalhits_nhit  = new TH1F( "h_mix_ecal_hcalhits_nhit", "h_mix_ecal_hcalhits_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_mix_ecal_hcalhits_dx = new TH1F( "h_mix_ecal_hcalhits_dx", "h_mix_ecal_hcalhits_dx", 60, -1*radius, radius ) ;
   TH1F* h_mix_ecal_hcalhits_dy = new TH1F( "h_mix_ecal_hcalhits_dy", "h_mix_ecal_hcalhits_dy", 60, -1*radius, radius ) ;
   TH1F* h_mix_ecal_hcalhits_dr = new TH1F( "h_mix_ecal_hcalhits_dr", "h_mix_ecal_hcalhits_dr", 60, 0., radius ) ;

   TH1F* h_mix_hcal_sum_e = new TH1F( "h_mix_hcal_sum_e", "h_mix_hcal_sum_e", 100, 0., 5. ) ;
   TH1F* h_mix_hcal_nhit  = new TH1F( "h_mix_hcal_nhit", "h_mix_hcal_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_mix_hcal_dx = new TH1F( "h_mix_hcal_dx", "h_mix_hcal_dx", 60, -1*radius, radius ) ;
   TH1F* h_mix_hcal_dy = new TH1F( "h_mix_hcal_dy", "h_mix_hcal_dy", 60, -1*radius, radius ) ;
   TH1F* h_mix_hcal_dr = new TH1F( "h_mix_hcal_dr", "h_mix_hcal_dr", 60, 0., radius ) ;
   TH1F* h_mix_hcal_dx_xgt0 = new TH1F( "h_mix_hcal_dx_xgt0", "h_mix_hcal_dx_xgt0", 60, -1*radius, radius ) ;
   TH1F* h_mix_hcal_dx_xlt0 = new TH1F( "h_mix_hcal_dx_xlt0", "h_mix_hcal_dx_xlt0", 60, -1*radius, radius ) ;


   TH1F* h_mix_hcal_ecalmip_sum_e = new TH1F( "h_mix_hcal_ecalmip_sum_e", "h_mix_hcal_ecalmip_sum_e", 100, 0., 5. ) ;
   TH1F* h_mix_hcal_ecalmip_nhit  = new TH1F( "h_mix_hcal_ecalmip_nhit", "h_mix_hcal_ecalmip_nhit", 21, -0.5, 20.5 ) ;
   TH1F* h_mix_hcal_ecalmip_dx = new TH1F( "h_mix_hcal_ecalmip_dx", "h_mix_hcal_ecalmip_dx", 60, -1*radius, radius ) ;
   TH1F* h_mix_hcal_ecalmip_dy = new TH1F( "h_mix_hcal_ecalmip_dy", "h_mix_hcal_ecalmip_dy", 60, -1*radius, radius ) ;
   TH1F* h_mix_hcal_ecalmip_dr = new TH1F( "h_mix_hcal_ecalmip_dr", "h_mix_hcal_ecalmip_dr", 60, 0., radius ) ;

   TH2F* h_mix_ecal_nhit_vs_sum_e = new TH2F( "h_mix_ecal_nhit_vs_sum_e", "h_mix_ecal_nhit_vs_sum_e", 100, 0., 2., 21, -0.5, 20.5 ) ;
   TH2F* h_mix_hcal_nhit_vs_sum_e = new TH2F( "h_mix_hcal_nhit_vs_sum_e", "h_mix_hcal_nhit_vs_sum_e", 100, 0., 5., 21, -0.5, 20.5 ) ;
   TH2F* h_mix_hcal_ecalmip_nhit_vs_sum_e = new TH2F( "h_mix_hcal_ecalmip_nhit_vs_sum_e", "h_mix_hcal_ecalmip_nhit_vs_sum_e", 100, 0., 5., 21, -0.5, 20.5 ) ;
   TH2F* h_mix_ecal_hcalhits_nhit_vs_sum_e = new TH2F( "h_mix_ecal_hcalhits_nhit_vs_sum_e", "h_mix_ecal_hcalhits_nhit_vs_sum_e", 100, 0., 2., 21, -0.5, 20.5 ) ;

   TH2F* h_mix_ecal_nhit_vs_dx = new TH2F( "h_mix_ecal_nhit_vs_dx", "h_mix_ecal_nhit_vs_dx", 60, -1*radius, radius, 21, -0.5, 20.5 ) ;


   Long64_t nentries = fChain->GetEntries();
   printf("\n\n Chain has %llu entries\n\n", nentries ) ;



   std::map<const char*,TH1F*> hist_map_same ;
   std::map<const char*,TH1F*> hist_map_mix ;

   hist_map_same["ecal_sum_e"] = h_ecal_sum_e ;
   hist_map_same["ecal_nhit"]  = h_ecal_nhit ;
   hist_map_same["ecal_dx"]    = h_ecal_dx ;
   hist_map_same["ecal_dy"]    = h_ecal_dy ;
   hist_map_same["ecal_dr"]    = h_ecal_dr ;

   hist_map_mix["ecal_sum_e"] = h_mix_ecal_sum_e ;
   hist_map_mix["ecal_nhit"]  = h_mix_ecal_nhit ;
   hist_map_mix["ecal_dx"]    = h_mix_ecal_dx ;
   hist_map_mix["ecal_dy"]    = h_mix_ecal_dy ;
   hist_map_mix["ecal_dr"]    = h_mix_ecal_dr ;


   hist_map_same["ecal_hcalhits_sum_e"] = h_ecal_hcalhits_sum_e ;
   hist_map_same["ecal_hcalhits_nhit"]  = h_ecal_hcalhits_nhit ;
   hist_map_same["ecal_hcalhits_dx"]    = h_ecal_hcalhits_dx ;
   hist_map_same["ecal_hcalhits_dy"]    = h_ecal_hcalhits_dy ;
   hist_map_same["ecal_hcalhits_dr"]    = h_ecal_hcalhits_dr ;

   hist_map_mix["ecal_hcalhits_sum_e"] = h_mix_ecal_hcalhits_sum_e ;
   hist_map_mix["ecal_hcalhits_nhit"]  = h_mix_ecal_hcalhits_nhit ;
   hist_map_mix["ecal_hcalhits_dx"]    = h_mix_ecal_hcalhits_dx ;
   hist_map_mix["ecal_hcalhits_dy"]    = h_mix_ecal_hcalhits_dy ;
   hist_map_mix["ecal_hcalhits_dr"]    = h_mix_ecal_hcalhits_dr ;



   hist_map_same["hcal_sum_e"] = h_hcal_sum_e ;
   hist_map_same["hcal_nhit"]  = h_hcal_nhit ;
   hist_map_same["hcal_dx"]    = h_hcal_dx ;
   hist_map_same["hcal_dy"]    = h_hcal_dy ;
   hist_map_same["hcal_dr"]    = h_hcal_dr ;
   hist_map_same["hcal_dx_xgt0"]    = h_hcal_dx_xgt0 ;
   hist_map_same["hcal_dx_xlt0"]    = h_hcal_dx_xlt0 ;

   hist_map_same["hcal_ecalmip_sum_e"] = h_hcal_ecalmip_sum_e ;
   hist_map_same["hcal_ecalmip_nhit"]  = h_hcal_ecalmip_nhit ;
   hist_map_same["hcal_ecalmip_dx"]    = h_hcal_ecalmip_dx ;
   hist_map_same["hcal_ecalmip_dy"]    = h_hcal_ecalmip_dy ;
   hist_map_same["hcal_ecalmip_dr"]    = h_hcal_ecalmip_dr ;

   hist_map_mix["hcal_sum_e"] = h_mix_hcal_sum_e ;
   hist_map_mix["hcal_nhit"]  = h_mix_hcal_nhit ;
   hist_map_mix["hcal_dx"]    = h_mix_hcal_dx ;
   hist_map_mix["hcal_dy"]    = h_mix_hcal_dy ;
   hist_map_mix["hcal_dr"]    = h_mix_hcal_dr ;
   hist_map_mix["hcal_dx_xgt0"]    = h_mix_hcal_dx_xgt0 ;
   hist_map_mix["hcal_dx_xlt0"]    = h_mix_hcal_dx_xlt0 ;

   hist_map_mix["hcal_ecalmip_sum_e"] = h_mix_hcal_ecalmip_sum_e ;
   hist_map_mix["hcal_ecalmip_nhit"]  = h_mix_hcal_ecalmip_nhit ;
   hist_map_mix["hcal_ecalmip_dx"]    = h_mix_hcal_ecalmip_dx ;
   hist_map_mix["hcal_ecalmip_dy"]    = h_mix_hcal_ecalmip_dy ;
   hist_map_mix["hcal_ecalmip_dr"]    = h_mix_hcal_ecalmip_dr ;



   std::map<const char*,TH2F*> hist_map2d_same ;
   std::map<const char*,TH2F*> hist_map2d_mix ;

   hist_map2d_same["ecal_nhit_vs_sum_e"] = h_ecal_nhit_vs_sum_e ;
   hist_map2d_same["hcal_nhit_vs_sum_e"] = h_hcal_nhit_vs_sum_e ;
   hist_map2d_same["hcal_ecalmip_nhit_vs_sum_e"] = h_hcal_ecalmip_nhit_vs_sum_e ;
   hist_map2d_same["ecal_hcalhits_nhit_vs_sum_e"] = h_ecal_hcalhits_nhit_vs_sum_e ;
   hist_map2d_same["ecal_nhit_vs_dx"] = h_ecal_nhit_vs_dx ;

   hist_map2d_mix["ecal_nhit_vs_sum_e"] = h_mix_ecal_nhit_vs_sum_e ;
   hist_map2d_mix["hcal_nhit_vs_sum_e"] = h_mix_hcal_nhit_vs_sum_e ;
   hist_map2d_mix["hcal_ecalmip_nhit_vs_sum_e"] = h_mix_hcal_ecalmip_nhit_vs_sum_e ;
   hist_map2d_mix["ecal_hcalhits_nhit_vs_sum_e"] = h_mix_ecal_hcalhits_nhit_vs_sum_e ;
   hist_map2d_mix["ecal_nhit_vs_dx"] = h_mix_ecal_nhit_vs_dx ;




   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if ( jentry % 1000 == 0 ) printf("  %9llu / %9llu\n", jentry, nentries ) ;

      if ( verbose ) printf("\n\n ========= Event %llu\n", jentry ) ;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      doTrackLoop( hist_map_same, hist_map2d_same,  false ) ;
      doTrackLoop( hist_map_mix, hist_map2d_mix, true ) ;

      if ( verbose ) {
            char answ = getchar() ;
            if ( answ == 'q' ) return ;
      }

      previous_event_rcProjEcalx->clear() ;
      previous_event_rcProjEcaly->clear() ;
      previous_event_rcProjHcalx->clear() ;
      previous_event_rcProjHcaly->clear() ;
      previous_event_rcNumFST->clear() ;
      for ( int ti=0; ti<rcN; ti++ ) {
         previous_event_rcProjEcalx -> push_back( rcProjEcalx->at(ti) ) ;
         previous_event_rcProjEcaly -> push_back( rcProjEcaly->at(ti) ) ;
         previous_event_rcProjHcalx -> push_back( rcProjHcalx->at(ti) ) ;
         previous_event_rcProjHcaly -> push_back( rcProjHcaly->at(ti) ) ;
         previous_event_rcNumFST -> push_back( rcNumFST->at(ti) ) ;
      } // ti



   } // jentry


   gDirectory -> ls() ;

   saveHist( "mip-search-hists.root", "h*" ) ;

} // Loop







