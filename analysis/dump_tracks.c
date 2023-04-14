#define dump_tracks_cxx
#include "dump_tracks.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void dump_tracks::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      printf("\n\n ====== Event %5llu\n\n", jentry ) ;

      printf("  Reconstructed tracks:  %3d\n", rcN ) ;

      for ( int ti=0; ti<rcN; ti++ ) {
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
         printf("   MC : Eta = %7.3f, Phi = %7.3f, Pt = %7.3f,  q = %2d,  Id = %3d, PID = %8d, dr = %7.3f\n",
                mcEta->at(mcind), mcPhi->at(mcind), mcPt->at(mcind), mcCharge->at(mcind), mcind, mcPid->at(mcind), mindr ) ;
         int tpi = rcTrackIdT->at(ti) ;
         for ( int tpi=0; tpi<tprojN; tpi++ ) {
            if ( tprojIdT->at(tpi) != rcTrackIdT->at(ti) ) continue ;
            if ( tprojIdD -> at(tpi) != 4 ) continue ;
            float pt = sqrt( tprojPx->at(tpi) * tprojPx->at(tpi) + tprojPy->at(tpi) * tprojPy->at(tpi) ) ;
            float theta = atan2( pt, tprojPz->at(tpi) ) ;
            float phi = atan2( tprojPy->at(tpi), tprojPx->at(tpi) ) ;
            float eta = -1 * log( tan( theta/2. ) ) ;
      ///   printf(" tproj:  %3d  idt = %5d,  x,y,z = (%7.2f, %7.2f, %7.2f)  Eta = %7.3f, phi = %7.3f, pt = %7.2f\n",
      ///      tpi, tprojIdT->at(tpi), tprojX->at(tpi), tprojY->at(tpi), tprojZ->at(tpi), eta, phi, pt ) ;
            printf(" tproj: Eta = %7.3f, Phi = %7.3f, Pt = %7.3f  idt = %5d,  x,y,z = (%7.2f, %7.2f, %7.2f)\n",
               eta, phi, pt, tprojIdT->at(tpi), tprojX->at(tpi), tprojY->at(tpi), tprojZ->at(tpi) ) ;
         } // tpi
         printf("\n") ;
      } // ti

///   for ( int tpi=0; tpi<tprojN; tpi++ ) {
///      if ( tprojIdD -> at(tpi) != 4 ) continue ;
///      float pt = sqrt( tprojPx->at(tpi) * tprojPx->at(tpi) + tprojPy->at(tpi) * tprojPy->at(tpi) ) ;
///      float theta = atan2( pt, tprojPz->at(tpi) ) ;
///      float phi = atan2( tprojPy->at(tpi), tprojPx->at(tpi) ) ;
///      float eta = -1 * log( tan( theta/2. ) ) ;
///      printf(" tproj:  %3d  idt = %5d,  x,y,z = (%7.2f, %7.2f, %7.2f)  Eta = %7.3f, phi = %7.3f, pt = %7.2f\n",
///         tpi, tprojIdT->at(tpi), tprojX->at(tpi), tprojY->at(tpi), tprojZ->at(tpi), eta, phi, pt ) ;
///   } // tpi



   } // jentry
}
