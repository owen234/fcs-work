
#include "histio.c"
#include "utils.c"

void draw_mip_plots3( const char* infile = "mip-search-hists-all.root" ) {

   gDirectory -> Delete("h*" ) ;

   loadHist( infile ) ;

   gDirectory -> ls() ;

   gStyle -> SetOptStat(0) ;

   TCanvas* can1 = get_canvas( "can1", "can1", 50, 50, 1100, 1100 ) ;

   can1 -> cd() ;
   can1 -> Clear() ;

   can1 -> Divide( 2, 2 ) ;

   TH1F* h_ecal_dx = get_hist( "h_ecal_dx" ) ;
   TH1F* h_mix_ecal_dx = get_hist( "h_mix_ecal_dx" ) ;

   h_ecal_dx -> SetFillColor( 30 ) ;
   h_mix_ecal_dx -> SetFillColor( 32 ) ;

   //h_ecal_dx -> SetMaximum( 1.2*(h_ecal_dx->GetBinContent(12) ) ) ;
   //h_mix_ecal_dx -> SetMaximum( 1.2*(h_mix_ecal_dx->GetBinContent(3) ) ) ;

   h_ecal_dx -> SetXTitle( "ECAL dx (calo-track)" ) ;
   h_mix_ecal_dx -> SetXTitle( "ECAL dx (calo-track)" ) ;

   h_ecal_dx -> SetTitle( "Tracking and FCS from same event") ;
   h_mix_ecal_dx -> SetTitle( "Uncorrelated tracking and FCS") ;




   TH1F* h_ecal_dy = get_hist( "h_ecal_dy" ) ;
   TH1F* h_mix_ecal_dy = get_hist( "h_mix_ecal_dy" ) ;

   h_ecal_dy -> SetFillColor( 30 ) ;
   h_mix_ecal_dy -> SetFillColor( 32 ) ;

   //h_ecal_dy -> SetMaximum( 1.2*(h_ecal_dy->GetBinContent(12) ) ) ;
   //h_mix_ecal_dy -> SetMaximum( 1.2*(h_mix_ecal_dy->GetBinContent(3) ) ) ;

   h_ecal_dy -> SetXTitle( "ECAL dy (calo-track)" ) ;
   h_mix_ecal_dy -> SetXTitle( "ECAL dy (calo-track)" ) ;

   h_ecal_dy -> SetTitle( "Tracking and FCS from same event") ;
   h_mix_ecal_dy -> SetTitle( "Uncorrelated tracking and FCS") ;



   can1 -> cd(1) ;
   h_ecal_dx -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_dx -> Draw("same axig") ;

   can1 -> cd(2) ;
   h_mix_ecal_dx -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_dx -> Draw("same axig") ;

   can1 -> cd(3) ;
   h_ecal_dy -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_dy -> Draw("same axig") ;

   can1 -> cd(4) ;
   h_mix_ecal_dy -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_dy -> Draw("same axig") ;



}





