
#include "histio.c"
#include "utils.c"

void draw_mip_plots1( const char* infile = "mip-search-hists-all.root" ) {

   gDirectory -> Delete("h*" ) ;

   loadHist( infile ) ;

   gDirectory -> ls() ;

   gStyle -> SetOptStat(0) ;

   TCanvas* can1 = get_canvas( "can1", "can1", 50, 50, 1100, 1100 ) ;

   can1 -> cd() ;
   can1 -> Clear() ;

   can1 -> Divide( 2, 2 ) ;

   TH1F* h_ecal_sum_e = get_hist( "h_ecal_sum_e" ) ;
   TH1F* h_mix_ecal_sum_e = get_hist( "h_mix_ecal_sum_e" ) ;

   h_ecal_sum_e -> SetFillColor( 30 ) ;
   h_mix_ecal_sum_e -> SetFillColor( 32 ) ;

   h_ecal_sum_e -> SetMaximum( 1.2*(h_ecal_sum_e->GetBinContent(12) ) ) ;
   h_mix_ecal_sum_e -> SetMaximum( 1.2*(h_mix_ecal_sum_e->GetBinContent(3) ) ) ;

   h_ecal_sum_e -> SetXTitle( "ECAL Energy sum (GeV)" ) ;
   h_mix_ecal_sum_e -> SetXTitle( "ECAL Energy sum (GeV)" ) ;

   h_ecal_sum_e -> SetTitle( "Tracking and FCS from same event") ;
   h_mix_ecal_sum_e -> SetTitle( "Uncorrelated tracking and FCS") ;




   TH1F* h_ecal_nhit = get_hist( "h_ecal_nhit" ) ;
   TH1F* h_mix_ecal_nhit = get_hist( "h_mix_ecal_nhit" ) ;

   h_ecal_nhit -> SetFillColor( 30 ) ;
   h_mix_ecal_nhit -> SetFillColor( 32 ) ;

   //////////h_ecal_nhit -> SetMaximum( 1.2*(h_ecal_nhit->GetBinContent(12) ) ) ;
   //////////h_mix_ecal_nhit -> SetMaximum( 1.2*(h_mix_ecal_nhit->GetBinContent(3) ) ) ;

   h_ecal_nhit -> SetXTitle( "ECAL Nhits (r<20 cm)" ) ;
   h_mix_ecal_nhit -> SetXTitle( "ECAL Nhits (r<20 cm)" ) ;

   h_ecal_nhit -> SetTitle( "Tracking and FCS from same event") ;
   h_mix_ecal_nhit -> SetTitle( "Uncorrelated tracking and FCS") ;



   can1 -> cd(1) ;
   h_ecal_sum_e -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_sum_e -> Draw("same axig") ;

   can1 -> cd(2) ;
   h_mix_ecal_sum_e -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_sum_e -> Draw("same axig") ;

   can1 -> cd(3) ;
   h_ecal_nhit -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_nhit -> Draw("same axig") ;

   can1 -> cd(4) ;
   h_mix_ecal_nhit -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_nhit -> Draw("same axig") ;



}





