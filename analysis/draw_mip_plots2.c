
#include "histio.c"
#include "utils.c"

void draw_mip_plots2( const char* infile = "mip-search-hists-all.root" ) {

   gDirectory -> Delete("h*" ) ;

   loadHist( infile ) ;

   gDirectory -> ls() ;

   gStyle -> SetOptStat(0) ;

   TCanvas* can1 = get_canvas( "can1", "can1", 50, 50, 800, 1100 ) ;

   can1 -> cd() ;
   can1 -> Clear() ;

   can1 -> Divide( 2, 4 ) ;

   TH2F* h_ecal_nhit_vs_sum_e = get_hist2d( "h_ecal_nhit_vs_sum_e" ) ;
   TH2F* h_mix_ecal_nhit_vs_sum_e = get_hist2d( "h_mix_ecal_nhit_vs_sum_e" ) ;

   TH1D* h_ecal_sum_e_nhit1 = h_ecal_nhit_vs_sum_e -> ProjectionX( "h_ecal_sum_e_nhit1", 2, 2 ) ;
   TH1D* h_mix_ecal_sum_e_nhit1 = h_mix_ecal_nhit_vs_sum_e -> ProjectionX( "h_mix_ecal_sum_e_nhit1", 2, 2 ) ;

   TH1D* h_ecal_sum_e_nhit2 = h_ecal_nhit_vs_sum_e -> ProjectionX( "h_ecal_sum_e_nhit2", 3, 3 ) ;
   TH1D* h_mix_ecal_sum_e_nhit2 = h_mix_ecal_nhit_vs_sum_e -> ProjectionX( "h_mix_ecal_sum_e_nhit2", 3, 3 ) ;

   TH1D* h_ecal_sum_e_nhit3 = h_ecal_nhit_vs_sum_e -> ProjectionX( "h_ecal_sum_e_nhit3", 4, 4 ) ;
   TH1D* h_mix_ecal_sum_e_nhit3 = h_mix_ecal_nhit_vs_sum_e -> ProjectionX( "h_mix_ecal_sum_e_nhit3", 4, 4 ) ;

   TH1D* h_ecal_sum_e_nhit4 = h_ecal_nhit_vs_sum_e -> ProjectionX( "h_ecal_sum_e_nhit4", 5, 5 ) ;
   TH1D* h_mix_ecal_sum_e_nhit4 = h_mix_ecal_nhit_vs_sum_e -> ProjectionX( "h_mix_ecal_sum_e_nhit4", 5, 5 ) ;



   h_ecal_sum_e_nhit1 -> SetFillColor( 30 ) ;
   h_mix_ecal_sum_e_nhit1 -> SetFillColor( 32 ) ;

   //h_ecal_sum_e_nhit1 -> SetMaximum( 1.2*(h_ecal_sum_e_nhit1->GetBinContent(3) ) ) ;
   //h_mix_ecal_sum_e_nhit1 -> SetMaximum( 1.2*(h_mix_ecal_sum_e_nhit1->GetBinContent(3) ) ) ;

   h_ecal_sum_e_nhit1 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;
   h_mix_ecal_sum_e_nhit1 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;

   h_ecal_sum_e_nhit1 -> SetTitle( "Tracking and FCS from same event, 1 hit") ;
   h_mix_ecal_sum_e_nhit1 -> SetTitle( "Uncorrelated tracking and FCS, 1 hit") ;



   h_ecal_sum_e_nhit2 -> SetFillColor( 30 ) ;
   h_mix_ecal_sum_e_nhit2 -> SetFillColor( 32 ) ;

   //h_ecal_sum_e_nhit2 -> SetMaximum( 1.2*(h_ecal_sum_e_nhit2->GetBinContent(3) ) ) ;
   //h_mix_ecal_sum_e_nhit2 -> SetMaximum( 1.2*(h_mix_ecal_sum_e_nhit2->GetBinContent(3) ) ) ;

   h_ecal_sum_e_nhit2 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;
   h_mix_ecal_sum_e_nhit2 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;

   h_ecal_sum_e_nhit2 -> SetTitle( "Tracking and FCS from same event, 2 hits") ;
   h_mix_ecal_sum_e_nhit2 -> SetTitle( "Uncorrelated tracking and FCS, 2 hits") ;



   h_ecal_sum_e_nhit3 -> SetFillColor( 30 ) ;
   h_mix_ecal_sum_e_nhit3 -> SetFillColor( 32 ) ;

   //h_ecal_sum_e_nhit3 -> SetMaximum( 1.2*(h_ecal_sum_e_nhit3->GetBinContent(3) ) ) ;
   //h_mix_ecal_sum_e_nhit3 -> SetMaximum( 1.2*(h_mix_ecal_sum_e_nhit3->GetBinContent(3) ) ) ;

   h_ecal_sum_e_nhit3 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;
   h_mix_ecal_sum_e_nhit3 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;

   h_ecal_sum_e_nhit3 -> SetTitle( "Tracking and FCS from same event, 3 hits") ;
   h_mix_ecal_sum_e_nhit3 -> SetTitle( "Uncorrelated tracking and FCS, 3 hits") ;



   h_ecal_sum_e_nhit4 -> SetFillColor( 30 ) ;
   h_mix_ecal_sum_e_nhit4 -> SetFillColor( 32 ) ;

   //h_ecal_sum_e_nhit4 -> SetMaximum( 1.2*(h_ecal_sum_e_nhit4->GetBinContent(3) ) ) ;
   //h_mix_ecal_sum_e_nhit4 -> SetMaximum( 1.2*(h_mix_ecal_sum_e_nhit4->GetBinContent(3) ) ) ;

   h_ecal_sum_e_nhit4 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;
   h_mix_ecal_sum_e_nhit4 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;

   h_ecal_sum_e_nhit4 -> SetTitle( "Tracking and FCS from same event, 4 hits") ;
   h_mix_ecal_sum_e_nhit4 -> SetTitle( "Uncorrelated tracking and FCS, 4 hits") ;







   can1 -> cd(1) ;
   h_ecal_sum_e_nhit1 -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_sum_e_nhit1 -> Draw("same axig") ;

   can1 -> cd(2) ;
   h_mix_ecal_sum_e_nhit1 -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_sum_e_nhit1 -> Draw("same axig") ;


   can1 -> cd(3) ;
   h_ecal_sum_e_nhit2 -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_sum_e_nhit2 -> Draw("same axig") ;

   can1 -> cd(4) ;
   h_mix_ecal_sum_e_nhit2 -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_sum_e_nhit2 -> Draw("same axig") ;


   can1 -> cd(5) ;
   h_ecal_sum_e_nhit3 -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_sum_e_nhit3 -> Draw("same axig") ;

   can1 -> cd(6) ;
   h_mix_ecal_sum_e_nhit3 -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_sum_e_nhit3 -> Draw("same axig") ;


   can1 -> cd(7) ;
   h_ecal_sum_e_nhit4 -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_sum_e_nhit4 -> Draw("same axig") ;

   can1 -> cd(8) ;
   h_mix_ecal_sum_e_nhit4 -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_sum_e_nhit4 -> Draw("same axig") ;




}





