
#include "histio.c"
#include "utils.c"

void draw_bur_plots1( const char* infile = "mip-search-hists-all.root" ) {

   gSystem -> Exec( "mkdir -p plots" ) ;

   gDirectory -> Delete("h*" ) ;
   gStyle -> SetPadLeftMargin(0.15) ;
   gStyle -> SetPadBottomMargin(0.15) ;

   loadHist( infile ) ;

   gDirectory -> ls() ;

   gStyle -> SetOptStat(0) ;

   TCanvas* can1 = get_canvas( "can1", "can1", 50, 50, 800, 800 ) ;

   can1 -> cd() ;
   can1 -> Clear() ;


   TH1F* h_ecal_sum_e = get_hist( "h_ecal_sum_e" ) ;
   TH1F* h_mix_ecal_sum_e = get_hist( "h_mix_ecal_sum_e" ) ;

   TH2F* h_ecal_nhit_vs_sum_e = get_hist2d( "h_ecal_nhit_vs_sum_e" ) ;
   TH2F* h_mix_ecal_nhit_vs_sum_e = get_hist2d( "h_mix_ecal_nhit_vs_sum_e" ) ;

   TH1D* h_ecal_sum_e_nhit12 = h_ecal_nhit_vs_sum_e -> ProjectionX( "h_ecal_sum_e_nhit12", 2, 3 ) ;
   TH1D* h_mix_ecal_sum_e_nhit12 = h_mix_ecal_nhit_vs_sum_e -> ProjectionX( "h_mix_ecal_sum_e_nhit12", 2, 3 ) ;

   h_ecal_sum_e_nhit12 -> SetFillColor(30) ;
   h_mix_ecal_sum_e_nhit12 -> SetFillColor(32) ;

   h_ecal_sum_e -> SetFillColor( 30 ) ;
   h_mix_ecal_sum_e -> SetFillColor( 32 ) ;

   h_ecal_sum_e -> SetMaximum( 1.2*(h_ecal_sum_e->GetBinContent(12) ) ) ;
   h_mix_ecal_sum_e -> SetMaximum( 1.2*(h_mix_ecal_sum_e->GetBinContent(3) ) ) ;

   h_ecal_sum_e -> SetXTitle( "ECAL Energy sum (GeV)" ) ;
   h_mix_ecal_sum_e -> SetXTitle( "ECAL Energy sum (GeV)" ) ;

   h_ecal_sum_e -> SetYTitle( "Tracks / 20 MeV" ) ;
   h_mix_ecal_sum_e -> SetYTitle( "Tracks / 20 MeV" ) ;

   h_ecal_sum_e -> SetTitle( "Tracking and FCS from same event") ;
   h_mix_ecal_sum_e -> SetTitle( "Uncorrelated tracking and FCS") ;


   h_ecal_sum_e_nhit12 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;
   h_mix_ecal_sum_e_nhit12 -> SetXTitle( "ECAL Energy sum (GeV)" ) ;

   h_ecal_sum_e_nhit12 -> SetYTitle( "Tracks / 20 MeV" ) ;
   h_mix_ecal_sum_e_nhit12 -> SetYTitle( "Tracks / 20 MeV" ) ;

   h_ecal_sum_e_nhit12 -> SetTitle( "Tracking and FCS from same event, N ECAL = 1 or 2") ;
   h_mix_ecal_sum_e_nhit12 -> SetTitle( "Uncorrelated tracking and FCS, N ECAL = 1 or 2") ;


   char a ;




   //----------

   h_ecal_sum_e -> SetTitleOffset( 1.2, "x" ) ;
   h_ecal_sum_e -> SetTitleOffset( 1.8, "y" ) ;

   h_ecal_sum_e -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_sum_e -> Draw("same axig") ;

   can1 -> SaveAs("plots/ecal-energy.pdf" ) ;
   can1 -> SaveAs("plots/ecal-energy.png" ) ;

   can1 -> Update() ;
   can1 -> Draw() ;
   gSystem -> ProcessEvents() ;

   a = getchar() ;

   //----------

   h_mix_ecal_sum_e -> SetTitleOffset( 1.2, "x" ) ;
   h_mix_ecal_sum_e -> SetTitleOffset( 1.8, "y" ) ;

   h_mix_ecal_sum_e -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_sum_e -> Draw("same axig") ;

   can1 -> SaveAs("plots/ecal-energy-uncorrelated-trk-fcs.pdf" ) ;
   can1 -> SaveAs("plots/ecal-energy-uncorrelated-trk-fcs.png" ) ;

   can1 -> Update() ;
   can1 -> Draw() ;
   gSystem -> ProcessEvents() ;

   a = getchar() ;



   //----------

   h_ecal_sum_e_nhit12 -> SetTitleOffset( 1.2, "x" ) ;
   h_ecal_sum_e_nhit12 -> SetTitleOffset( 1.8, "y" ) ;
   h_ecal_sum_e_nhit12 -> GetXaxis() -> SetRange(1,50) ;

   h_ecal_sum_e_nhit12 -> Draw() ;
   gPad->SetGridx(1) ;
   h_ecal_sum_e_nhit12 -> Draw("same axig") ;

   can1 -> SaveAs("plots/ecal-energy-nemchit12.pdf" ) ;
   can1 -> SaveAs("plots/ecal-energy-nemchit12.png" ) ;

   can1 -> Update() ;
   can1 -> Draw() ;
   gSystem -> ProcessEvents() ;

   a = getchar() ;

   //----------

   h_mix_ecal_sum_e_nhit12 -> SetTitleOffset( 1.2, "x" ) ;
   h_mix_ecal_sum_e_nhit12 -> SetTitleOffset( 1.8, "y" ) ;
   h_mix_ecal_sum_e_nhit12 -> GetXaxis() -> SetRange(1,50) ;

   h_mix_ecal_sum_e_nhit12 -> Draw() ;
   gPad->SetGridx(1) ;
   h_mix_ecal_sum_e_nhit12 -> Draw("same axig") ;

   can1 -> SaveAs("plots/ecal-energy-nemchit12-uncorrelated-trk-fcs.pdf" ) ;
   can1 -> SaveAs("plots/ecal-energy-nemchit12-uncorrelated-trk-fcs.png" ) ;

   can1 -> Update() ;
   can1 -> Draw() ;
   gSystem -> ProcessEvents() ;

   a = getchar() ;






}





