// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim_pythia6.C

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class St_db_Maker;
St_db_Maker *db_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *primaryMaker = 0;

class StarPythia6;
StarPythia6 *pythia6 = 0;

class StarFilterMaker;
class FcsDYFilter;
class FcsJetFilter;
class FcsMuonFilter;
StarFilterMaker *dyfilter =0;
StarFilterMaker *dybgfilter=0;
StarFilterMaker *jetfilter=0;
StarFilterMaker *jpsifilter=0;
StarFilterMaker *muonfilter=0;

TString LHAPDF_DATA_PATH="/star/u/akio/lhapdf";

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag;
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}
// ----------------------------------------------------------------------------
void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{
  for ( Int_t i=0; i<n; i++ ) {
    printf("\n ==== owen: event %5d\n", i ) ; fflush(stdout) ;
    chain->Clear();
    chain->Make();
  }
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Pythia6( TString mode="pp:DY", Int_t tune=370)
{
  
  //  gSystem->Load( "libStarGeneratorPoolPythia6_4_23.so" );
  if ( LHAPDF_DATA_PATH.Contains("afs") ) {
     cout << "WARNING: LHAPDF_DATA_PATH points to an afs volume" << endl << endl;
     cout << "         You are advised to copy the PDF files you need into a local" << endl;
     cout << "         directory and set the LHAPDF_DATA_PATH to point to it."      << endl;
  }

  gSystem->Setenv("LHAPDF_DATA_PATH", LHAPDF_DATA_PATH.Data() );

  gSystem->Load( "/opt/star/$STAR_HOST_SYS/lib/libLHAPDF.so");
  gSystem->Load( "libPythia6_4_28.so");

  pythia6 = new StarPythia6("pythia6");
  pythia6->SetFrame("CMS", 510.0 );
  pythia6->SetBlue("proton");
  pythia6->SetYell("proton");
  pythia6->PyTune(tune);
  PySubs_t &pysubs = pythia6->pysubs();
  if ( mode == "pp:minbias" ) {
      pysubs.msel = 2;
  }else if( mode == "pp:qcd" ) {
      pysubs.msel = 1;
  }else if( mode == "pp:DY" ) {
      pysubs.msel = 0;
      pysubs.msub(1)=1; //qq->ll & qq->Z  
  }else if( mode == "pp:JPsi" ) {
      pysubs.msel = 0;
      pysubs.msub(86)=1;
      pysubs.msub(106)=1;
      pysubs.msub(107)=1;
  }

  pythia6->Init();
  PyPars_t &pypars = pythia6->pypars();
  printf("PARP(90) was %f, replacing it with 0.213\n",pypars.parp(90));
  pypars.parp(90)=0.213;

  // Set pi0 (id=102), pi+ (id=106) and eta (id=109) stable
  PyDat3_t &pydat3 = pythia6->pydat3();
  printf("MDCY(102,1) was %f, replacing it with 0\n",pydat3.mdcy(102,1));

  pydat3.mdcy(102,1) = 0; //PI0
  pydat3.mdcy(106,1) = 0; //PI+
  pydat3.mdcy(109,1) = 0; // ETA
  pydat3.mdcy(116,1) = 0; // K+
  pydat3.mdcy(112,1) = 0; // K_SHORT
  pydat3.mdcy(105,1) = 0; // K_LONG
  pydat3.mdcy(164,1) = 0; // ! LAMBDA0 3122
  pydat3.mdcy(167,1) = 0; // ! SIGMA0 3212
  pydat3.mdcy(162,1) = 0; // ! SIGMA- 3112
  pydat3.mdcy(169,1) = 0; // ! SIGMA+ 3222
  pydat3.mdcy(172,1) = 0; // ! Xi- 3312
  pydat3.mdcy(174,1) = 0; // ! Xi0 3322
  pydat3.mdcy(176,1) = 0; // ! OMEGA- 3334

  primaryMaker->AddGenerator(pythia6);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void runPythia( Int_t nevents=10, Int_t run=1, bool use_muon_filter=true, char* particle="mb", float vz=0.0, Int_t tune=370){ 
  cout<<"Random seed : "<<run<<endl;

 //--- owen: trying to add this to see if it will fix it (similar to fwd_tracking.C).
  gSystem->Load( "libStarRoot.so" );
  gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

  gROOT->ProcessLine(".L bfc.C");{
    TString simple = "y2017 geant gstar agml usexgeom ";
    bfc(0, simple );
  }
  
  gSystem->Load( "libVMC.so");
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );
  gSystem->Load( "StarGeneratorFilt.so");
  gSystem->Load( "libMathMore.so"   );  
  gSystem->Load( "xgeometry.so"     );


  db_maker = (St_db_Maker *)chain->GetMaker("db");
  if(db_maker){
      db_maker->SetAttr("blacklist", "tpc");
      db_maker->SetAttr("blacklist", "svt");
      db_maker->SetAttr("blacklist", "ssd");
      db_maker->SetAttr("blacklist", "ist");
      db_maker->SetAttr("blacklist", "pxl");
      db_maker->SetAttr("blacklist", "pp2pp");
      db_maker->SetAttr("blacklist", "ftpc");
      db_maker->SetAttr("blacklist", "emc");
      db_maker->SetAttr("blacklist", "eemc");
      db_maker->SetAttr("blacklist", "mtd");
      db_maker->SetAttr("blacklist", "pmd");
      db_maker->SetAttr("blacklist", "tof");
  }

  // Setup RNG seed and map all ROOT TRandom here
  StarRandom::seed( run+7000 );
  ///int rn_seed = gRandom -> Integer( 1e6 ) ;
  ///printf("\n\n ===== owen : setting random seed to %d\n\n", rn_seed ) ;
  ///StarRandom::seed( rn_seed );
  StarRandom::capture();
  
  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  //  StarPrimaryMaker *
  primaryMaker = new StarPrimaryMaker();
  {
      if ( use_muon_filter ) {
         primaryMaker -> SetFileName( Form("pythia_%s_vz%d_run%d_mufilter.root",particle,int(vz),run));
      } else {
         primaryMaker -> SetFileName( Form("pythia_%s_vz%d_run%d.root",particle,int(vz),run));
      }
      chain -> AddBefore( "geant", primaryMaker );
  }

  //
  // Setup an event generator
  //
  TString proc(particle);
  if(proc.Contains("mb")){
      Pythia6("pp:minbias", tune);
      //-- owen : add forward muon filter
      if ( use_muon_filter ) {
         printf("\n\n Adding muon filter.\n\n") ; fflush(stdout) ;
         muonfilter = (StarFilterMaker*) new FcsMuonFilter() ;
         primaryMaker->AddFilter(muonfilter) ;
         printf("\n\n Done adding muon filter.\n\n") ;
      }
  }else if(proc.Contains("jet")){  //JET
      Pythia6("pp:qcd", tune);
      jetfilter = (StarFilterMaker *)new FcsJetFilter();
      primaryMaker->AddFilter(jetfilter);
  }else if( proc.Contains("dy") && !proc.Contains("dybg") ){ // DY signal
      Pythia6("pp:DY", tune);
      dyfilter = (StarFilterMaker *)new FcsDYFilter();
      primaryMaker->AddFilter(dyfilter);
  }else if(proc.Contains("dybg") && !proc.Contains("dybgSingle") ){ //DY background via QCD
      Pythia6("pp:qcd", tune);
      dybgfilter = new FcsDYBGFilter();
      primaryMaker->AddFilter(dybgfilter);
  }else if(proc.Contains("dybgSingle")){ //DY background via QCD
      Pythia6("pp:qcd", tune);
      dybgSinglefilter = new FcsDYBGFilterSingle();
      primaryMaker->AddFilter(dybgSinglefilter);
  }else if( proc.Contains("JPsi")){ // Jpsi signal
      Pythia6("pp:JPsi", tune);
      jpsifilter = (StarFilterMaker *)new FcsJPsiFilter();
      primaryMaker->AddFilter(jpsifilter);
  }

  //
  // Setup cuts on which particles get passed to geant for
  // simulation.  
  //
  // If ptmax < ptmin indicates an infinite ptmax.
  // ptmin will always be the low pT cutoff.
  //
  //                    ptmin  ptmax
  primaryMaker->SetPtRange  (0.0,  -1.0);         // GeV
  //
  // If etamax < etamin, there is no cut in eta.
  // otherwise, particles outside of the specified range are cut.
  //
  //                    etamin etamax
  primaryMaker->SetEtaRange ( 0.0, +10.0 );
  //
  //  phirange will be mapped into 0 to 2 pi internally.
  //
  //                    phimin phimax
  primaryMaker->SetPhiRange ( 0., TMath::TwoPi() );
  
    // Setup a realistic z-vertex distribution:
  //   x = 0 gauss width = 1mm
  //   y = 0 gauss width = 1mm
  //   z = 0 gauss width = 30cm
  // 
  primaryMaker->SetVertex( 0., 0., vz );
  //primaryMaker->SetSigma( 0.1, 0.1, 30.0 );
  primaryMaker->SetSigma( 0., 0., 0. );

  //
  // Initialize primary event generator and all sub makers
  //
  primaryMaker -> Init();

  //chenging pypars(90) 
  if(tune==370){
      PyPars_t &pypars = pythia6->pypars();
      printf("PARP(90) was %f, replacing it with 0.213\n",pypars.parp(90));
      pypars.parp(90)=0.213;
  }
  //
  // Setup geometry and set starsim to use agusread for input
  //
  //geometry("fwddev1a");
  //geometry("ftsref6a");
  geometry("dev2022");
  //geometry("sitrver0");
  command("gkine -4 0");
  if ( use_muon_filter ) {
     command(Form("gfile o pythia_%s_vz%d_run%d_mufilter.fzd",particle,int(vz),run));
  } else {
     command(Form("gfile o pythia_%s_vz%d_run%d.fzd",particle,int(vz),run));
  }
  
  //
  // Trigger on nevents
  //
  trig( nevents );

  pythia6->PyStat(1);

  command("call agexit");  // Make sure that STARSIM exits properly
}
// ----------------------------------------------------------------------------

