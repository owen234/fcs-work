//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable, but with only one arg


// Run very fast fwd tracking
// generate some input data using genfzd 

TFile *output = 0;

void fwd_tracking(      int n = 500,
                const char *inFile =  "simu/seed.fzd",
                std::string configFile = "simu/seed.xml",
                const char *geom = "dev2022") {
    TString _geom = geom;
    bool SiIneff = false;
    bool useConstBz = false;
    bool useFCS = true;

    // Setup the chain for reading an FZD
    TString _chain;
    if ( useFCS )
        _chain = Form("fzin %s sdt20211016 fstFastSim fcsSim fcsWFF fcsCluster fwdTrack MakeEvent StEvent ReverseField agml usexgeom bigbig  evout cmudst tree", _geom.Data());
    else 
        _chain = Form("fzin %s sdt20211016 MakeEvent StEvent ReverseField agml usexgeom bigbig fstFastSim fcsSim fwdTrack evout cmudst tree", _geom.Data());

    //_chain = Form("fzin %s sdt20230202 MakeEvent StEvent ReverseField agml usexgeom bigbig evout cmudst tree", _geom.Data());

    gSystem->Load( "libStarRoot.so" );
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);

    printf("\n\n *** Owen : turning debug on inside fzin.\n\n") ;
    St_geant_Maker* fzin_maker = (St_geant_Maker*) chain->GetMaker( "geant" ) ;
    if ( fzin_maker == 0x0 ) { printf("\n\n *** Null pointer!\n\n" ) ; return ; }
    ///////fzin_maker -> SetDebug(2) ;

    if ( useConstBz )
        StarMagField::setConstBz(true);


   gSystem->Load( "libStFttSimMaker" );
   StFttSlowSimMaker *fttSlowSim = new StFttSlowSimMaker();
    cout << "Adding StFttSlowSimMaker to chain" << endl;
   chain->AddMaker(fttSlowSim);
//goto chain_loop;

    if (useFCS) {

        StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");  
        cout << "fcsdbmkr="<<fcsdbmkr<<endl;
        StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");  
        cout << "fcsdb="<<fcsdb<<endl;    
        //fcsdbmkr->setDbAccess(1);

        // Configure FCS simulator
        StFcsFastSimulatorMaker *fcssim = (StFcsFastSimulatorMaker*) chain->GetMaker("fcsSim");
 ////// fcssim->setDebug(1);
        //fcssim->setLeakyHcal(0);


        StFcsWaveformFitMaker *fcsWFF= (StFcsWaveformFitMaker*) chain->GetMaker("StFcsWaveformFitMaker");
        fcsWFF->setEnergySelect(0);

        StFcsClusterMaker *fcsclu = (StFcsClusterMaker*) chain->GetMaker("StFcsClusterMaker");
 ////// fcsclu->setDebug(1);    
    }

    // Configure FTT FastSim
        // StFttFastSimMaker *fttFastSim = (StFttFastSimMaker*) chain->GetMaker( "fttSim" );
        // cout << "Adding StFttFastSimMaker to chain" << endl;
        // chain->AddMaker(fttFastSim);

    



    // Configure FST FastSim
        TString qaoutname(gSystem->BaseName(inFile));
        qaoutname.ReplaceAll(".fzd", ".FastSimu.QA.root");
        StFstFastSimMaker *fstFastSim = (StFstFastSimMaker*) chain->GetMaker( "fstFastSim" );;

        if (SiIneff)
            fstFastSim->SetInEfficiency(0.1); // inefficiency of Si 

        fstFastSim->SetQAFileName(qaoutname);

        cout << "Adding StFstFastSimMaker to chain" << endl;
        chain->AddMaker(fstFastSim);


    // Configure the Forward Tracker
        StFwdTrackMaker * fwdTrack = (StFwdTrackMaker*) chain->GetMaker( "fwdTrack" );;
        
        // config file set here overides chain opt
        cout << "Running FwdTracking with config: " << configFile << endl;
        fwdTrack->SetConfigFile( configFile );
        fwdTrack->SetGenerateTree( true );
        fwdTrack->SetGenerateHistograms( true );
 ////// fwdTrack->SetDebug();

        // chain->AddMaker(fwdTrack);

    if (useFCS) {
       // FwdTrack and FcsCluster assciation
       gSystem->Load("StFcsTrackMatchMaker");
       StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
	   match->setMaxDistance(6,10);
       match->setFileName("fcstrk.root");
////// match->SetDebug();
    }

    StMuDstMaker * muDstMaker = (StMuDstMaker*)chain->GetMaker( "MuDst" );
    chain->AddAfter( "FcsTrkMatch", muDstMaker );

chain_loop:
	chain->Init();

    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {

        cout << "--------->START EVENT: " << i << endl;

        chain->Clear();
        if (kStOK != chain->Make())
            break;


        // StMuDst * mds = muDstMaker->muDst();
        // StMuFwdTrackCollection * ftc = mds->muFwdTrackCollection();
        // cout << "Number of StMuFwdTracks: " << ftc->numberOfFwdTracks() << endl;
        // for ( size_t iTrack = 0; iTrack < ftc->numberOfFwdTracks(); iTrack++ ){
        //     StMuFwdTrack * muFwdTrack = ftc->getFwdTrack( iTrack );
        //     cout << "muFwdTrack->mPt = " << muFwdTrack->momentum().Pt() << endl;

        // }

        cout << "<---------- END EVENT" << endl;
    }

}
