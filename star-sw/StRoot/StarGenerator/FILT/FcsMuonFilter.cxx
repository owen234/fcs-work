#include "FcsMuonFilter.h"

#include "StarGenerator/EVENT/StarGenParticle.h"
#include "StarGenerator/EVENT/StarGenEvent.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

static const float ZFCS    = 710.16 + 13.9 + 15.0;
static const float XFCSMin = 16.69 + 1.5*5.542;
static const float XFCSMax = 16.69 + 22.0 * 5.542 - 1.5*5.542;
static const float YFCS    = 34.0/2.0 * 5.542 - 1.5*5.542;
static const float ETCUT   = 1.0;
static const float MASSCUT = 2.5;

static int ntot=0;
static int ngood=0;

using namespace std;
//_______________________________________________________________
FcsMuonFilter::FcsMuonFilter():StarFilterMaker("fcsMuonFilter"){
    cout<<"FCS Muon filter is used!!!"<<endl;
}
//_______________________________________________________________
Int_t FcsMuonFilter::Filter( StarGenEvent *mEvent){
    ntot++;
    // Get a reference to the current event 
    StarGenEvent& event = *mEvent;
    
    //event.Print();
    if(event.GetNumberOfParticles() <= 0) {return kError;}
    
    // apply Muon conditions for events
    TIter Iterator = event.IterAll();
    StarGenParticle *p = 0;
    vector<StarGenParticle*> forwardParticles;
    if ( ntot%1000 == 0 ) printf("  FcsMuonFilter::Filter: ntot %9d  ngood %9d\n", ntot, ngood ) ; fflush(stdout) ;
    while( ( p = (StarGenParticle*)Iterator.Next() ) ){
    // if ( ntot%1000 == 0 ) {
    //    printf("      pid = %8d status=%2d  px,py,pz = (%7.2f, %7.2f, %7.2f)\n", p->GetId(), p->GetStatus(), p->GetPx(), p->GetPy(), p->GetPz() ) ;
    // }
       if(p->GetStatus() != 1)continue;
       int pid = abs(p->GetId());
       if(pid!=13) continue; //  muons only
       float pt = p->pt() ;
       float theta = atan2( pt, p->GetPz() ) ;
       if ( theta == 0 ) theta = 1e-6 ;
       float eta = -1 * log( tan( theta/2. ) ) ;
       printf("  FcsMuonFilter::Filter: ntot = %9d, ngood = %9d : found a muon  pid = %8d status=%2d  px,py,pz = (%7.2f, %7.2f, %7.2f)   pt = %7.2f, eta = %7.3f\n", ntot, ngood, p->GetId(), p->GetStatus(), p->GetPx(), p->GetPy(), p->GetPz(), pt, eta ) ; fflush(stdout) ;

       if ( !(eta > 2.2 && eta < 4.5 && pt > 0.2 ) ) continue ;

    // if(p->GetPz()<0.0) continue; // +z direction only
    // if(p->pt()<ETCUT) continue;
    // //simple box cut around FCS  (vertex here is in [mm])
    // float x = fabs (p->GetVx()/10.0 + p->GetPx() / p->GetPz() * (ZFCS - p->GetVz()/10.0));
    // if(x<XFCSMin || x>XFCSMax) continue;
    // float y = fabs (p->GetVy()/10.0 + p->GetPy() / p->GetPz() * (ZFCS - p->GetVz()/10.0));
    // if(y>YFCS) continue;	

       forwardParticles.push_back(p);
    }
    unsigned int size=forwardParticles.size();
    if(size<1) return StarGenEvent::kReject;
    printf("  FcsMuonFilter : accept event\n") ; fflush(stdout) ;
    ngood ++ ;
    return (StarGenEvent::kAccept | 0x08);
}
