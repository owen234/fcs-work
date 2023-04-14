/*! \class FcsMuonFilter

  This class is just used to filter non forward muon events
  Mostly copied from FcsJPsiFilter.h
*/

#ifndef STAR_FcsMuonFilter
#define STAR_FcsMuonFilter

#include <vector>
#include <string>
#include "StarGenerator/FILT/StarFilterMaker.h"
#include "StarGenerator/EVENT/StarGenEvent.h"

class StarGenParticleMaster;
class StarGenParticle;
class StarGenEvent;

class FcsMuonFilter : public StarFilterMaker
{
 public:
  FcsMuonFilter(); ///constructor
  virtual ~FcsMuonFilter(){;};///destructor
  Int_t Filter( StarGenEvent *mEvent );

  ClassDef(FcsMuonFilter,1);  
};

#endif
