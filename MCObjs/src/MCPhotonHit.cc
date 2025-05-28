#include "MCObjs/MCPhotonHit.hh"

ClassImp(MCPhotonHit)

MCPhotonHit::MCPhotonHit()
  : TObject()
{
  fTime = 0;
  fKE = 0;
  fProcessName = ""; // added by kmlee
}

MCPhotonHit::MCPhotonHit(const MCPhotonHit & photon)
  : TObject()
{
  fTime = photon.GetTime();
  fKE = photon.GetKineticEnergy();
  fProcessName = photon.GetProcessName(); // added by kmlee
}

MCPhotonHit::~MCPhotonHit() {}
