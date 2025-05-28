#ifndef MCPhotonHit_hh
#define MCPhotonHit_hh

#include "TObject.h"

class MCPhotonHit : public TObject {
public:
  MCPhotonHit();
  MCPhotonHit(const MCPhotonHit & photon);
  ~MCPhotonHit();

  void SetTime(float t);
  void SetKineticEnergy(float KE);
  void SetProcessName(const std::string processName) { fProcessName = processName; } // added by kmlee

  float GetTime() const;
  float GetKineticEnergy() const;
  float GetWavelength() const;
  std::string GetProcessName() const; // added by kmlee

private:
  float fTime;
  float fKE;

  std::string fProcessName; // Name of the process that created this hit, added by kmlee

  ClassDef(MCPhotonHit, 1)
};

inline void MCPhotonHit::SetTime(float t) { fTime = t; }
inline void MCPhotonHit::SetKineticEnergy(float KE) { fKE = KE; }
inline float MCPhotonHit::GetTime() const { return fTime; }
inline float MCPhotonHit::GetKineticEnergy() const { return fKE; }
inline float MCPhotonHit::GetWavelength() const { return 1239.83968E-6 / fKE; }
inline std::string MCPhotonHit::GetProcessName() const { return fProcessName; } // added by kmlee

#endif
