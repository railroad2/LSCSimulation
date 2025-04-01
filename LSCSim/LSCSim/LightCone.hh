#ifndef __LightCone_hh__
#define __LightCone_hh__

#include "G4LogicalVolume.hh"

class LightCone : public G4LogicalVolume {
public:
    LightCone(G4String ifname); 
    ~LightCone() {}

    void Construct_LightCone(G4String ifname);

}

#endif
