#ifndef __LightCon_hh__
#define __LightCon_hh__

#include <unistd.h>
#include "G4LogicalVolume.hh"

class LightCon {
public:
    LightCon() {};
    ~LightCon() {}

    G4LogicalVolume* Construct_LightCon(G4String ifname);

};

#endif
