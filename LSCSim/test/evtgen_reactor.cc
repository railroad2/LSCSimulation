#include <iostream>

#include "LSCSim/LSCEventGen.hh"

using namespace std;

int doit(int nevent)
{
    LSCEventGen_IBD* evtgen = new LSCEventGen_IBD();

    double Ev = 5.0;
    double theta = -1;
    G4ThreeVector uv(1, 0, 0);
    
    G4String fn_geom = "/home/kmlee/opt/lscsim/LSCSimulation/LSCSim/data/geometry.dat";
    evtgen->ReadGeometryFile(fn_geom);

    for (int i=0; i<nevent; i++) {
        evtgen->GeneratePosition();
        evtgen->GenerateEvent(Ev, uv, theta); 
        evtgen->SetFormat_HEPEvt();
        evtgen->Print_HEPEvt();
    }

    return 0;
}

int main(int argc, char** argv)
{
    int nevent = 1;
    if (argc > 1) {
        nevent = atoi(argv[1]);
    }
    doit(nevent); 

    return 0;
}
