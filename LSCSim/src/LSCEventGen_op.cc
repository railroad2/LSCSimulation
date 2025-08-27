#include <iostream>

#include "G4RandomTools.hh"

#include "LSCSim/LSCEventGen.hh"

#define DEBUG  0

#define IDop   0

using namespace std;

LSCEventGen_op::LSCEventGen_op(const char *arg_dbname)
    : LSCEventGen(arg_dbname)
{
    SetSeed(-1); 
}


void LSCEventGen_op::GenerateEvent(int nphotons, double E)
{
    double phi = -1;

    Form_HEPEvt tmp;
    vector<Form_HEPEvt> evt;

    double to_GeV = 1e-3;
    G4ThreeVector dir; 
    G4ThreeVector mom; 

    for (int i=0; i<nphotons; i++) {
        dir = G4RandomDirection();
        mom = dir * E;
        // optical photon
        tmp = {1, IDop, 0, 0, 
               mom.x()*to_GeV, mom.y()*to_GeV, mom.z()*to_GeV, // momentum in GeV
               E*to_GeV,                      // mass in GeV
               0,                             // dt in ns
               _pos.x(), _pos.y(), _pos.z(),  // vertex _position in mm
               0, 0, 0};                      // polarization
        evt.push_back(tmp); 
    }

    _evt = evt;
}



