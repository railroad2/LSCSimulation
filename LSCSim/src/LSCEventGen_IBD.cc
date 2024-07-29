#include <iostream>

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "LSCEvtGen/LSCEvtGen.hh"
#include "GLG4Sim/GLG4PosGen.hh"

using namespace std;

LSCEventGen_IBD::LSCEventGen_IBD(const char *arg_dbname)
    : GLG4VVertexGen(arg_dbname), 
      _v0(-12),
      _p0(2212),
      _e1(-11),
      _n1(2112)
{}
      

void LSCEventGen_IBD::GeneratePrimaryVertex(G4Event * argEvent)
{
}

void LSCEventGen_IBD::GenerateIBD(double Ev, G4ThreeVector uv, double theta=-1) 
{
    _v0.Set4Momentum(uv.x(), uv.y(), uv.z(), Ev); // xyzt
    if (theta == -1) {
        theta = G4UniformRand() * pi; 
    }

    cout << _v0 << endl;
    cout << theta << endl;
}

