//

#ifndef __LSCEvtGen_hh__
#define __LSCEvtGen_hh__ 

#include "GLG4VertexGen.hh"

class LSCEventGen : public GLG4VVertexGen {
public:
    virtual LSCEventGen(const char * arg_dbname) : _dbname(arg_dbname) {}
    virtual ~LSCEventGen() {}
    
    virtual void GeneratePrimaryVertex(G4Event * argEvent) = 0;
    virtual void SetState(G4String newValues) = 0;
    virtual G4String GetState() = 0;

    void Print_HEPEvt();

    // set functions
    void SetSeed(double seed) { _rseed = seed; }

    // get functions
    double GetSeed() const { return _rseed; }

private:
    double _rseed = 0;
};


class LSCEventGen_IBD : public LSCEventGen {
public:
    virtual LSCEventGen(const char * arg_dbname);
    virtual ~LSCEventGen();
    virtual void GeneratePrimaryVertex(G4Event * argEvent);
    virtual void SetState(G4String newValues); 

    void GenerateIBD(double Ev);
private:
    G4PrimaryParticle* _v0;  // incoming neutrino
    G4PrimaryParticle* _p0;  // proton
    G4PrimaryParticle* _e1;  // positron
    G4PrimaryParticle* _n1;  // neutron
};


class LSCEventGen_ve : public LSCEventGen {
public:
    virtual LSCEventGen(const char * arg_dbname);
    virtual ~LSCEventGen();
    virtual void GeneratePrimaryVertex(G4Event * argEvent);
    virtual void SetState(G4String newValues);

private:
    G4PrimaryParticle* _v0; // incoming neutrino
    G4PrimaryParticle* _e0; // electron
    G4PrimaryParticle* _v1; // neutrino
    G4PrimaryParticle* _e1; // electron
};


