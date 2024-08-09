#include <iostream>
#include <getopt.h>

#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCEventGen.hh"

using namespace std;

G4String fn_geom = "/home/kmlee/opt/lscsim/LSCSimulation/LSCSim/data/geometry.dat";
G4String fn_source = "/home/kmlee/opt/lscsim/LSCSimulation/LSCSim/data/source.dat";

int reactor(int nevent, double Ev);
int source(int nevent, double Ev);
int solar(int nevent, double Ev);
void PrintHelp();

int main(int argc, char** argv)
{
    if (argc < 2) {
        PrintHelp();
        return -1;
    }

    int opt;

    int nevent = 1;
    double Ev = -1;  // neutrino energy

    G4String evt_type = "reactor";
    //G4String fn_geom = "geometry.dat";
    G4String fn_setting;

    while ((opt = getopt(argc, argv, "t:n:g:e:s:h")) != -1) {
        switch (opt) {
            case 't': evt_type = G4String(optarg); break;
            case 'n': nevent = atoi(optarg); break;
            case 'g': fn_geom = G4String(optarg); break;
            case 'e': Ev = atof(optarg); break;
            case 's': fn_source = G4String(optarg); break;
            case 'h': PrintHelp(); break;
            default: PrintHelp();
        }
    }

    if      (evt_type == "reactor") reactor(nevent, Ev);
    else if (evt_type == "source")  source (nevent, Ev);
    else if (evt_type == "solar")   solar  (nevent, Ev);
    else PrintHelp();

    return 0;
}

int reactor(int nevent, double Ev)
{
    LSCEventGen_IBD* evtgen = new LSCEventGen_IBD();
    
    bool rndEv = false;
    if (Ev == -1) rndEv = true;

    double theta = -1;
    G4ThreeVector uv(1, 0, 0); // direction of incoming neutrino
    
    evtgen->ReadGeometryFile(fn_geom);

    for (int i=0; i<nevent; i++) {
        if (rndEv) Ev = G4UniformRand() * 9.2 + 1.8; // reactor neutrino 1.8 ~ 10 MeV

        evtgen->GeneratePosition();
        evtgen->GenerateEvent(Ev, uv, theta); 
        evtgen->Print_HEPEvt();
    }

    return 0;
}


int source(int nevent, double Ev)
{
    LSCEventGen_IBD* evtgen = new LSCEventGen_IBD();

    bool rndEv = false;
    if (Ev == -1) rndEv = true;

    double theta = -1;
    G4ThreeVector uv(1, 0, 0); 
    G4ThreeVector pos_src;
    GLG4param& src_db = GLG4param::GetDB();

    evtgen->ReadGeometryFile(fn_geom);
    src_db.ReadFile(fn_source);

    double emin = src_db["emin"];
    double emax = src_db["emax"];

    pos_src = G4ThreeVector(src_db["sourcex"], src_db["sourcey"], src_db["sourcez"]);
    pos_src *= 10; // cm -> mm;

    for (int i=0; i<nevent; i++) {
        if (rndEv) Ev = G4UniformRand() * (emax - emin) + emin; // source neutrino 1.8 ~ 10 MeV

        evtgen->GeneratePosition();
        uv = evtgen->GetPosition() - pos_src;
        evtgen->GenerateEvent(Ev, uv, theta); 
        evtgen->Print_HEPEvt();
    }

    return 0;
}


int solar(int nevent, double Ev)
{
    LSCEventGen_ve* evtgen = new LSCEventGen_ve();

    bool rndEv = false;
    if (Ev == -1) rndEv = true;

    double theta = -1;
    G4ThreeVector uv(0, 0, 1);
    
    evtgen->ReadGeometryFile(fn_geom);

    for (int i=0; i<nevent; i++) {
        if (rndEv) Ev = G4UniformRand() * 20 + 0.01; // solar neutrino 0.01 ~ 20.01 MeV

        evtgen->GeneratePosition();
        evtgen->GenerateEvent(Ev, uv, theta); 
        evtgen->Print_HEPEvt();
    }

    delete (evtgen);

    return 0;
}


void PrintHelp()
{
    cout << endl;
    cout << "Usage: evtgen [-t event type] [-n # of event] [-g geometryfile]" << endl;
    cout << "              [-e neutrino energy] [-s sourcefile]" << endl;
    //cout << "              [-e neutrino energy] <setting file>" << endl;
    cout << "       event type - reactor, source, solar " << endl;
}

