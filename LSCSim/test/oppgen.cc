#include <iostream>
#include <sstream>
#include <getopt.h>

#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCEventGen.hh"

using namespace std;

G4String detshape = "sphere";

G4ThreeVector string2dir(string input);
int opticalphotons(int nevent, double nphotons, double Ev, G4ThreeVector pos);
bool exists(const G4String &fname);
void PrintHelp();

int main(int argc, char** argv)
{
    if (argc < 2) {
        PrintHelp();
        return -1;
    }

    int opt;

    int nevent = 1;
    double nphotons = 1;  // number of photons per event
    double Eop = -1;       // optical photon energy

    G4String evt_type = "reactor";
    G4String fn_geom;
    G4ThreeVector pos(0, 0, 0);

    while ((opt = getopt(argc, argv, "n:g:e:d:k:p:h")) != -1) {
        switch (opt) {
            case 'n': nevent = atoi(optarg); break;
            case 'g': fn_geom = G4String(optarg); break;
            case 'e': Eop = atof(optarg); break;
            case 'k': nphotons = atof(optarg); break;
            case 'd': detshape= G4String(optarg); break;
            case 'p': pos = string2dir(optarg); break;
            case 'h': PrintHelp(); break;
            default: PrintHelp();
        }
    }

    opticalphotons(nevent, nphotons, Eop, pos);

    return 0;
}

int opticalphotons(int nevent, double nphotons, double Ev, G4ThreeVector pos)
{
    LSCEventGen_op* evtgen = new LSCEventGen_op();
    evtgen->SetPosition(pos);
    
    bool rndEv = false;
    if (Ev == -1) rndEv = true;

    /*
    if (!exists(fn_geom)) {
        cout << "File does not exist: " << fn_geom<< endl;
        return -1;
    }
    evtgen->ReadGeometryFile(fn_geom);

    GLG4param& src_db = GLG4param::GetDB();

    pos_src = G4ThreeVector(src_db["x"], src_db["y"], src_db["z"]);
    */

    for (int i=0; i<nevent; i++) {
        if (rndEv) Ev = 3e-6; // TODO make random energy generator for optical photon in LS

        evtgen->GenerateEvent(nphotons, Ev);
        evtgen->Print_HEPEvt();
    }

    return 0;
}

G4ThreeVector string2dir(string input) 
{
    char separator = ',';    
    istringstream iss(input);
    string str_buf;
    vector<double> vec;
    while (getline(iss, str_buf, separator)) {
        vec.push_back(stof(str_buf));
    }

    return G4ThreeVector(vec[0], vec[1], vec[2]);
}


void PrintHelp()
{
    cout << endl;
    cout << "Usage: oppgen [-n # of event] [-k # of photons per event]" << endl;
    cout << "              [-e photon energy] " << endl;
}

bool exists(const G4String &fname)
{
    ifstream f(fname);
    return f.good();
}
