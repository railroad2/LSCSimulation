#include <iostream>
#include "TMath.h"

#include "LSCSim/LSCEventGen.hh"

using namespace std;

int doit()
{
    LSCEventGen_ve* evtgen = new LSCEventGen_ve("");

    double Ev = 2.0;
    double theta = -1;
    G4ThreeVector uv(0, 0, 1);
    
    G4String fn_geom = "/home/kmlee/opt/lscsim/LSCSimulation/LSCSim/data/geometry.dat";
    evtgen->ReadGeometryFile(fn_geom);

    evtgen->GeneratePosition();
    evtgen->GenerateEvent(Ev, uv, theta); 
    evtgen->SetFormat_HEPEvt();
    evtgen->Print_HEPEvt();

    delete (evtgen);

    return 0;
}

int test_dilog()
{
    LSCEventGen_ve evtgen("");
    double x = 0;;
    for (int i=0; i<30; i++) {
        x = -1.5+i*0.1;
        double a = TMath::DiLog(x);
        double b = evtgen.DiLog(x);
        cout << "--------------------------------" << endl;
        cout << "x\t" << x << endl;
        cout << "ROOT::TMath::DiLog(x)\t" << a << endl;
        cout << "LSCEventgen_ve::DiLog(x)\t"<< b << endl;
        cout << "Difference\t" << a - b << endl;
        cout << "--------------------------------" << endl;
    }
    return 0;
}

int main(int argc, char** argv)
{
    doit(); 
    //test_dilog();

    return 0;
}

