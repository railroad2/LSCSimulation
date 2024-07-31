#include <iostream>

#include "LSCSim/LSCEventGen.hh"

using namespace std;

int doit()
{
    LSCEventGen_IBD* evtgen = new LSCEventGen_IBD("");

    double Ev = 2.0;
    double theta = -1;
    G4ThreeVector uv(1, 0, 0);
    
    evtgen->GenerateEvent(Ev, uv, theta); 
    evtgen->SetFormat_HEPEvt();
    evtgen->Print_HEPEvt();

    return 0;
}

int main(int argc, char** argv)
{
    doit(); 

    return 0;
}
