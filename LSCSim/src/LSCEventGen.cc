#include "LSCSim/LSCEventGen.hh"

#include <iostream>
#include <sys/time.h>
#include <ctime>

using namespace std;

void LSCEventGen::SetSeed(long seed=-1)
{
    if (seed < 0) {
        struct timeval tv; 
        gettimeofday(&tv, NULL);
        seed = tv.tv_sec*1000 + tv.tv_usec;
    }

    _rseed = seed;
    G4Random::setTheSeed(_rseed);
}

void LSCEventGen::Print_HEPEvt()
{
#if DEBUG
    cerr << "#DEBUG: Print_HEPEvt()" << endl;
#endif

    int print_pol = 0;

    if (_evt[0].PLX == 0 && _evt[i].PLY == 0 && _evt[0].PLZ == 0) 
        print_pol = 0;
    else
        print_pol = 1;

    cout << _evt.size() << endl; // number of primary particles

    for (int i=0; i<_evt.size(); i++) {
        cout << _evt[i].ISTHEP  << '\t';
        cout << _evt[i].IDHEP   << '\t';
        cout << _evt[i].JDAHEP1 << '\t';
        cout << _evt[i].JDAHEP2 << '\t';
        cout << scientific;
        cout << _evt[i].PHEP1   << '\t';
        cout << _evt[i].PHEP2   << '\t';
        cout << _evt[i].PHEP3   << '\t';
        cout << _evt[i].PHEP5   << '\t';
        cout << _evt[i].DT      << '\t';
        cout << _evt[i].X       << '\t';
        cout << _evt[i].Y       << '\t';
        cout << _evt[i].Z       << '\t';
        if (print_pol) {
            cout << _evt[i].PLX << '\t';
            cout << _evt[i].PLY << '\t';
            cout << _evt[i].PLZ << '\t';
        }
        cout << endl;
    }

    return;
}

