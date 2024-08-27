#include <iostream>
#include <cmath>

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "LSCSim/LSCEventGen.hh"
#include "GLG4Sim/GLG4PosGen.hh"

#define DEBUG   0

#define Mp      proton_mass_c2 
#define Mn      neutron_mass_c2
#define Me      electron_mass_c2
#define M       ((Mp + Mn) / 2)
#define Delta   (Mn - Mp)

using namespace std;

LSCEventGen_IBD::LSCEventGen_IBD(const char *arg_dbname)
    : LSCEventGen(arg_dbname)
{
    SetSeed(-1);
    _pv0 = G4LorentzVector(0, 0, 0, 0);
    _pp0 = G4LorentzVector(0, 0, 0, 0);
    _pe1 = G4LorentzVector(0, 0, 0, 0);
    _pn1 = G4LorentzVector(0, 0, 0, 0);
}


void LSCEventGen_IBD::GeneratePrimaryVertex(G4Event * argEvent)
{ return; }


void LSCEventGen_IBD::GenerateEvent(double Ev, G4ThreeVector uv, double theta)
{
    double phi = -1;

    G4ThreeVector uv0 = uv; // normalize
    double norm = sqrt(uv.dot(uv));

    if (norm > 0) {
        uv0 /= norm;
    }
    else {
        uv0 = _pos - _pos_source;
        norm = sqrt(uv0.dot(uv0));
        uv0 /= norm;
    }    

    _pv0.set(uv0*Ev, Ev);
    _pp0.set(0, 0, 0, Mp);


    if (theta == -1) {
        SetPDF(Ev);
        G4RandGeneral rngtheta(_PDF, 721);
        theta = rngtheta() * pi;
    }

    if (phi == -1) {
        phi = G4UniformRand() * 2 * pi;
    }

    double Ee1 = GetEe1(Ev, theta);
    double pe1 = sqrt(Ee1*Ee1 - Me*Me);

    _pe1.setRThetaPhi(pe1, theta, phi);
    _pe1.setE(Ee1);
    _pe1.rotateUz(uv0);
    _pn1 = _pv0 + _pp0 - _pe1;

    _pn1.setE(sqrt(_pn1.vect() * _pn1.vect() + Mn * Mn));
            
#if DEBUG > 1
    cerr << "#DEBUG: Positron momentum calculation" << endl;
    cerr << "#DEBUG:   theta\t" << theta << endl;
    cerr << "#DEBUG:     phi\t" << phi << endl;
    cerr << "#DEBUG:     pe1\t" << pe1 << endl;
    cerr << "#DEBUG:     Ee1\t" << Ee1 << endl;
    cerr << "#DEBUG:    _pe1\t" << _pe1 << endl;
#endif

#if DEBUG 
    cerr << "#DEBUG: Four momentum vectors" << endl;
    cerr << "#DEBUG:    _pv0\t" << _pv0 << endl;
    cerr << "#DEBUG:    _pp0\t" << _pp0 << endl;
    cerr << "#DEBUG:    _pe1\t" << _pe1 << endl;
    cerr << "#DEBUG:    _pn1\t" << _pn1 << endl;
    cerr << "#DEBUG:   total\t" << _pn1 + _pe1 - _pv0 - _pp0 << endl;
#endif

    return;
} 

void LSCEventGen_IBD::SetForm_HEPEvt()
{
    // neutrino
    Form_HEPEvt tmp;
    vector<Form_HEPEvt> evt;
    double to_GeV = 1e-3;

    tmp = {1, -12, 0, 0, 
           _pv0.x()*to_GeV, _pv0.y()*to_GeV, _pv0.z()*to_GeV, // momentum in GeV
           _pv0.m()*to_GeV,               // mass in GeV
           0,                             // dt in ns
           _pos.x(), _pos.y(), _pos.z(),  // vertex _position in mm
           0, 0, 0};                      // polarization
    evt.push_back(tmp); 

    // _positron
    tmp = {1, -11, 0, 0, 
           _pe1.x()*to_GeV, _pe1.y()*to_GeV, _pe1.z()*to_GeV, //momentum in GeV
           _pe1.m()*to_GeV,               // mass in GeV
           0,                             // dt in ns
           _pos.x(), _pos.y(), _pos.z(),  // vertex _position in mm
           0, 0, 0};                      // polarization
    evt.push_back(tmp);

    // neutron
    tmp = {1, 2112, 0, 0, 
           _pn1.x()*to_GeV, _pn1.y()*to_GeV, _pn1.z()*to_GeV, //momentum in GeV
           _pn1.m()*to_GeV,               // mass in GeV
           0,                             // dt in ns
           _pos.x(), _pos.y(), _pos.z(),  // vertex _position in mm
           0, 0, 0};                      // polarization
    evt.push_back(tmp);

    _evt = evt;
}


double LSCEventGen_IBD::GetEe1(double E, double theta)
{
    /*-----------------------------------------------------------------
    Computes energy transfer from neutrino to positron 
    masses and energies are in MeV.
    E is neutrino energy and theta is scattered angle of the positron. 
    -----------------------------------------------------------------*/

    double costheta = cos(theta);
    double Ee0      = E - Delta;

    if (Ee0 < Me) { return 0; }

    double pe0 = sqrt(Ee0*Ee0 - Me*Me);
    double ve0 = pe0 / Ee0;
    double d1  = E/M * (1 - ve0 * costheta);
    double d2  = (Delta*Delta - Me*Me) / (2*M);
    double Ee1 = Ee0 * (1 - d1) - d2;

    if (Ee1*Ee1 - Me*Me < 0) { return 0; }
    else { return Ee1; }
}


double LSCEventGen_IBD::GetDifCrossSection(double E, double theta)
{
    return GetDifCrossSection_costheta(E, cos(theta));
}

double LSCEventGen_IBD::GetDifCrossSection_costheta(double E, double costheta)
{
    /*-----------------------------------------------------------------
    Computes differential cross section 
    for the given neutrino energy, E,
    and the positron scattering angle, theta.
    Vogel P., Beacom J. F. (1999) Phys Rev D 60 53003 
    -----------------------------------------------------------------*/

    double Ee0   = E - Delta;

    if (Ee0 < Me) { return 0; }

    double pe0   = sqrt(Ee0*Ee0 - Me*Me);
    double ve0   = pe0/Ee0;

    double f     = 1.0;
    double g     = 1.26;
    double f2    = 3.706;

    double c1    = 2*(f + f2)*g;
    double c2    = f*f + g*g;
    double c3    = f*f + 3*g*g;
    double c4    = f*f - g*g;

    double d1    = E/M * (1 - ve0*costheta);
    double d2    = (Delta*Delta - Me*Me)/2/M;
    double Ee1   = Ee0*(1 - d1) - d2;

    if (Ee1*Ee1 - Me*Me < 0) { return 0; }

    double pe1   = sqrt(Ee1*Ee1 - Me*Me);
    double ve1   = pe1/Ee1;

    double gamma; 
    gamma  = c1*((2*Ee0 + Delta)*(1 - ve0*costheta) - Me*Me/Ee0);
    gamma += c2*(Delta*(1 + ve0*costheta) + Me*Me/Ee0);
    gamma += c3*((Ee0 + Delta)*(1 - ve0*costheta) - Delta);
    gamma += c4*((Ee0 + Delta)*(1 - ve0*costheta) - Delta) * ve0*costheta;

    double GF = 1.16637e-11;
    double Vud = 0.97373;
    double delta_inner = 0.024;
    double hbarc = 197.327053;

    double sig_0 = GF*GF * Vud*Vud / pi * (1. + delta_inner) * hbarc * hbarc * 1e-26;
    double dsig = (c3 + c4*ve1*costheta)*Ee1*pe1 - gamma*Ee0*pe0/M;                         
    dsig *= sig_0 / 2;

    return dsig;
}

void LSCEventGen_IBD::SetPDF(double E)
{
    double theta = 0;
    double dtheta = pi / 720;

    for (int i=0; i<721; i++) {
        theta = i * dtheta; 
        _PDF[i] = GetDifCrossSection(E, theta)*1e45;
    }
}

int LSCEventGen_IBD::GeneratePosition_sourceMC(double L)
{
    double Lmin = GetLmin();

    double L1 = L - Lmin;
    
    double factor = Lmin * Lmin / L / L;

#if DEBUG > 2
    cerr << "Lmin " << Lmin << endl;
    cerr << "L " << L << endl;
    cerr << "factor "  << factor << endl;
#endif

    if (G4UniformRand() < factor) 
        return 1;
    else 
        return 0;
}

double LSCEventGen_IBD::GetLmin()
{
    double Lmin = 0;
    double targetR = _geom_db["target_radius"];
    double targetH = _geom_db["target_height"];
    double targetT = _geom_db["target_thickness"];

    targetR *= 10;
    targetH *= 10;
    targetT *= 10;

    double LS_R = targetR - targetT;
    double LS_Hmin = -(targetH/2 - targetT);
    double LS_Hmax = targetH/2 - targetT;

    G4ThreeVector pos_src = GetPositionSource();

    double src_R = sqrt(pow(pos_src.x(), 2) + pow(pos_src.y(), 2));
    double src_Z = pos_src.z();

#if DEBUG > 2
    cerr << "src_R = " << src_R << endl;
    cerr << "src_Z = " << src_Z << endl;
#endif

    if (src_R > LS_R) {
        if (src_Z >= LS_Hmin && src_Z <= LS_Hmax)
            Lmin = src_R - LS_R;
        else if(src_Z < LS_Hmin) 
            Lmin = sqrt(pow(src_R - LS_R, 2) + pow(src_Z - LS_Hmin, 2));
        else if(src_Z > LS_Hmax) 
            Lmin = sqrt(pow(src_R - LS_R, 2) + pow(LS_Hmax - src_Z, 2));
    }
    else {
        if (src_Z >= LS_Hmin && src_Z <= LS_Hmax)
            Lmin = 0;
        else if(src_Z < LS_Hmin) 
            Lmin = sqrt(pow(src_Z - LS_Hmin, 2));
        else if(src_Z > LS_Hmax) 
            Lmin = sqrt(pow(LS_Hmax - src_Z, 2));
    }

    return Lmin; 
}
