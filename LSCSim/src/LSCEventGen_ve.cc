#include <iostream>
#include <cmath>

#include "TMath.h"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "LSCSim/LSCEventGen.hh"

#define DEBUG  0

#define Me     electron_mass_c2
#define IDe    11
#define IDve   12
#define IDvmu  14
#define IDvtau 16

using namespace std;

LSCEventGen_ve::LSCEventGen_ve(const char *arg_dbname)
    : LSCEventGen(arg_dbname)
{
    SetSeed(-1); 
    _pv0 = G4LorentzVector(0, 0, 0, 0);
    _pe0 = G4LorentzVector(0, 0, 0, 0);
    _pv1 = G4LorentzVector(0, 0, 0, 0);
    _pe1 = G4LorentzVector(0, 0, 0, 0);
}

void LSCEventGen_ve::GeneratePrimaryVertex(G4Event * argEvent)
{ ; }


void LSCEventGen_ve::GenerateEvent(double Ev, G4ThreeVector uv, double theta)
{
    double phi = -1;

    double norm = sqrt(uv.dot(uv));
    G4ThreeVector uv0 = uv;

    if (norm > 0) {
        uv0 = uv / norm;
    }

    _pv0.set(uv0*Ev, Ev);
    _pe0.set(0, 0, 0, Me);

    G4String vflavour = "e";

    if (theta == -1) {
        SetPDF(Ev, vflavour);
        G4RandGeneral rngtheta(_PDF, 721);
        theta = G4UniformRand() * pi / 2; 
    }

    if (phi == -1) {
        phi = G4UniformRand() * 2 * pi;
    }

    double Ee1 = GetEe1(Ev, theta);
    double pe1 = sqrt(Ee1*Ee1 - Me*Me);
    
    _pe1.setRThetaPhi(pe1, theta, phi);
    _pe1.setE(Ee1);
    _pe1.rotateUz(uv0);
    _pv1 = _pv0 + _pe0 - _pe1;
    _pv1.setE(sqrt(_pv1.vect()*_pv1.vect()));

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
    cerr << "#DEBUG:    _pe0\t" << _pe0 << endl;
    cerr << "#DEBUG:    _pv1\t" << _pv1 << endl;
    cerr << "#DEBUG:    _pe1\t" << _pe1 << endl;
    cerr << "#DEBUG:   total\t" << _pv1 + _pe1 - _pv0 - _pe0 << endl;
#endif

}


void LSCEventGen_ve::SetForm_HEPEvt()
{
    Form_HEPEvt tmp;
    vector<Form_HEPEvt> evt;
    int IDv = IDve;

    if      (_vflavour == "e")   IDv = IDve;
    else if (_vflavour == "mu")  IDv = IDvmu;
    else if (_vflavour == "tau") IDv = IDvtau;

    double to_GeV = 1e-3;

    // incoming neutrino
    tmp = {1, IDv, 0, 0, 
           _pv0.x()*to_GeV, _pv0.y()*to_GeV, _pv0.z()*to_GeV, // momentum in GeV
           _pv0.m()*to_GeV,               // mass in GeV
           0,                             // dt in ns
           _pos.x(), _pos.y(), _pos.z(),  // vertex _position in mm
           0, 0, 0};                      // polarization
    evt.push_back(tmp); 

    // electron
    tmp = {1, IDe, 0, 0, 
           _pe1.x()*to_GeV, _pe1.y()*to_GeV, _pe1.z()*to_GeV, // momentum in GeV
           _pe1.m()*to_GeV,               // mass in GeV
           0,                             // dt in ns
           _pos.x(), _pos.y(), _pos.z(),  // vertex _position in mm
           0, 0, 0};                      // polarization
    evt.push_back(tmp);

    // neutrino
    tmp = {1, IDv, 0, 0, 
           _pv1.x()*to_GeV, _pv1.y()*to_GeV, _pv1.z()*to_GeV, // momentum in GeV
           _pv1.m()*to_GeV,               // mass in GeV
           0,                             // dt in ns
           _pos.x(), _pos.y(), _pos.z(),  // vertex _position in mm
           0, 0, 0};                      // polarization
    evt.push_back(tmp);

    _evt = evt;
}



double LSCEventGen_ve::GetEe1(double E, double theta)
{
    double costheta = cos(theta);
        
    double T = 2 * Me * pow(E * costheta, 2);
    T /= pow((Me + E), 2) - pow(E * costheta, 2);
    
    double Ee1 = T + Me;
    
    return Ee1;
}

double LSCEventGen_ve::DiLog(double y)
{
    // source from G4LivermoreGammaConversionModelRC.cc
    y = 1-y;
    double fdilog = 0.0;
    if (y <= 0.5) {
        fdilog = pow(pi,2)/6. + (1.-y)*(log(1-y)-1.)+pow((1.-y),2)*((1./2.)*log(1.-y)-1./4.)
        +pow((1.-y),3)*((1./3.)*log(1.-y)-1./9.)+pow((1.-y),4)*((1./4.)*log(1.-y)-1./16.);
    }
    if (0.5 < y && y < 2.) {
        fdilog = 1.-y+pow((1.-y),2)/4.+pow((1.-y),3)/9.+pow((1.-y),4)/16.+
        pow((1.-y),5)/25.+pow((1.-y),6)/36.+pow((1.-y),7)/49.;
    }
    if (y >= 2.) {
        fdilog = -pow(log(y),2)/2. - pow(pi,2)/6. + (log(y)+1.)/y +
        (log(y)/2.+1./4.)/pow(y,2) + (log(y)/3.+1./9.)/pow(y,3);
    }
    return fdilog;
}

double LSCEventGen_ve::GetDifCrossSection_T(double E, double T, G4String vflavour) 
{
    _vflavour = vflavour;
    // Bahcall - Neutrino electron scattering and solar neutrino experiments

    double m_e = Me; // MeV
    static double sigma_0 = 1.7232659e-44; // cm^2 / MeV,  2 * G_F * G_F * m_e / pi
    static double SsqthW = 0.23121; // (4), PDG

    static double rho_NC_vl = 1.0126;
    double x = sqrt(1 + 2 * m_e / T);
    double I = (1 / 3. + (3. - x * x) * (0.5 * x * log((x + 1) / (x - 1)) - 1)) / 6.;
    double kappa_nue_e = 0.9791 + 0.0097 * I;
    double kappa_numu_e = 0.9970 - 0.00037 * I;

    double g_L_nue_e = rho_NC_vl * (0.5 - kappa_nue_e * SsqthW) - 1;
    double g_R_nue_e = -rho_NC_vl * kappa_nue_e * SsqthW;

    double g_L_numu_e = rho_NC_vl * (0.5 - kappa_numu_e * SsqthW);
    double g_R_numu_e = -rho_NC_vl * kappa_numu_e * SsqthW;

    double g_L, g_R;
    if (vflavour == "e")
    {
        g_L = g_L_nue_e;
        g_R = g_R_nue_e;
    }
    else if (vflavour == "mu" || vflavour == "tau")
    {
        g_L = g_L_numu_e;
        g_R = g_R_numu_e;
    }
    else
    {
        std::cerr << " Wrong vflavour value : " << vflavour << ", should be e, mu, or tau." << std::endl;
        return 0;
    }

    double q = E;
    double z = T / q;
    double E_e = T + m_e;
    double l = sqrt(E_e * E_e - m_e * m_e);
    double beta = l / E_e;
    double L_z = -DiLog(z); // L(z), L(x) = -Li_2 (x)                                                                                                                      
    double L_beta = -DiLog(beta); // L(beta)                                                                                                                               
    double Logz = log(z);

    static double alpha = 7.2973525693e-3; // (11) fine structure constant
    double term1 = E_e / l * log((E_e + l) / m_e) - 1;
    double term2 = 2 * log(1 - z - m_e / (E_e + l));

    double f_minus = term1 * (term2 - log(1 - z) - 0.5 * Logz - 5. / 12.)
              + 0.5 * (L_z - L_beta) - 0.5 * log(1 - z) * log(1 - z) - (11. / 12. + z / 2) * log(1 - z) 
              + z * (Logz + 0.5 * log(2 * q / m_e)) - (31. / 18. + Logz / 12.) * beta - 11. * z / 12. + z * z / 24.;

    double f_plus = (term1 * ((1 - z) * (1 - z) * (term2 - log(1 - z) - 0.5 * Logz - 2. / 3.) - 0.5 * (z * z * Logz + 1 - z)) 
              - 0.5 * (1 - z) * (1 - z) * (log(1 - z) * log(1 - z) + beta * (-DiLog(1 - z) - Logz * log(1 - z))) 
              + log(1 - z) * (0.5 * z * z * Logz + (1 - z) / 3. * (2 * z - 0.5)) + 0.5 * z * z * DiLog(1 - z) - z * (1 - 2 * z) / 3. * Logz - z * (1 - z) / 6.  
              - beta / 12. * (Logz + (1 - z) * ((115. - 109. * z) / 6.))) / (1 - z) / (1 - z);

    double f_plusminus = term1 * term2;

    double Lterm = g_L * g_L * (1 + alpha / pi * f_minus);
    double Rterm = g_R * g_R * (1 - z) * (1 - z) * (1 + alpha / pi * f_plus);
    double LRterm = g_L * g_R * m_e / q * z * (1 + alpha / pi * f_plusminus);

    double dsig = sigma_0 * (Lterm + Rterm - LRterm);
    double T_max = E / (1 + m_e / 2 / E);

    if (T > T_max) dsig = 0;
    if (dsig < 0) dsig = 0;

    return dsig;
}

double LSCEventGen_ve::GetDifCrossSection_costheta(double E, double costheta, G4String vflavour)
{
    double E0 = Me + E;
    double Ecostheta2 = pow(E * costheta, 2);
    double T = 2 * Me * Ecostheta2 / ( E0*E0 - Ecostheta2 );
    double conv = 4 * E0*E0 * Ecostheta2 / pow(E0*E0 - Ecostheta2, 2);

    return conv * GetDifCrossSection_T(E, T, vflavour);
}

double LSCEventGen_ve::GetDifCrossSection(double E, double theta, G4String vflavour)
{
    double costheta = cos(theta);

    return GetDifCrossSection_costheta(E, costheta, vflavour);
}


void LSCEventGen_ve::SetPDF(double E, G4String vflavour)
{
    double theta = 0;
    double dtheta = pi / 720;

    for (int i=0; i<721; i++) {
        theta = i*dtheta;
        _PDF[i] = GetDifCrossSection(E, theta, vflavour)*1e45;
    }
}

