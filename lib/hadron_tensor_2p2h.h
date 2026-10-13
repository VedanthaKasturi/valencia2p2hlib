//
//  hadron_tensor_2p2h.h
//  Pions
//
//  Created by Joanna Sobczyk on 20/01/16.
//  Copyright (c) 2016 Joanna Sobczyk. All rights reserved.
//

#ifndef __Pions__hadron_tensor_2p2h__
#define __Pions__hadron_tensor_2p2h__

#include <stdio.h>
#include "1p1h1pion.h"
#include "utils.h"

#include "event.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"




extern bool FERMI;
extern int NR_PART;
enum CHANNEL {CC, CCanti};
extern CHANNEL CH_2p2h;
extern bool Global;



void ulind(double QZR,double Q, double XKF,std::complex<double> &CUFUN, std::complex<double> &cudel);
double DUR (double q0, double dq, double pf1, double pf2);
std::complex<double> delta_lind (std::complex<double> q_zero,double q_mod, double rho,double k_fermi);
std::complex<double> ulind_crossed(double QZR, double Q, double XKF);
void delta_lind_pieces(std::complex<double> q_zero,double q_mod, double rho,double k_fermi, std::complex<double>& direct, std::complex<double>& crossed);



/**********************
 Isospin bookkeeping for asymmetric nuclei (Z != N)

 Every diagram is a W vertex on a nucleon hole N1 that emits a pion, which is absorbed on a second nucleon
 N2 (a particle-hole excitation). A nucleon-hole bubble hole -> part has its own Fermi momenta: kF of the hole
 species fixes the range of the hole energy, kF of the particle species the Pauli blocking of the final nucleon.
 The weights are labelled by the final-state charges (N1', N2'): for neutrinos pp, np, pn and, as the p<->n
 mirror, nn, pn, np for antineutrinos. The code works with the neutrino isospin tables for both; for
 antineutrinos the species below are mirrored.
 **********************/
enum Nucleon {PROT, NEUT};

inline Nucleon mirror(Nucleon n) { return n==PROT ? NEUT : PROT; }
inline double kF_of(const FermiGas& fg, Nucleon n) { return n==PROT ? fg.kF_p : fg.kF_n; }

struct BubbleSpecies {
    Nucleon hole;
    Nucleon part;
};

inline BubbleSpecies bubble_species(Nucleon hole, Nucleon part) {
    BubbleSpecies b;
    b.hole = CH_2p2h==CCanti ? mirror(hole) : hole;
    b.part = CH_2p2h==CCanti ? mirror(part) : part;
    return b;
}

// Struck-nucleon bubble (N1) of the isospin slots W_1, W_2, W_3 (A[0], A[1], A[2]) of the one-body diagrams.
// slot 0: p -> p pi+ (ch 1);  slot 1: n -> n pi+ (ch 2);  slot 2: n -> p pi0 (ch 3)
inline BubbleSpecies vertex_bubble(int slot) {
    switch (slot) {
        case 0: return bubble_species(PROT, PROT);
        case 1: return bubble_species(NEUT, NEUT);
        default: return bubble_species(NEUT, PROT);
    }
}

// Terms of the interference diagrams (W on N1 interfering with W on N2): the slot they feed and the bubbles of the
// two nucleon lines (first: momentum transfer q-k, second: momentum transfer k).
struct InterferenceTerm {
    int slot;
    BubbleSpecies first;
    BubbleSpecies second;
};

inline InterferenceTerm interference_term(int t) {
    InterferenceTerm it;
    switch (t) {
        case 0:  it.slot=0; it.first=bubble_species(NEUT,PROT);  it.second=bubble_species(PROT,PROT);  break; // tr3_mu tr1_nu
        case 1:  it.slot=0; it.first=bubble_species(PROT,PROT);   it.second=bubble_species(NEUT,PROT);  break; // tr1_mu tr3_nu
        case 2:  it.slot=1; it.first=bubble_species(NEUT,PROT);  it.second=bubble_species(NEUT,NEUT); break; // tr3_mu tr2_nu
        default: it.slot=2; it.first=bubble_species(NEUT,NEUT); it.second=bubble_species(NEUT,PROT);  break; // tr2_mu tr3_nu
    }
    return it;
}

// Species of the intermediate nucleon when the W is absorbed first (direct nucleon pole) or after the pion
// emission (crossed nucleon pole). It is a Pauli-blocked state if its momentum is below its kF.
inline Nucleon intermediate_direct()  { return CH_2p2h==CCanti ? NEUT : PROT; }
inline Nucleon intermediate_crossed() { return CH_2p2h==CCanti ? PROT : NEUT; }

// Range of the hole energy of a bubble with energy transfer q0 and momentum transfer dq
inline void bubble_limits(double q0, double dq, double kF_hole, double kF_part, double& Emin, double& Emax) {
    Emax = sqrt(M*M+kF_hole*kF_hole);
    double pom = std::max(M, sqrt(M*M+kF_part*kF_part)-q0);
    Emin = std::max(pom, 0.5*(-q0+dq*sqrt(1.-4.*M*M/(q0*q0-dq*dq))));
}

// Union over the three bubbles of the isospin slots of the hole-energy range. 'any' tells whether at least one
// bubble is open (and, if requireLindhard, has a non-zero Lindhard function); lindhard is that of the first open one.
struct BubbleBox {
    double lo, hi;
    double lindhard;
    bool any;
};

inline BubbleBox bubble_box(const FermiGas& fg, double q0, double dq, bool requireLindhard) {
    BubbleBox box;
    box.any = false;
    for (int slot=0; slot<3; slot++) {
        BubbleSpecies b = vertex_bubble(slot);
        double kFh = kF_of(fg, b.hole), kFp = kF_of(fg, b.part);
        double lind = requireLindhard ? DUR(q0, dq, kFh, kFp) : 1.;
        double Emin, Emax;
        bubble_limits(q0, dq, kFh, kFp, Emin, Emax);
        bool open = !(std::abs(lind)<pow(10,-20)) && !(Emin>Emax);
        if (slot==0) { box.lo=Emin; box.hi=Emax; box.lindhard=lind; } // value used if no bubble is open
        if (!open) continue;
        if (!box.any) { box.lo=Emin; box.hi=Emax; box.lindhard=lind; box.any=true; }
        else { if (Emin<box.lo) box.lo=Emin; if (Emax>box.hi) box.hi=Emax; }
    }
    return box;
}


// Nucleon-hole (N) and Delta-hole (D) Lindhard functions of the pion absorbed in the one-body diagrams.
// _ch: charged pion (pi+ for neutrinos, pi- for antineutrinos); _0: neutral pion; N_0abs: the neutral-pion nucleon
// bubble of its absorption on a proton / a neutron (index PROT / NEUT).
struct PionBubbles {
    std::complex<double> N_ch, D_ch;
    std::complex<double> N_0, D_0;
    std::complex<double> N_0abs[2];
};
void pion_bubbles(const FermiGas& fg, double Epi, double kk, PionBubbles& pb, bool withDelta);

// Whether the pion can excite a nucleon-hole pair at all (second bubble)
inline bool pion_bubbles_open(const PionBubbles& pb, double thr) {
    return !(std::abs(pb.N_ch.imag())<thr && std::abs(pb.N_0.imag())<thr);
}

void hadron_2p2h_part1_Epi_thetapi_phipi_kpi_phip_Ep_isospin(double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, const FermiGas& fg, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33);


#endif /* defined(__Pions__hadron_tensor_2p2h__) */
