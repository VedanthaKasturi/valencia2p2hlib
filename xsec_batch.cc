//
//  xsec_batch.cc
//
//  Batch evaluator for the Valencia 2p2h differential cross section.
//
//  Unlike main.cc, which only runs the Monte-Carlo event loop, this program
//  evaluates d(sigma)/dx^10 at phase-space points supplied by the caller.  It
//  is the backend used by NuGenFlow's `cross_section.valencia_2p2h.Valencia2p2h`
//  to train a normalizing flow on this model.
//
//  The ten variables are the ones the generator throws in saveEvents_2p2h
//  (src/2p2h_MC.cpp), in the same order, so that every derived quantity
//  (rho, k_fermi, w, qq, Emin/Emax) can be recomputed here exactly as the
//  generator does.
//
//  Normalisation
//  -------------
//  saveEvents_2p2h stores weight_ch * coeff, with
//    coeff = w * sin(theta_l) * (Emax-Emin) * (qq+3kF) * (Ein-m_l-Q)
//            * Pi * 2Pi * 2Pi * Pi * rmax
//  so that sigma = mean(weight_stored).  Writing V for the product of the ten
//  sampling ranges, coeff / V = sin(theta_l) / (Emax_2 - Emin_2), hence
//    d(sigma)/dx^10 = weight_bare * sin(theta_l) / (Emax_2 - Emin_2)
//  We additionally apply the 1e38/cm2 * 2Pi factor from total_xsec_2p2h so the
//  result is in 1e-38 cm^2 per unit of the ten variables (the 2Pi is the
//  lepton azimuth, which is not one of the ten).
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>

#include "generatormt.h"
#include "1p1h1pion.h"
#include "hadron_tensor_2p2h.h"
#include "hadron_tensor_2p2h_2b.h"
#include "hadron_tensor_2p2h_rho.h"
#include "2p2h_MC.h"
#include "event.h"
#include "utils.h"


static void usage() {
    std::cerr << "Usage: ./xsec_batch <kinematics_file> <result_file> "
                 "<mode_el> <ReDelta> <oxygen> <Global>\n";
    std::cerr << "  kinematics_file: whitespace-separated rows of\n";
    std::cerr << "      Enu r theta_lepton E_lep E_pi p_pi theta_pi phi_pi E_nuc E_nuc2 phi_nuc\n";
    std::cerr << "      (MeV, fm, radians -- the generator's own units)\n";
    std::cerr << "  result_file: one output row per input row,\n";
    std::cerr << "      dsig_pp dsig_pn dsig_np k_fermi qq w valid\n";
    std::cerr << "  mode_el: 0 - CC, 1 - CCanti\n";
    std::cerr << "  ReDelta: 0 - 0 MeV, 1 - -50 rho/rho0 MeV, 2 - 30 rho/rho0 MeV\n";
    std::cerr << "  oxygen : 0 - C12, 1 - O16\n";
    std::cerr << "  Global : 0 - local Fermi gas (default), 1 - global Fermi gas\n";
}


int main(int argc, char *argv[]) {
    if (argc < 7) {
        usage();
        return 1;
    }

    const std::string kin_path = argv[1];
    const std::string res_path = argv[2];
    const int mode_el = std::stoi(argv[3]);
    const int rede    = std::stoi(argv[4]);
    const int oxygen  = std::stoi(argv[5]);
    const bool global_fg = std::stoi(argv[6]) != 0;

    // ---- globals, exactly as main.cc sets them -----------------------------
    M = 940*MeV;
    mpi = 139.57*MeV;

    mpip = mpi;
    mpi0 = mpi;
    mpi2 = mpi*mpi;
    mpi4 = mpi2*mpi2;
    DP_const  = delusive*fstar/mpip*sqrt3;
    CDP_const = delusive*fstar/mpip/sqrt3;

    lepton_mass = mmu;

    CH_2p2h = CC;
    if (mode_el == 1) CH_2p2h = CCanti;

    switch (rede) {
        case 0:  ReDelta = 0;   break;
        case 1:  ReDelta = -50; break;
        case 2:  ReDelta = 30;  break;
        default: ReDelta = 0;   break;
    }

    // ---- nucleus, as saveEvents_2p2h selects it ----------------------------
    FERMI = true;
    Global = global_fg;

    nucleus_str nuc;
    if (oxygen == 1 && !Global)      { NR_PART = 16; nuc = O16;    }
    else if (oxygen == 1 && Global)  { NR_PART = 16; nuc = GFG_16; }
    else if (oxygen == 0 && !Global) { NR_PART = 12; nuc = C12;    }
    else                             { NR_PART = 12; nuc = GFG_12; }

    std::ifstream fin(kin_path);
    if (!fin) {
        std::cerr << "xsec_batch: cannot open " << kin_path << "\n";
        return 2;
    }
    std::ofstream fout(res_path);
    if (!fout) {
        std::cerr << "xsec_batch: cannot open " << res_path << " for writing\n";
        return 2;
    }
    fout << std::scientific << std::setprecision(10);

    std::complex<double> cunuc(0.,0.);
    std::complex<double> udel(0.,0.);

    const double norm = pow(10,38)/cm2*2*Pi;   // matches total_xsec_2p2h

    std::string line;
    long long n_rows = 0;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);

        double Ein, r, theta_lepton, E_lep, Epi, kk, thetapi, phi_pi, Ep, Ep2, phi_p;
        if (!(iss >> Ein >> r >> theta_lepton >> E_lep >> Epi >> kk
                  >> thetapi >> phi_pi >> Ep >> Ep2 >> phi_p)) {
            std::cerr << "xsec_batch: malformed row " << n_rows << "\n";
            fout << "0 0 0 0 0 0 0\n";
            ++n_rows;
            continue;
        }
        ++n_rows;

        bool result_zero = false;

        // --- density and Fermi momentum at r (2p2h_MC.cpp:56-62) ------------
        double rho = getRho(0, &nuc, true)/fermi3 + getRho(0, &nuc, false)/fermi3;
        if (!Global) {
            rho = getRho(r, &nuc, true)/fermi3 + getRho(r, &nuc, false)/fermi3;
        }
        double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);

        // --- lepton kinematics (2p2h_MC.cpp:64-80) --------------------------
        const double corr_in_medium = nuc.Q_2p2h;
        const double w = Ein - E_lep - corr_in_medium;

        double kout2 = E_lep*E_lep - lepton_mass*lepton_mass;
        if (kout2 < 0) { fout << "0 0 0 " << k_fermi << " 0 " << w << " 0\n"; continue; }
        const double kout = sqrt(kout2);
        const double qq = sqrt(Ein*Ein + kout*kout - 2*kout*Ein*cos(theta_lepton));

        if (qq < w) result_zero = true;
        if (qq > 1.2*GeV) result_zero = true;

        // --- pion sector and Lindhard cuts (2p2h_MC.cpp:88-114) -------------
        ulind(Epi/mpi, kk/mpi, k_fermi/mpi, cunuc, udel);
        cunuc *= pow(mpi,2);
        if (std::abs(cunuc.imag()) < pow(10,-20)) result_zero = true;

        const double cosi = cos(thetapi);
        const double sini = sqrt(std::max(1.-cosi*cosi, 0.));
        const double kPi1 = kk*sini*cos(phi_pi);
        const double kPi2 = kk*sini*sin(phi_pi);
        const double kPi3 = kk*cosi;

        const double q_minus_k = sqrt(kPi1*kPi1 + kPi2*kPi2 + (kPi3-qq)*(kPi3-qq));
        const double lindhard = DUR(w-Epi, q_minus_k, k_fermi, k_fermi);
        if (std::abs(lindhard) < pow(10,-20)) result_zero = true;

        // --- nucleon energy limits (2p2h_MC.cpp:120-131) --------------------
        const double Emax = sqrt(M*M + k_fermi*k_fermi);
        const double pom  = std::max(M, Emax-(w-Epi));
        const double Emin = std::max(pom,
            0.5*(-(w-Epi) + q_minus_k*sqrt(1.-4.*M*M/((w-Epi)*(w-Epi)-q_minus_k*q_minus_k))));

        const double Emax_2 = sqrt(M*M + k_fermi*k_fermi);
        const double pom_2  = std::max(M, Emax_2-Epi);
        const double Emin_2 = std::max(pom_2,
            0.5*(-(Epi) + kk*sqrt(1.-4.*M*M/(Epi*Epi-kk*kk))));

        if (!(Emin < Emax))     result_zero = true;
        if (!(Emin_2 < Emax_2)) result_zero = true;

        const double E_range = Emax_2 - Emin_2;

        if (result_zero || !std::isfinite(E_range) || E_range <= 0) {
            fout << "0 0 0 " << k_fermi << " " << qq << " " << w << " 0\n";
            continue;
        }

        // --- the weight itself ---------------------------------------------
        event e;
        event_2p2h_rho_gaussian(r, Ein, theta_lepton, Ep2, Ep, phi_p, kk, phi_pi,
                                thetapi, Epi, w, qq, rho, 1, E_range,
                                &e, &nuc, corr_in_medium);

        // coeff / V = sin(theta_l) / (Emax_2 - Emin_2)
        const double to_density = sin(theta_lepton) / E_range * norm;

        double dsig_pp = e.weight_pp * to_density;
        double dsig_pn = e.weight_pn * to_density;
        double dsig_np = e.weight_np * to_density;

        if (!std::isfinite(dsig_pp)) dsig_pp = 0;
        if (!std::isfinite(dsig_pn)) dsig_pn = 0;
        if (!std::isfinite(dsig_np)) dsig_np = 0;

        fout << dsig_pp << " " << dsig_pn << " " << dsig_np << " "
             << k_fermi << " " << qq << " " << w << " 1\n";
    }

    fout.close();
    return 0;
}
