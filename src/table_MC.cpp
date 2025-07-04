//
// Created by Vedantha Srinivas Kasturi on 01.07.2025.
//

#include "table_MC.h"


int isospin_channel = 0;


void hadron_tensor_2p2h_part1_part2_rho_MC ( double w, double qq, nucleus_str *nuc, double* w_tensor, double* w_tensor_2, double* w_tensor_3)
{
    unsigned long int Ntrials = 1000000;
    unsigned long int itr = 0;
    unsigned long int sucess_itr = 0;

    init_deltaFF(-w*w + qq*qq, w);

    D4V<double> q(w, 0, 0, qq);
    D4V<DM> axial_nuc, vector_nuc, scalar_nuc;
    Wnn(q, vector_nuc, axial_nuc, scalar_nuc);

    w_tensor[0] = w_tensor[1] = w_tensor[2] = w_tensor[3] = w_tensor[4] = 0;
    w_tensor_2[0] = w_tensor[1] = w_tensor_2[2] = w_tensor_2[3] = w_tensor_2[4] = 0;
    w_tensor_3[0] = w_tensor_3[1] = w_tensor_3[2] = w_tensor_3[3] = w_tensor_3[4] = 0;

    std::complex<double> A00 = 0, A11 = 0, A12 = 0, A03 = 0, A33 = 0;
    std::complex<double> C00 = 0, C11 = 0, C12 = 0, C03 = 0, C33 = 0;
    std::complex<double> E00 = 0, E11 = 0, E12 = 0, E03 = 0, E33 = 0;

    double rmax = nuc->rmax;
    double factor_rho = 1./2./M_PI;

    bool result_zero = false;
    bool result_zero_lind = false;
    bool result_zero_cunuc = false;

    while (itr < Ntrials) {
        std::complex<double> B00=0, B11=0, B12=0, B03=0, B33=0;
        std::complex<double> D00=0, D11=0, D12=0, D03=0, D33=0;
        std::complex<double> F00=0, F11=0, F12=0, F03=0, F33=0;

        double r = rmax * frandom();
        double rho = getRho(r, nuc, true)/fermi3 + getRho(r, nuc, false)/fermi3;
        double k_fermi = pow(3 * M_PI * M_PI * rho/2, 1.0/3.0);

        double Epi = w *frandom();
        double thetapi = Pi * frandom();
        double phipi  = 2*Pi * frandom();
        double kpi = /*(qq + 3*k_fermi) */ 2*GeV * frandom();
        double phi_p = 2*M_PI* frandom();
        double cosi = cos(thetapi);

        std::complex<double> cunuc(0.,0.);
        std::complex<double> udel(0.,0.);
        std::complex<double> utot(0.,0.);

        ulind(Epi/mpi, kpi/mpi, k_fermi/mpi,cunuc, udel);
        cunuc*=pow(mpi,2);//*fpNN2/mpi/mpi*gprime;

        if (std::abs(cunuc.imag())<pow(10,-20)) result_zero_cunuc = true;


        D4V<double> kPi(Epi, kpi*sqrt(1 - cosi*cosi)*cos(phipi), kpi*sqrt(1-cosi*cosi)*sin(phipi), kpi*cosi);

        double q_minus_k = sqrt(pow(kPi(1), 2) + pow(kPi(2), 2) + pow(kPi(3) - qq, 2));

        double lindhard = DUR(q(0) - kPi(0), sqrt(pow(kPi(1),2) + pow(kPi(2),2) + pow(kPi(3) - qq ,2) ), k_fermi, k_fermi);
        double lindhard_kpi = DUR(kPi(0), kpi, k_fermi, k_fermi);

        if ( ( std::abs(lindhard) < pow(10, -21) ) || ( std::abs(lindhard_kpi) < pow(10, -21) ) ) result_zero_lind = true;

        double Emax = sqrt( M*M + k_fermi*k_fermi);
        double pom = std::max(M, Emax - (q(0) - kPi(0)));
        double Emin = std::max(pom, 0.5 * ( -(q(0) - kPi(0)) + q_minus_k * sqrt(1.-4.*M*M/ ((q(0) - kPi(0))*(q(0) - kPi(0)) - q_minus_k * q_minus_k ))) );

        double Ep1 = Emin + (Emax - Emin) * frandom();

        if ( Emin > Emax ) result_zero = true;

        double Emax2=sqrt( M*M + k_fermi*k_fermi );
        double pom2 = std::max(M, Emax2 - kPi(0) );
        double Emin2=std::max(pom2, 0.5*(-kPi(0)+kpi*sqrt(1.-4.*M*M/(kPi(0)*kPi(0)-kpi*kpi))));

        double Ep2 = Emin2 + (Emax2 - Emin2) * frandom();

        if (Emin2 > Emax2 ) result_zero = true;

        if (!result_zero && !result_zero_lind && !result_zero_cunuc) {
            //sucess_itr;

            std::complex<double> B00_t[3], D00_t[3], F00_t[3];
            std::complex<double> B11_t[3], D11_t[3], F11_t[3];
            std::complex<double> B12_t[3], D12_t[3], F12_t[3];
            std::complex<double> B03_t[3], D03_t[3], F03_t[3];
            std::complex<double> B33_t[3], D33_t[3], F33_t[3];

            /* Pion Exchange */
            hadron_2p2h_part1_Epi_thetapi_phipi_kpi_phip_Ep_isospin ( Ep1, phi_p, kpi, phipi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, B00_t, B11_t, B12_t, B03_t, B33_t);

            /* Rho Exchange */
            hadron_2p2h_rho_Epi_thetapi_phipi_kpi_phip_Ep_isospin ( Ep1, phi_p, kpi, phipi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, F00_t, F11_t, F12_t, F03_t, F33_t);

            /* 2 different Bubbles */
            hadron_2p2h_part2_Epi_thetapi_phipi_kpi_E1_E2_isospin (Ep2, Ep1, phi_p, kpi, phipi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, D00_t, D11_t, D12_t, D03_t, D33_t);



            if (isospin_channel == 0) {

                B00 = B00_t[0];
                B03 = B03_t[0];
                B11 = B11_t[0];
                B12 = B12_t[0];
                B33 = B33_t[0];

                D00 = D00_t[0];
                D03 = D03_t[0];
                D11 = D11_t[0];
                D12 = D12_t[0];
                D33 = D33_t[0];

                F00 = F00_t[0];
                F03 = F03_t[0];
                F11 = F11_t[0];
                F12 = F12_t[0];
                F33 = F33_t[0];

            } else if (isospin_channel == 1) {
                B00 = B00_t[1];
                B03 = B03_t[1];
                B11 = B11_t[1];
                B12 = B12_t[1];
                B33 = B33_t[1];

                D00 = D00_t[1];
                D03 = D03_t[1];
                D11 = D11_t[1];
                D12 = D12_t[1];
                D33 = D33_t[1];


                F00 = F00_t[1];
                F03 = F03_t[1];
                F11 = F11_t[1];
                F12 = F12_t[1];
                F33 = F33_t[1];
            } else if (isospin_channel == 2) {
                B00 = B00_t[2];
                B03 = B03_t[2];
                B11 = B11_t[2];
                B12 = B12_t[2];
                B33 = B33_t[2];

                D00 = D00_t[2];
                D03 = D03_t[2];
                D11 = D11_t[2];
                D12 = D12_t[2];
                D33 = D33_t[2];


                F00 = F00_t[2];
                F03 = F03_t[2];
                F11 = F11_t[2];
                F12 = F12_t[2];
                F33 = F33_t[2];

            } else {
                for (int channel_number = 0; channel_number < 3; channel_number++ ) {
                    B00 += B00_t[channel_number];
                    B03 += B03_t[channel_number];
                    B11 += B11_t[channel_number];
                    B12 += B12_t[channel_number];
                    B33 += B33_t[channel_number];

                    D00 += D00_t[channel_number];
                    D03 += D03_t[channel_number];
                    D11 += D11_t[channel_number];
                    D12 += D12_t[channel_number];
                    D33 += D33_t[channel_number];

                    F00 += F00_t[channel_number];
                    F03 += F03_t[channel_number];
                    F11 += F11_t[channel_number];
                    F12 += F12_t[channel_number];
                    F33 += F33_t[channel_number];
                }
            }

            double coeff_pion_exchange = rmax * w * Pi * 2*Pi /* (qq + 3*k_fermi) */ * 2*GeV * 2*M_PI * (Emax - Emin) ;
            double coeff_rho_exchange = rmax * w * Pi * 2*Pi /* (qq + 3*k_fermi) */ * 2*GeV * 2*M_PI * (Emax - Emin) ;
            double coeff_2b = rmax * w * Pi * 2*Pi /* (qq + 3*k_fermi) */ * 2*GeV * (Emax - Emin) * (Emax2 - Emin2) ;

            B00 *= coeff_pion_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            B03 *= coeff_pion_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            B11 *= coeff_pion_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            B12 *= coeff_pion_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            B33 *= coeff_pion_exchange * r * r * pow(fermi, 3) * 4 * M_PI;

            D00 *= coeff_2b * r * r * pow(fermi, 3) * 4 * M_PI;
            D03 *= coeff_2b * r * r * pow(fermi, 3) * 4 * M_PI;
            D11 *= coeff_2b * r * r * pow(fermi, 3) * 4 * M_PI;
            D12 *= coeff_2b * r * r * pow(fermi, 3) * 4 * M_PI;
            D33 *= coeff_2b * r * r * pow(fermi, 3) * 4 * M_PI;

            F00 *= coeff_rho_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            F03 *= coeff_rho_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            F11 *= coeff_rho_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            F12 *= coeff_rho_exchange * r * r * pow(fermi, 3) * 4 * M_PI;
            F33 *= coeff_rho_exchange * r * r * pow(fermi, 3) * 4 * M_PI;

        }

        A00 += B00;
        A03 += B03;
        A11 += B11;
        A12 += B12;
        A33 += B33;

        C00 += D00;
        C03 += D03;
        C11 += D11;
        C12 += D12;
        C33 += D33;

        E00 += F00;
        E03 += F03;
        E11 += F11;
        E12 += F12;
        E33 += F33;

        result_zero = false;
        result_zero_lind = false;
        result_zero_cunuc = false;
        itr++;
    }

    w_tensor[0] = factor_rho * A00.real() / itr;
    w_tensor[1] = factor_rho * A03.real() / itr;
    w_tensor[2] = factor_rho * A11.real() / itr;
    w_tensor[3] = factor_rho * A12.imag() / itr;
    w_tensor[4] = factor_rho * A33.real() / itr;

    w_tensor_2[0] = factor_rho * C00.real() / itr;
    w_tensor_2[1] = factor_rho * C03.real() / itr;
    w_tensor_2[2] = factor_rho * C11.real() / itr;
    w_tensor_2[3] = factor_rho * C12.imag() / itr;
    w_tensor_2[4] = factor_rho * C33.real() / itr;

    w_tensor_3[0] = factor_rho * E00.real() / itr;
    w_tensor_3[1] = factor_rho * E03.real() / itr;
    w_tensor_3[2] = factor_rho * E11.real() / itr;
    w_tensor_3[3] = factor_rho * E12.imag() / itr;
    w_tensor_3[4] = factor_rho * E33.real() / itr;



}

void hadron_tensor_3p3h_MC ( double w, double qq, nucleus_str *nuc, double* w_tensor )
{
    unsigned long int Ntrials = 100000;
    unsigned long int itr = 0;
    init_deltaFF(-w*w + qq*qq, w);

    w_tensor[0] = w_tensor[1] = w_tensor[2] = w_tensor[3] = w_tensor[4] = 0;

    std::complex<double> A00 = 0, A11 = 0, A12 = 0, A03 = 0, A33 = 0;

    double rmax = nuc->rmax;
    double factor_rho = 1./2./M_PI;

    while(itr<Ntrials) {
        std::complex<double> B00=0, B11=0, B12=0, B03=0, B33=0;
        double r = rmax * frandom();
        double rho = getRho(r, nuc, true)/fermi3 + getRho(r, nuc, false)/fermi3;
        double k_fermi = pow(3 * M_PI * M_PI * rho/2, 1.0/3.0);

        double p_nuc = frandom() * k_fermi;
        double theta_nuc = frandom() * M_PI;
        double phi_nuc = frandom() * 2 * M_PI;


        hadron_delta_p_phi_theta ( theta_nuc, phi_nuc, p_nuc, w,  qq, rho, true, false, false, B00, B03, B11, B12, B33);

        double coeff = rmax * k_fermi * M_PI * 2 * M_PI;

        B00 *= coeff * r * r * pow(fermi, 3) * 4 * M_PI;
        B03 *= coeff * r * r * pow(fermi, 3) * 4 * M_PI;
        B11 *= coeff * r * r * pow(fermi, 3) * 4 * M_PI;
        B12 *= coeff * r * r * pow(fermi, 3) * 4 * M_PI;
        B33 *= coeff * r * r * pow(fermi, 3) * 4 * M_PI;

        A00 += B00;
        A03 += B03;
        A11 += B11;
        A12 += B12;
        A33 += B33;

        itr++;
    }

    w_tensor[0] = factor_rho * A00.real() / Ntrials;
    w_tensor[1] = factor_rho * A03.real() / Ntrials;
    w_tensor[2] = factor_rho * A11.real() / Ntrials;
    w_tensor[3] = factor_rho * A12.imag() / Ntrials;
    w_tensor[4] = factor_rho * A33.real() / Ntrials;

}
