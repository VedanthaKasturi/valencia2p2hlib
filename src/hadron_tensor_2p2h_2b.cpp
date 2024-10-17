#include "hadron_tensor_2p2h_2b.h"



/*
 full ISOSPIN
 */
void hadron_2p2h_part2_Epi_thetapi_phipi_kpi_E1_E2_isospin(double Ep2, double Ep1, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, double rho, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33) {
	
	
	//INTEGRATION PARAMETERS
	
	double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);

	double mom_nuc=sqrt(3./5.)*k_fermi;
	
	double cosi=cos(thetapi);
	
	double cos_theta_nuc = 0;
	
	double sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	D4V<double> p(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,0,mom_nuc*sin_theta_nuc,mom_nuc*cos_theta_nuc); //nucleon
	D4V<double> p_med(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,0,mom_nuc*sin_theta_nuc,mom_nuc*cos_theta_nuc); //nucleon
	D4V<double> l_med(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,0,mom_nuc*sin_theta_nuc,mom_nuc*cos_theta_nuc); //nucleon
	D4V<double> q(w,0,0,qq);
	D4V<double> kPi(Epi, kk*sqrt(1-cosi*cosi)*cos(phi_pi), kk*sqrt(1-cosi*cosi)*sin(phi_pi), kk*cosi);
	
	double q_minus_k = sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2));
	double lindhard = DUR(q(0)-kPi(0), q_minus_k, k_fermi, k_fermi);
	double lindhard_kpi = DUR(kPi(0), kk, k_fermi, k_fermi);
	
	if (std::abs(lindhard)<pow(10,-21) || std::abs(lindhard_kpi)<pow(10,-21)){
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	
	mom_nuc = sqrt(Ep1*Ep1-M*M);
	cos_theta_nuc = ( pow(q(0)-kPi(0),2)-q_minus_k*q_minus_k + 2.*Ep1*(q(0)-kPi(0)) )/(2.*mom_nuc*q_minus_k);
	if (std::abs(cos_theta_nuc)>1) std::cerr<<"Error! 1 bubble CosTh="<<cos_theta_nuc<<"\n";
	sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	p(0)=Ep1-k_fermi*k_fermi/2./M*FERMI;
	p(1)=mom_nuc*sin_theta_nuc*cos(phi_p);
	p(2)=mom_nuc*sin_theta_nuc*sin(phi_p);
	p(3)=mom_nuc*cos_theta_nuc;
	p_med(0)=Ep1-k_fermi*k_fermi/2./M*FERMI;
	p_med(1)=mom_nuc*sin_theta_nuc*cos(phi_p);
	p_med(2)=mom_nuc*sin_theta_nuc*sin(phi_p);
	p_med(3)=mom_nuc*cos_theta_nuc;
	rotate(p,q-kPi);
	rotate(p_med,q-kPi);
    
    
	mom_nuc = sqrt(Ep2*Ep2-M*M);
    
	
	cos_theta_nuc = ( pow(kPi(0),2)-kk*kk + 2.*Ep2*kPi(0) )/(2.*mom_nuc*kk);
	if (std::abs(cos_theta_nuc)>1) std::cerr<<"Error! 2 bubble CosTh="<<cos_theta_nuc<<"\n";
	sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	

	l_med(0)=Ep2-k_fermi*k_fermi/2./M*FERMI;
    l_med(1)=0;
    l_med(2)=mom_nuc*sin_theta_nuc;
	l_med(3)=mom_nuc*cos_theta_nuc;
	rotate(l_med, kPi);

	D4V<double> kPiminusq = q-kPi;
	
	
	DM pslash_med=p_med*dirac;
	DM lslash_med=l_med*dirac;
	
	DM qslash=q*dirac;
	DM kPislash = kPi*dirac;
	
	channel = 0; //we calculate everything with CG coefficients = 1
	/*first bubble*/
	D4V<DM> delta_cur_vec, delta_cur_ax, delta_cur_cross_vec, delta_cur_cross_ax, nucl_cur_vec, nucl_cur_ax,nucl_cur_s,  nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s, vector_contact, axial_contact, pion_pole_cur, pion_in_flight_cur;
	
	delta_current_medium(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_vec, delta_cur_ax);
	delta_current_crossed_medium(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_cross_vec, delta_cur_cross_ax);
	
	pion_pole(q,kPi,qslash, pion_pole_cur);
	
	contact_term(q,kPi, vector_contact, axial_contact);
	pion_in_flight(q,kPi, pion_in_flight_cur);
	
	D4V<double> p_inter=p+q;
	double p_inter_mom = sqrt(p_inter(1)*p_inter(1)+p_inter(2)*p_inter(2)+p_inter(3)*p_inter(3));
    D4V<double> p_inter2=l_med+kPi-q;
    double p_inter_mom2 =  sqrt(p_inter2(1)*p_inter2(1)+p_inter2(2)*p_inter2(2)+p_inter2(3)*p_inter2(3));
	
	
	if (p_inter_mom>k_fermi) {
		nucleon_cross_current_medium(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s,k_fermi*k_fermi/2./M*FERMI);
		nucleon_current_medium(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_vec, nucl_cur_ax, nucl_cur_s,k_fermi*k_fermi/2./M*FERMI);
	}
	
	/*second bubble*/
	D4V<DM> delta_cur_vec_2, delta_cur_ax_2, delta_cur_cross_vec_2, delta_cur_cross_ax_2, nucl_cur_vec_2, nucl_cur_ax_2, nucl_cur_s_2,  nucl_cur_cross_vec_2, nucl_cur_cross_ax_2, nucl_cur_cross_s_2, vector_contact_2, axial_contact_2, pion_pole_cur_2, pion_in_flight_cur_2;
	
	delta_current_medium(rho, l_med, q, kPiminusq, lslash_med, qslash, delta_cur_vec_2, delta_cur_ax_2);
    delta_current_crossed_medium(rho, l_med, q, kPiminusq, lslash_med, qslash, delta_cur_cross_vec_2, delta_cur_cross_ax_2);
	
	pion_pole(q,kPiminusq,qslash, pion_pole_cur_2);
	contact_term(q,kPiminusq, vector_contact_2, axial_contact_2);
	pion_in_flight(q,kPiminusq, pion_in_flight_cur_2);
 	
	if (p_inter_mom2>k_fermi) {
		nucleon_cross_current_medium(kPiminusq, l_med, q, lslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_cross_vec_2, nucl_cur_cross_ax_2, nucl_cur_cross_s_2,k_fermi*k_fermi/2./M*FERMI);
		nucleon_current_medium(kPiminusq, l_med, q, lslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_vec_2, nucl_cur_ax_2, nucl_cur_s_2,k_fermi*k_fermi/2./M*FERMI);
	}
	
	
	double f_pi = (lambda_pi*lambda_pi-mpi*mpi)/(lambda_pi*lambda_pi-kPi*kPi);
	double f_pi_kpi_q =(lambda_pi*lambda_pi-mpi*mpi)/(lambda_pi*lambda_pi-kPiminusq*kPiminusq);
	double onShell = kPi*kPi-mpi*mpi;
	double onShell_kpi_q = kPiminusq*kPiminusq-mpi*mpi;
	
	
	//CHARGE CURRENT
	//calculate 3 different isospin channels:
	//  n --> p + pi0  channel = 3
	//  n --> n + pi+  channel = 2
	//  p --> p + pi+  channel = 1
	
	//constant that we multiply Delta, DeltaCrossed, Nucleon, NucleonCrossed currents
	
	double eff_long=1;
	double eff_long_kq = 1;
	double gprime1 = gprime_l(kPi*kPi, kk*kk);
    double gprime2 = gprime_l(kPiminusq*kPiminusq, q_minus_k*q_minus_k);
    eff_long = 1.+ gprime1*onShell/((kk*kk)*f_pi*f_pi);
    eff_long_kq = 1.+ gprime2*onShell_kpi_q/((q_minus_k*q_minus_k)*f_pi_kpi_q*f_pi_kpi_q);
	
	D4V<DM> j1,j2,j3,j4, j1_b, j2_b, j3_b, j4_b;
	
	int canal[4]={1,2,3,4};
    
	
	int ch= canal[0];
	
	j1 =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax;
	if (p_inter_mom>k_fermi)
        j1 +=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s+ (CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s;
	
	j1*=eff_long;
	j1+= (CGCT_v[ch]/CGCT_v[0])*vector_contact + axial_contact*(CGCT_a[ch]/CGCT_a[0])+ pion_pole_cur*(CGPP_a[ch]/CGPP_a[0])+pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]);
    
	
	j1_b =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec_2 + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax_2 +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec_2 + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax_2;
	if (p_inter_mom2>k_fermi)
        j1_b+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec_2 + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax_2+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s_2+(CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec_2+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax_2+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s_2;
	
	
	j1_b*=eff_long_kq;
	j1_b+=pion_in_flight_cur_2*(CGPF_v[ch]/CGPF_v[0])+(CGCT_v[ch]/CGCT_v[0])*vector_contact_2 + axial_contact_2*(CGCT_a[ch]/CGCT_a[0]) + pion_pole_cur_2*(CGPP_a[ch]/CGPP_a[0]);
	
	ch= canal[1];
	
	j2 =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax;
	if (p_inter_mom>k_fermi)
        j2+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s+ (CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s;
	
	j2*=eff_long;
	j2+= (CGCT_v[ch]/CGCT_v[0])*vector_contact + axial_contact*(CGCT_a[ch]/CGCT_a[0]) + pion_pole_cur*(CGPP_a[ch]/CGPP_a[0])+pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]);
     
	j2_b =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec_2 + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax_2 +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec_2 + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax_2;
	if (p_inter_mom2>k_fermi)
        j2_b+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec_2 + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax_2+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s_2+(CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec_2+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax_2+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s_2;
	
	j2_b*=eff_long_kq;
	j2_b+=pion_in_flight_cur_2*(CGPF_v[ch]/CGPF_v[0]) +(CGCT_v[ch]/CGCT_v[0])*vector_contact_2 + axial_contact_2*(CGCT_a[ch]/CGCT_a[0]) +  pion_pole_cur_2*(CGPP_a[ch]/CGPP_a[0]);
	ch= canal[2]; //3
	
	j3 =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax;
	if (p_inter_mom>k_fermi)
        j3+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s+ (CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s;
	
	j3*=eff_long;
	j3+=(CGCT_v[ch]/CGCT_v[0])*vector_contact + axial_contact*(CGCT_a[ch]/CGCT_a[0]) +  pion_pole_cur*(CGPP_a[ch]/CGPP_a[0])+pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]);
	
	j3_b =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec_2 + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax_2 +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec_2 + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax_2;
	if (p_inter_mom2>k_fermi)
        j3_b+= (CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec_2 + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax_2+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s_2+(CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec_2+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax_2+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s_2;
	
	j3_b*=eff_long_kq;
	j3_b+= pion_in_flight_cur_2*(CGPF_v[ch]/CGPF_v[0]) + (CGCT_v[ch]/CGCT_v[0])*vector_contact_2 + axial_contact_2*(CGCT_a[ch]/CGCT_a[0]) +pion_pole_cur_2*(CGPP_a[ch]/CGPP_a[0]);
 
	
	DM trace1 = (pslash_med+M_N*identi)*(kPislash-qslash)*gamma5*(pslash_med+qslash-kPislash+identi*M_N);  //first bubble
	DM trace2 = gamma0*(lslash_med+kPislash+identi*M_N)*kPislash*gamma5*(lslash_med+M_N*identi)*gamma0;		//second bubble
	
	//1: nu p --> l- p pi+ \n
	//2: nu n --> l- n pi+ \n
	//3: nu n --> l- p pi0 \n
	//1 bubble: 1x3
	//2 bubble: 3x1
	//3 bubble: 3x2
	//4 bubble: 2x3
	
	
	
	D4V<std::complex<double>> tr1_nu, tr2_nu, tr3_nu, tr4_nu, tr1_mu, tr2_mu, tr3_mu, tr4_mu;
	
	
	tr1_mu = j1.CT((trace1*j1));
	tr2_mu = j2.CT((trace1*j2));
	tr3_mu = j3.CT((trace1*j3));
	tr4_mu = j4.CT((trace1*j4));
	
	
	j1_b.hermit();
	j2_b.hermit();
	j3_b.hermit();
	j4_b.hermit();
	
	tr1_nu = j1_b.CT(trace2*j1_b);
	tr2_nu = j2_b.CT(trace2*j2_b);
	tr3_nu = j3_b.CT(trace2*j3_b);
	tr4_nu = j4_b.CT(trace2*j4_b);
	
	D4T<std::complex<double>> Wmunu_1(tr3_mu,tr1_nu); //pp
	Wmunu_1.Add(tr1_mu,tr3_nu); //pp
	D4T<std::complex<double>> Wmunu_3(tr2_mu,tr3_nu*(-1.)); //pn
	D4T<std::complex<double>> Wmunu_2(tr3_mu,tr2_nu*(-1.)); //np

	
	
	double Emax=sqrt(M*M+k_fermi*k_fermi);
	
	double pom = std::max(M,Emax-kPi(0));
	double Emin=std::max(pom, 0.5*(-kPi(0)+kk*sqrt(1.-4.*M*M/(kPi(0)*kPi(0)-kk*kk))));
	
	double Emax2=sqrt(M*M+k_fermi*k_fermi);
	

	double Emin2=std::max(pom, 0.5*(-(q(0)-kPi(0))+q_minus_k*sqrt(1.-4.*M*M/((q(0)-kPi(0))*(q(0)-kPi(0))-q_minus_k*q_minus_k))));
	
	double coeff =-pow(f_pi*f_pi_kpi_q,2)*lindhard*lindhard_kpi
	
	*pow(gA/2./fpi,2)/onShell/onShell_kpi_q
	
	/pow(M,4)/2./sqrt(2.)*sqrt(1-cosi*cosi)
	
	*kk*kk/pow(2*M_PI,4)/4./(Emax-Emin)/(Emax2-Emin2); //  /4 from two traces
	//dim: 2-2+2  -2 +2 + 2 = 4
	
	
	A00[0]=coeff*Wmunu_1(0,0);
	A11[0]=coeff*Wmunu_1(1,1);
	A03[0]=coeff*Wmunu_1(0,3);
	A12[0]=coeff*Wmunu_1(1,2);
	A33[0]=coeff*Wmunu_1(3,3);
	
	A00[1]=coeff*Wmunu_2(0,0);
	A11[1]=coeff*Wmunu_2(1,1);
	A03[1]=coeff*Wmunu_2(0,3);
	A12[1]=coeff*Wmunu_2(1,2);
	A33[1]=coeff*Wmunu_2(3,3);
	
	A00[2]=coeff*Wmunu_3(0,0);
	A11[2]=coeff*Wmunu_3(1,1);
	A03[2]=coeff*Wmunu_3(0,3);
	A12[2]=coeff*Wmunu_3(1,2);
	A33[2]=coeff*Wmunu_3(3,3);
	
	
}

