#include "hadron_tensor_2p2h_2b.h"



/*
 full ISOSPIN
 */
void hadron_2p2h_part2_Epi_thetapi_phipi_kpi_E1_E2_isospin(double Ep2, double Ep1, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, const FermiGas& fg, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33) {
	
	
	//INTEGRATION PARAMETERS
	
	const double rho = fg.rho;
	
	double cosi=cos(thetapi);
	
	D4V<double> q(w,0,0,qq);
	D4V<double> kPi(Epi, kk*sqrt(1-cosi*cosi)*cos(phi_pi), kk*sqrt(1-cosi*cosi)*sin(phi_pi), kk*cosi);
	
	double q_minus_k = sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2));
	
	//the four terms have their own bubbles (hole and particle species, hence their own Fermi momenta)
	InterferenceTerm term[4];
	double lindhard[4], lindhard_kpi[4], Emax_t[4], Emin_t[4], Emax2_t[4], Emin2_t[4];
	bool open[4];
	bool any_open=false;
	for (int t=0; t<4; t++) {
		term[t] = interference_term(t);
		double kFh1=kF_of(fg, term[t].first.hole), kFp1=kF_of(fg, term[t].first.part);
		double kFh2=kF_of(fg, term[t].second.hole), kFp2=kF_of(fg, term[t].second.part);
		lindhard[t] = DUR(q(0)-kPi(0), q_minus_k, kFh1, kFp1);
		lindhard_kpi[t] = DUR(kPi(0), kk, kFh2, kFp2);
		
		//hole-energy ranges of the two bubbles: they are the integration ranges of Ep1 and Ep2
		double Emin1, Emax1, Emin2b, Emax2b;
		bubble_limits(q(0)-kPi(0), q_minus_k, kFh1, kFp1, Emin1, Emax1);
		bubble_limits(kPi(0), kk, kFh2, kFp2, Emin2b, Emax2b);
		
		//normalisation of the energy integrals: the width of the hole-energy range of the bubble in which each energy is integrated
		Emax_t[t]=Emax2b;
		Emin_t[t]=Emin2b;
		Emax2_t[t]=Emax1;
		Emin2_t[t]=Emin1;
		
		open[t] = !(std::abs(lindhard[t])<pow(10,-21)) && !(std::abs(lindhard_kpi[t])<pow(10,-21))
			&& !(Ep1<Emin1) && !(Ep1>Emax1) && !(Ep2<Emin2b) && !(Ep2>Emax2b);
		if (open[t]) any_open=true;
	}
	
	if (!any_open){
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	
	double mom_nuc1 = sqrt(Ep1*Ep1-M*M);
	double cos_theta_nuc1 = ( pow(q(0)-kPi(0),2)-q_minus_k*q_minus_k + 2.*Ep1*(q(0)-kPi(0)) )/(2.*mom_nuc1*q_minus_k);
	if (std::abs(cos_theta_nuc1)>1) std::cerr<<"Error! 1 bubble CosTh="<<cos_theta_nuc1<<"\n";
	double sin_theta_nuc1=sqrt(1.-pow(cos_theta_nuc1,2));
	
	double mom_nuc2 = sqrt(Ep2*Ep2-M*M);
	double cos_theta_nuc2 = ( pow(kPi(0),2)-kk*kk + 2.*Ep2*kPi(0) )/(2.*mom_nuc2*kk);
	if (std::abs(cos_theta_nuc2)>1) std::cerr<<"Error! 2 bubble CosTh="<<cos_theta_nuc2<<"\n";
	double sin_theta_nuc2=sqrt(1.-pow(cos_theta_nuc2,2));
	
	D4V<double> kPiminusq = q-kPi;
	
	DM qslash=q*dirac;
	DM kPislash = kPi*dirac;
	
	channel = 0; //we calculate everything with CG coefficients = 1
	
	//currents that do not depend on the nucleons
	D4V<DM> vector_contact, axial_contact, pion_pole_cur, pion_in_flight_cur;
	D4V<DM> vector_contact_2, axial_contact_2, pion_pole_cur_2, pion_in_flight_cur_2;
	pion_pole(q,kPi,qslash, pion_pole_cur);
	contact_term(q,kPi, vector_contact, axial_contact);
	pion_in_flight(q,kPi, pion_in_flight_cur);
	pion_pole(q,kPiminusq,qslash, pion_pole_cur_2);
	contact_term(q,kPiminusq, vector_contact_2, axial_contact_2);
	pion_in_flight(q,kPiminusq, pion_in_flight_cur_2);
	
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
	
	// Same for neutrinos and antineutrinos: the mirror p<->n is in the species of the terms, not in the tables
	int canal[3]={1,2,3};
	
	//1: nu p --> l- p pi+ \n
	//2: nu n --> l- n pi+ \n
	//3: nu n --> l- p pi0 \n
	
	const bool same_species = (fg.kF_p==fg.kF_n); //then both species share the currents
	//the Pauli blocking of the intermediate nucleon of the direct and of the crossed nucleon pole depends on its species
	bool nucl_direct1=false, nucl_crossed1=false, nucl_direct2=false, nucl_crossed2=false;
	
	/*first bubble: W vertex on N1. Its traces depend on the species of the hole*/
	struct Line {
		D4V<std::complex<double>> tr[3]; //for the three isospin channels
	};
	Line line1[2], line2[2];
	bool built1[2]={false,false}, built2[2]={false,false};
	
	auto build_line1 = [&](Nucleon hole, Line& L) {
		double shift = kF_of(fg,hole)*kF_of(fg,hole)/2./M*FERMI;
		D4V<double> p(Ep1-shift,mom_nuc1*sin_theta_nuc1*cos(phi_p),mom_nuc1*sin_theta_nuc1*sin(phi_p),mom_nuc1*cos_theta_nuc1);
		D4V<double> p_med(Ep1-shift,mom_nuc1*sin_theta_nuc1*cos(phi_p),mom_nuc1*sin_theta_nuc1*sin(phi_p),mom_nuc1*cos_theta_nuc1);
		rotate(p,q-kPi);
		rotate(p_med,q-kPi);
		DM pslash_med=p_med*dirac;
		
		D4V<DM> delta_cur_vec, delta_cur_ax, delta_cur_cross_vec, delta_cur_cross_ax, nucl_cur_vec, nucl_cur_ax,nucl_cur_s,  nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s;
		delta_current_medium(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_vec, delta_cur_ax);
		delta_current_crossed_medium(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_cross_vec, delta_cur_cross_ax);
		
		D4V<double> p_inter=p+q;
		double p_inter_mom = sqrt(p_inter(1)*p_inter(1)+p_inter(2)*p_inter(2)+p_inter(3)*p_inter(3));
		nucl_direct1  = p_inter_mom>kF_of(fg, intermediate_direct());
		nucl_crossed1 = p_inter_mom>kF_of(fg, intermediate_crossed());
		if (nucl_crossed1)
			nucleon_cross_current_medium(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s,shift);
		if (nucl_direct1)
			nucleon_current_medium(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_vec, nucl_cur_ax, nucl_cur_s,shift);
		
		DM trace1 = (pslash_med+M_N*identi)*(kPislash-qslash)*gamma5*(pslash_med+qslash-kPislash+identi*M_N);  //first bubble
		
		for (int c=0; c<3; c++) {
			int ch = canal[c];
			D4V<DM> j1;
			j1 =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax;
			if (nucl_direct1 && nucl_crossed1)
				j1 +=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s+ (CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s;
			else if (nucl_direct1)
				j1 +=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s;
			else if (nucl_crossed1)
				j1 +=(CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s;
			
			j1*=eff_long;
			j1+= (CGCT_v[ch]/CGCT_v[0])*vector_contact + axial_contact*(CGCT_a[ch]/CGCT_a[0])+ pion_pole_cur*(CGPP_a[ch]/CGPP_a[0])+pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]);
			
			L.tr[c] = j1.CT((trace1*j1));
		}
	};
	
	/*second bubble: W vertex on N2. Its traces depend on the species of the hole*/
	auto build_line2 = [&](Nucleon hole, Line& L) {
		double shift = kF_of(fg,hole)*kF_of(fg,hole)/2./M*FERMI;
		D4V<double> l_med(Ep2-shift,0,mom_nuc2*sin_theta_nuc2,mom_nuc2*cos_theta_nuc2);
		rotate(l_med, kPi);
		DM lslash_med=l_med*dirac;
		
		D4V<DM> delta_cur_vec_2, delta_cur_ax_2, delta_cur_cross_vec_2, delta_cur_cross_ax_2, nucl_cur_vec_2, nucl_cur_ax_2, nucl_cur_s_2,  nucl_cur_cross_vec_2, nucl_cur_cross_ax_2, nucl_cur_cross_s_2;
		delta_current_medium(rho, l_med, q, kPiminusq, lslash_med, qslash, delta_cur_vec_2, delta_cur_ax_2);
		delta_current_crossed_medium(rho, l_med, q, kPiminusq, lslash_med, qslash, delta_cur_cross_vec_2, delta_cur_cross_ax_2);
		
		D4V<double> p_inter2=l_med+kPi-q;
		double p_inter_mom2 =  sqrt(p_inter2(1)*p_inter2(1)+p_inter2(2)*p_inter2(2)+p_inter2(3)*p_inter2(3));
		nucl_direct2  = p_inter_mom2>kF_of(fg, intermediate_direct());
		nucl_crossed2 = p_inter_mom2>kF_of(fg, intermediate_crossed());
		if (nucl_crossed2)
			nucleon_cross_current_medium(kPiminusq, l_med, q, lslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_cross_vec_2, nucl_cur_cross_ax_2, nucl_cur_cross_s_2,shift);
		if (nucl_direct2)
			nucleon_current_medium(kPiminusq, l_med, q, lslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_vec_2, nucl_cur_ax_2, nucl_cur_s_2,shift);
		
		DM trace2 = gamma0*(lslash_med+kPislash+identi*M_N)*kPislash*gamma5*(lslash_med+M_N*identi)*gamma0;		//second bubble
		
		for (int c=0; c<3; c++) {
			int ch = canal[c];
			D4V<DM> j1_b;
			j1_b =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec_2 + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax_2 +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec_2 + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax_2;
			if (nucl_direct2 && nucl_crossed2)
				j1_b+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec_2 + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax_2+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s_2+(CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec_2+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax_2+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s_2;
			else if (nucl_direct2)
				j1_b+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec_2 + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax_2+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s_2;
			else if (nucl_crossed2)
				j1_b+=(CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec_2+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax_2+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s_2;
			
			j1_b*=eff_long_kq;
			j1_b+=pion_in_flight_cur_2*(CGPF_v[ch]/CGPF_v[0])+(CGCT_v[ch]/CGCT_v[0])*vector_contact_2 + axial_contact_2*(CGCT_a[ch]/CGCT_a[0]) + pion_pole_cur_2*(CGPP_a[ch]/CGPP_a[0]);
			
			j1_b.hermit();
			L.tr[c] = j1_b.CT(trace2*j1_b);
		}
	};
	
	//the traces of the species that the open terms need
	auto tr_mu = [&](Nucleon hole, int c) -> D4V<std::complex<double>>& {
		int idx = same_species ? 0 : (hole==PROT ? 0 : 1);
		if (!built1[idx]) { build_line1(hole, line1[idx]); built1[idx]=true; }
		return line1[idx].tr[c];
	};
	auto tr_nu = [&](Nucleon hole, int c) -> D4V<std::complex<double>>& {
		int idx = same_species ? 0 : (hole==PROT ? 0 : 1);
		if (!built2[idx]) { build_line2(hole, line2[idx]); built2[idx]=true; }
		return line2[idx].tr[c];
	};
	
	//1 bubble: 1x3
	//2 bubble: 3x1
	//3 bubble: 3x2
	//4 bubble: 2x3
	//term t: channel of the W vertex on N1 (trace tr_mu) and on N2 (trace tr_nu), see interference_term()
	const int ch1_of[4]={2,0,2,1}; //index into canal: channels 3,1,3,2
	const int ch2_of[4]={0,2,1,2}; //channels 1,3,2,3
	
	D4T<std::complex<double>> Wterm[4]={D4T<std::complex<double>>(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),D4T<std::complex<double>>(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),D4T<std::complex<double>>(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),D4T<std::complex<double>>(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)};
	double coeff[4]={0,0,0,0};
	for (int t=0; t<4; t++) {
		if (!open[t]) continue;
		D4V<std::complex<double>>& a = tr_mu(term[t].first.hole, ch1_of[t]);
		D4V<std::complex<double>>& b = tr_nu(term[t].second.hole, ch2_of[t]);
		if (t==2 || t==3) Wterm[t] = D4T<std::complex<double>>(a, b*(-1.));
		else Wterm[t] = D4T<std::complex<double>>(a, b);
		
		coeff[t] =-pow(f_pi*f_pi_kpi_q,2)*lindhard[t]*lindhard_kpi[t]
		
		*pow(gA/2./fpi,2)/onShell/onShell_kpi_q
		
		/pow(M,4)/2./sqrt(2.)*sqrt(1-cosi*cosi)
		
		*kk*kk/pow(2*M_PI,4)/4./(Emax_t[t]-Emin_t[t])/(Emax2_t[t]-Emin2_t[t]); //  /4 from two traces
		//dim: 2-2+2  -2 +2 + 2 = 4
	}
	
	//slot 0 gets the two terms 0 and 1 (W_1 = tr3_mu tr1_nu + tr1_mu tr3_nu)
	auto slot0 = [&](int i, int j) -> std::complex<double> {
		if (open[0] && open[1] && coeff[0]==coeff[1]) {
			D4T<std::complex<double>> W1(tr_mu(term[0].first.hole, ch1_of[0]), tr_nu(term[0].second.hole, ch2_of[0])); //pp
			W1.Add(tr_mu(term[1].first.hole, ch1_of[1]), tr_nu(term[1].second.hole, ch2_of[1])); //pp
			return coeff[0]*W1(i,j);
		}
		return coeff[0]*Wterm[0](i,j)+coeff[1]*Wterm[1](i,j);
	};
	
	A00[0]=slot0(0,0);
	A11[0]=slot0(1,1);
	A03[0]=slot0(0,3);
	A12[0]=slot0(1,2);
	A33[0]=slot0(3,3);
	
	A00[1]=coeff[2]*Wterm[2](0,0);
	A11[1]=coeff[2]*Wterm[2](1,1);
	A03[1]=coeff[2]*Wterm[2](0,3);
	A12[1]=coeff[2]*Wterm[2](1,2);
	A33[1]=coeff[2]*Wterm[2](3,3);
	
	A00[2]=coeff[3]*Wterm[3](0,0);
	A11[2]=coeff[3]*Wterm[3](1,1);
	A03[2]=coeff[3]*Wterm[3](0,3);
	A12[2]=coeff[3]*Wterm[3](1,2);
	A33[2]=coeff[3]*Wterm[3](3,3);
	
	
}

