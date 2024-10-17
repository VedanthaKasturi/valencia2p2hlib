#include "hadron_tensor_2p2h_rho.h"




void hadron_2p2h_rho_Epi_thetapi_phipi_kpi_phip_Ep_isospin(double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, double rho, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33) {
	
	double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);
	
	
	std::complex<double> cunuc(0.,0.);
	std::complex<double> udel(0.,0.);
	std::complex<double> utot(0.,0.);
	
	ulind(Epi/mpi, kk/mpi, k_fermi/mpi,cunuc, udel);
	cunuc*=pow(mpi,2);
	udel*=pow(mpi,2);
	udel=delta_lind (Epi*fermi,kk*fermi, rho*fermi3, k_fermi*fermi);
	udel/=(fermi*fermi);
	
    
	if (std::abs(cunuc.imag())<pow(10,-15)) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	
	
	//INTEGRATION PARAMETERS
	double mom_nuc=sqrt(3./5.)*k_fermi;
	double cosi=cos(thetapi);
	double cos_theta_nuc=0;
	double sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	D4V<double> p(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,0,mom_nuc*sin_theta_nuc,mom_nuc*cos_theta_nuc); //nucleon
	D4V<double> p_med(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,0,mom_nuc*sin_theta_nuc,mom_nuc*cos_theta_nuc); //nucleon
	D4V<double> q(w,0,0,qq);
	D4V<double> kPi(Epi, kk*sqrt(1-cosi*cosi)*cos(phi_pi), kk*sqrt(1-cosi*cosi)*sin(phi_pi), kk*cosi);
	
	double lindhard = DUR(q(0)-kPi(0), sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2)),k_fermi, k_fermi);
	double q_minus_k = sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2));
	
	if (std::abs(lindhard)<pow(10,-15)) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	double aver_en_1st = Ep;
	mom_nuc = sqrt(aver_en_1st*aver_en_1st-M*M);
	cos_theta_nuc = ( pow(q(0)-kPi(0),2)-q_minus_k*q_minus_k + 2.*aver_en_1st*(q(0)-kPi(0)) )/(2.*mom_nuc*q_minus_k);
	
	if (std::abs(cos_theta_nuc)>1) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;
		std::cerr<<"Error RHO! CosTh="<<cos_theta_nuc<<" Ep="<<Ep<<"\n";
	}
	
	sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	p(0)=aver_en_1st-k_fermi*k_fermi/2./M*FERMI;
	p(1)=mom_nuc*sin_theta_nuc*cos(phi_p);
	p(2)=mom_nuc*sin_theta_nuc*sin(phi_p);
	p(3)=mom_nuc*cos_theta_nuc; //nukleon
	rotate(p, q-kPi);
	p_med(0)=aver_en_1st-k_fermi*k_fermi/2./M*FERMI;
	p_med(1)=mom_nuc*sin_theta_nuc*cos(phi_p);
	p_med(2)=mom_nuc*sin_theta_nuc*sin(phi_p);
	p_med(3)=mom_nuc*cos_theta_nuc; //nukleon
	rotate(p_med, q-kPi);
	
	//	}
	
	
	D4V<double> pprime=p+q-kPi;
	
	
	DM pslash=p*dirac;
	DM pslash_med=p_med*dirac;
	
	DM qslash=q*dirac;
	
	channel = 0; //we calculate everything with CG coefficients = 1
	D4T<DM> delta_cur_vec, delta_cur_ax, delta_cur_cross_vec, delta_cur_cross_ax,   nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s, nucl_cur_vec, nucl_cur_ax,nucl_cur_s;
	delta_current_medium_rho(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_vec, delta_cur_ax);
	delta_current_crossed_medium_rho(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_cross_vec, delta_cur_cross_ax);
	
	D4V<double> p_inter=p+q;
	double p_inter_mom = sqrt(p_inter(1)*p_inter(1)+p_inter(2)*p_inter(2)+p_inter(3)*p_inter(3));
	
	if (p_inter_mom>k_fermi) {
		nucleon_cross_current_rho(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s,k_fermi*k_fermi/2./M*FERMI);
		nucleon_current_medium_rho(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_vec, nucl_cur_ax, nucl_cur_s,k_fermi*k_fermi/2./M*FERMI);
	}
	
	double f_rho2 = pow((lambda_rho*lambda_rho-mrho2)/(lambda_rho*lambda_rho-kPi*kPi),2);
	double onShell = kPi*kPi-mrho2;
	
	//CHARGE CURRENT
	//calculate 3 different isospin channels:
	//  n --> p + pi0  channel = 3
	//  n --> n + pi+  channel = 2
	//  p --> p + pi+  channel = 1
	
	//constant that we multiply Delta, DeltaCrossed, Nucleon, NucleonCrossed currents
	
	double eff_trans = 1;
	
    double gprime1 = gprime_t(kPi*kPi,kk*kk);
    eff_trans=1.+ gprime1*onShell/((kk*kk)*f_rho2*crho);
	
	DM pprimeslash=pprime*dirac+M_N*identi;
	
	
	
	D4T<DM> j1,j2, np1, np2, del1, del2;
	//	D4V<DM> j1,j2;
	D4T<std::complex<double>> Wmunu_1(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	D4T<std::complex<double>> Wmunu_2(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	D4T<std::complex<double>> Wmunu_3(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	
	int min_ch=1;
	int max_ch=4;

    if (CH_2p2h==CCanti) {
        min_ch=4;
        max_ch=7;
    }
    

	for (int ch= min_ch; ch<max_ch; ch++) {
		
		j1 =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax;
		
		if (p_inter_mom>k_fermi)
			j1+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+ (CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s;
		
		j1*=eff_trans;
		
		
		j2 = j1;
		
		//	j1.transp();
		j1 &=pprimeslash;
		j1 *=(pslash+M_N*identi);
		
		
		j2.hermit();
		j2 &=gamma0;
		j2 *=gamma0;
		
		
		D4T<DM> Amunu=j1*j2;
		//now we subtract NP-NP and D-D
        
       if (p_inter_mom>k_fermi) {
			np1= (CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s;
			np1*=eff_trans;
			np2=np1;
			np1 &=pprimeslash;
			np1 *=(pslash+M_N*identi);
			np2.hermit();
			np2 &=gamma0;
			np2 *=gamma0;
			
			Amunu+=(np1*np2*(-1.));
		}
        
		
		
		D4T<std::complex<double>> Bmunu(Amunu.CT(Amunu));
        
		if (ch==min_ch) Wmunu_1+=Bmunu;
		else if (ch==min_ch+1) Wmunu_2+=Bmunu;
		else  Wmunu_3+=Bmunu;
	}
	
	
	utot=cunuc+udel;
    gprime1 = gprime_t(kPi*kPi,kk*kk);
	double vt = pow(gA/2./fpi,2)*(crho*f_rho2*kk*kk/onShell +gprime1);
	
	double rhoInc = cunuc.imag();
	
    rhoInc/=pow(std::abs(1.-vt*utot),2);
	
	double Emax=sqrt(M*M+k_fermi*k_fermi);
	
	double pom = std::max(M,Emax-(q(0)-kPi(0)));
	double Emin=std::max(pom, 0.5*(-(q(0)-kPi(0))+q_minus_k*sqrt(1.-4.*M*M/((q(0)-kPi(0))*(q(0)-kPi(0))-q_minus_k*q_minus_k))));
	
	double coeff =-pow(f_rho2/onShell,2)*pow(gA/2./fpi,2)*rhoInc*lindhard*crho*kk*kk
	
	/M/M
	
	*kk*kk/pow(2*M_PI,4)/2.*sqrt(1-cosi*cosi)/(Emax-Emin)/(2.*M_PI); //the last /2 is from trace 2Pi from averaging over cosTheta_nuc
	//dim: 2-2+2  -2 +2 + 2 = 4
	
	
	A00[0]=coeff*(Wmunu_1(0,0)+Wmunu_3(0,0)/2.);
	A11[0]=coeff*(Wmunu_1(1,1)+Wmunu_3(1,1)/2.);
	A03[0]=coeff*(Wmunu_1(0,3)+Wmunu_3(0,3)/2.);
	A12[0]=coeff*(Wmunu_1(1,2)+Wmunu_3(1,2)/2.);
	A33[0]=coeff*(Wmunu_1(3,3)+Wmunu_3(3,3)/2.);
	
	A00[1]=coeff*Wmunu_2(0,0);
	A11[1]=coeff*Wmunu_2(1,1);
	A03[1]=coeff*Wmunu_2(0,3);
	A12[1]=coeff*Wmunu_2(1,2);
	A33[1]=coeff*Wmunu_2(3,3);
	
	A00[2]=coeff*Wmunu_3(0,0)/2.;
	A11[2]=coeff*Wmunu_3(1,1)/2.;
	A03[2]=coeff*Wmunu_3(0,3)/2.;
	A12[2]=coeff*Wmunu_3(1,2)/2.;
	A33[2]=coeff*Wmunu_3(3,3)/2.;
}





/**********************
 
 
 CURRENTS FOR RHO
 
 
 **********************/

void nucleon_current_medium_rho(D4V<double>&kPi, D4V<double>&p, D4V<double>&q, DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4T<DM>& nucleon_vec, D4T<DM>& nucleon_a, D4T<DM>& nucleon_sc, double shift){
	
	D4V<DM> ksigma = ((kPi*dirac)*dirac - dirac*(kPi*dirac))*delusive*0.5;
	double s = (p+q)*(p+q);
    //shift = 0;
	nucleon_vec.Replace(ksigma*(pslash+qslash+M_N), vector_nuc);
	nucleon_vec *= NP_const*CGNP_v[channel]*sqrt(crho)/(s-pow(M-shift,2)+delusive*M*width);
	nucleon_a.Replace(ksigma*(pslash+qslash+M_N), axial_nuc);
	nucleon_a *= NP_const*CGNP_a[channel]*sqrt(crho)/(s-pow(M-shift,2)+delusive*M*width);
	nucleon_sc.Replace(ksigma*(pslash+qslash+M_N), scalar_nuc);
	nucleon_sc *= NP_const*CGNP_s[channel]*sqrt(crho)/(s-pow(M-shift,2)+delusive*M*width);
	
	nucleon_vec.transp(); // \mu \alpha
	nucleon_a.transp();
	nucleon_sc.transp();
}


void nucleon_cross_current_rho(D4V<double>&kPi, D4V<double>&p, D4V<double>&q,DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4T<DM>& nucleon_vec, D4T<DM>& nucleon_a, D4T<DM>& nucleon_sc,double shift){
	
	D4V<double> pprime = p+q-kPi;
	DM pprimeslash = pprime*dirac;
	D4V<DM> ksigma = ((kPi*dirac)*dirac - dirac*(kPi*dirac))*delusive*0.5;
	double s = (p-kPi)*(p-kPi);
	
    //shift = 0;
	nucleon_vec.Replace(vector_nuc*(pprimeslash-qslash+M_N),ksigma*(1./(s-pow(M-shift,2)+delusive*M*width_cross)));
	nucleon_vec *= CNP_const*CGNPC_v[channel]*sqrt(crho);
	nucleon_a.Replace(axial_nuc*(pprimeslash-qslash+M_N),ksigma*(1./(s-pow(M-shift,2)+delusive*M*width_cross)));
	nucleon_a *= CNP_const*CGNPC_a[channel]*sqrt(crho);
	nucleon_sc.Replace(scalar_nuc*(pprimeslash-qslash+M_N),ksigma*(1./(s-pow(M-shift,2)+delusive*M*width_cross)));
	nucleon_sc *= CNP_const*CGNPC_s[channel]*sqrt(crho);
	
}


void delta_current_medium_rho(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4T<DM>& delta_vec, D4T<DM>& delta_ax){
	
	D4V<double> pdelta=p+q;
	DM kPislash = kPi*dirac;
	
	D4T<DM> projD, del_corr, projD1;
	projection_delta(pdelta, projD); //OK
	
	
	
	D4T<DM> verD_a, verD_v;
	vertex_delta(p, q, verD_a, verD_v); //OK
	
	double s = pdelta*pdelta;
	if (s<0) s=0;
	
	std::complex<double> denom = denominator_shifted(pdelta*pdelta,rho);
	delta_32(pdelta, projD);
    projD*=denom;
	
	
	projD1 = projD;
	
	projD1&=(gamma5*kPislash);
	D4T<DM> projD2(gamma5*dirac, kPi*projD);
	
	projD = projD1-projD2;
	
	
	delta_vec = projD*verD_v*DP_const*CGDP_v[channel]*delusive*sqrt(crho);
	delta_ax = projD*verD_a*DP_const*CGDP_a[channel]*delusive*sqrt(crho);
	
	delta_vec.transp(); // \mu \alpha
	delta_ax.transp();
}

void delta_current_crossed_medium_rho(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4T<DM>& delta_vec, D4T<DM>& delta_ax){
	
	D4V<double> pdelta=p-kPi;
	D4V<double> pprime=pdelta+q;
	DM kPislash = kPi*dirac;
	
	D4T<DM> projD, del_corr, projD1;
	projection_delta(pdelta, projD); //OK
	
	D4T<DM> verD_a, verD_v;
	vertex_delta(pprime, q*(-1.), verD_a, verD_v); //OK
	
	verD_v.hermit();
	verD_v&=gamma0;
	verD_v*=gamma0;
	
	verD_a.hermit();
	verD_a&=gamma0;
	verD_a*=gamma0;
	
	
	std::complex<double> denom = denominator_shifted(pdelta*pdelta,rho);
	
	delta_32(pdelta, projD);
    projD*=denom;
	
	projD1 = projD*(gamma5*kPislash);
	D4T<DM> projD2(projD*kPi,gamma5*dirac);
	
	delta_vec = (verD_v*(projD1-projD2))*CDP_const*CGDPC_v[channel]*sqrt(crho)*delusive*(-1.);
	delta_ax =  (verD_a*(projD1-projD2))*CDP_const*CGDPC_a[channel]*sqrt(crho)*delusive*(-1.);
	
	
	
}

