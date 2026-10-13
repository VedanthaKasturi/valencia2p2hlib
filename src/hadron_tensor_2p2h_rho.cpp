#include "hadron_tensor_2p2h_rho.h"




void hadron_2p2h_rho_Epi_thetapi_phipi_kpi_phip_Ep_isospin(double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, const FermiGas& fg, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33) {
	
	const double rho = fg.rho;
	
	//Lindhard functions of the pion that is absorbed (second bubble), for its charge
	PionBubbles pb;
	pion_bubbles(fg, Epi, kk, pb, true);
	
    
	if (!pion_bubbles_open(pb, pow(10,-15))) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	
	
	//INTEGRATION PARAMETERS
	double cosi=cos(thetapi);
	double cos_theta_nuc=0;
	double sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	D4V<double> q(w,0,0,qq);
	D4V<double> kPi(Epi, kk*sqrt(1-cosi*cosi)*cos(phi_pi), kk*sqrt(1-cosi*cosi)*sin(phi_pi), kk*cosi);
	
	double q_minus_k = sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2));
	
	//the three isospin slots have their own bubble (hole and particle species, hence their own Fermi momenta)
	BubbleSpecies vb[3];
	double kF_hole[3], lindhard[3], Emin[3], Emax[3];
	bool open[3];
	bool any_open=false;
	for (int k=0; k<3; k++) {
		vb[k] = vertex_bubble(k);
		kF_hole[k] = kF_of(fg, vb[k].hole);
		lindhard[k] = DUR(q(0)-kPi(0), sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2)), kF_hole[k], kF_of(fg, vb[k].part));
		bubble_limits(q(0)-kPi(0), q_minus_k, kF_hole[k], kF_of(fg, vb[k].part), Emin[k], Emax[k]);
		open[k] = !(std::abs(lindhard[k])<pow(10,-15)) && !(Ep<Emin[k]) && !(Ep>Emax[k]);
		if (open[k]) any_open=true;
	}
	
	if (!any_open) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	double aver_en_1st = Ep;
	double mom_nuc = sqrt(aver_en_1st*aver_en_1st-M*M);
	cos_theta_nuc = ( pow(q(0)-kPi(0),2)-q_minus_k*q_minus_k + 2.*aver_en_1st*(q(0)-kPi(0)) )/(2.*mom_nuc*q_minus_k);
	
	if (std::abs(cos_theta_nuc)>1) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;
	}
	
	sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	channel = 0; //we calculate everything with CG coefficients = 1
	DM qslash=q*dirac;
	
	//the Pauli blocking of the intermediate nucleon of the direct and of the crossed nucleon pole depends on its species
	bool nucl_direct=false, nucl_crossed=false;
	
	//currents of the W vertex on a hole; they depend on the species of the hole through its removal-energy shift only
	struct HoleCurrents {
		DM pslash, pslash_med, pprimeslash;
		D4T<DM> delta_cur_vec, delta_cur_ax, delta_cur_cross_vec, delta_cur_cross_ax, nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s, nucl_cur_vec, nucl_cur_ax, nucl_cur_s;
	};
	HoleCurrents hole_cur[2];
	bool hole_built[2]={false,false};
	const bool same_species = (fg.kF_p==fg.kF_n); //then both species share the currents
	
	for (int k=0; k<3; k++) {
		if (!open[k]) continue;
		int idx = same_species ? 0 : (vb[k].hole==PROT ? 0 : 1);
		if (hole_built[idx]) continue;
		hole_built[idx]=true;
		HoleCurrents& H = hole_cur[idx];
		double shift = kF_hole[k]*kF_hole[k]/2./M*FERMI;
		
		D4V<double> p(aver_en_1st-shift,mom_nuc*sin_theta_nuc*cos(phi_p),mom_nuc*sin_theta_nuc*sin(phi_p),mom_nuc*cos_theta_nuc); //nukleon
		rotate(p, q-kPi);
		D4V<double> p_med(aver_en_1st-shift,mom_nuc*sin_theta_nuc*cos(phi_p),mom_nuc*sin_theta_nuc*sin(phi_p),mom_nuc*cos_theta_nuc); //nukleon
		rotate(p_med, q-kPi);
		
		D4V<double> pprime=p+q-kPi;
		H.pslash=p*dirac;
		H.pslash_med=p_med*dirac;
		H.pprimeslash=pprime*dirac+M_N*identi;
		
		delta_current_medium_rho(rho, p_med, q, kPi, H.pslash_med, qslash, H.delta_cur_vec, H.delta_cur_ax);
		delta_current_crossed_medium_rho(rho, p_med, q, kPi, H.pslash_med, qslash, H.delta_cur_cross_vec, H.delta_cur_cross_ax);
		
		D4V<double> p_inter=p+q;
		double p_inter_mom = sqrt(p_inter(1)*p_inter(1)+p_inter(2)*p_inter(2)+p_inter(3)*p_inter(3));
		nucl_direct  = p_inter_mom>kF_of(fg, intermediate_direct());
		nucl_crossed = p_inter_mom>kF_of(fg, intermediate_crossed());
		
		if (nucl_crossed)
			nucleon_cross_current_rho(kPi, p_med, q, H.pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, H.nucl_cur_cross_vec, H.nucl_cur_cross_ax, H.nucl_cur_cross_s,shift);
		if (nucl_direct)
			nucleon_current_medium_rho(kPi, p_med, q, H.pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, H.nucl_cur_vec, H.nucl_cur_ax, H.nucl_cur_s,shift);
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
	
	
	
	D4T<DM> j1,j2, np1, np2, del1, del2;
	//	D4V<DM> j1,j2;
	D4T<std::complex<double>> Wmunu_1(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	D4T<std::complex<double>> Wmunu_2(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	D4T<std::complex<double>> Wmunu_3(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	D4T<std::complex<double>>* Wmunu[3]={&Wmunu_1,&Wmunu_2,&Wmunu_3};
	
	// Isospin channels entering the three slots W_1, W_2, W_3 (A[0], A[1], A[2]). The antineutrino channels are
	// the p<->n mirror of the neutrino ones (4==2, 5==1, 6==-3 in constants.h), so slot 1 takes channel 5
	// (n -> n pi-, final nn), slot 2 channel 4 (p -> p pi-, final pn) and slot 3 channel 6 (p -> n pi0).
	int chans[3]={1,2,3};
	if (CH_2p2h==CCanti) {
		chans[0]=5;
		chans[1]=4;
		chans[2]=6;
	}
    

	for (int k=0; k<3; k++) {
		if (!open[k]) continue;
		int ch=chans[k];
		HoleCurrents& H = hole_cur[same_species ? 0 : (vb[k].hole==PROT ? 0 : 1)];
		
		j1 =(CGDP_v[ch]/CGDP_v[0])*H.delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*H.delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*H.delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*H.delta_cur_cross_ax;
		
		if (nucl_direct && nucl_crossed)
			j1+=(CGNP_v[ch]/CGNP_v[0])*H.nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*H.nucl_cur_ax+ (CGNPC_v[ch]/CGNPC_v[0])*H.nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*H.nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*H.nucl_cur_cross_s+(CGNP_s[ch]/CGNP_s[0])*H.nucl_cur_s;
		else if (nucl_direct)
			j1+=(CGNP_v[ch]/CGNP_v[0])*H.nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*H.nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*H.nucl_cur_s;
		else if (nucl_crossed)
			j1+=(CGNPC_v[ch]/CGNPC_v[0])*H.nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*H.nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*H.nucl_cur_cross_s;
		
		j1*=eff_trans;
		
		
		j2 = j1;
		
		//	j1.transp();
		j1 &=H.pprimeslash;
		j1 *=(H.pslash+M_N*identi);
		
		
		j2.hermit();
		j2 &=gamma0;
		j2 *=gamma0;
		
		
		D4T<DM> Amunu=j1*j2;
		//now we subtract NP-NP and D-D
        
       if (nucl_direct) {
			np1= (CGNP_v[ch]/CGNP_v[0])*H.nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*H.nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*H.nucl_cur_s;
			np1*=eff_trans;
			np2=np1;
			np1 &=H.pprimeslash;
			np1 *=(H.pslash+M_N*identi);
			np2.hermit();
			np2 &=gamma0;
			np2 *=gamma0;
			
			Amunu+=(np1*np2*(-1.));
		}
        
		
		
		D4T<std::complex<double>> Bmunu(Amunu.CT(Amunu));
        
		(*Wmunu[k])+=Bmunu;
	}
	
	
    gprime1 = gprime_t(kPi*kPi,kk*kk);
	double vt = pow(gA/2./fpi,2)*(crho*f_rho2*kk*kk/onShell +gprime1);
	
	//dressing of the absorbed rho: ONLY NUCLEON IN THE NOMINATOR, nucleon and Delta in the RPA denominator
	auto rho_inc = [&](const std::complex<double>& Nabs, const std::complex<double>& Ntot, const std::complex<double>& Dtot) -> double {
		std::complex<double> utot=Ntot+Dtot;
		double rhoInc = Nabs.imag();
		rhoInc/=pow(std::abs(1.-vt*utot),2);
		return rhoInc;
	};
	//charged pion (slots 0 and 1) and neutral pion absorbed on a proton or on a neutron. The pi0 of slot 0 is
	//absorbed on the nucleon of its final state N2', the pi0 of slot 2 on the other one.
	double rhoInc_ch = rho_inc(pb.N_ch, pb.N_ch, pb.D_ch);
	double rhoInc_0[2] = {rho_inc(pb.N_0abs[PROT], pb.N_0, pb.D_0), rho_inc(pb.N_0abs[NEUT], pb.N_0, pb.D_0)};
	Nucleon abs0 = vb[0].part;
	Nucleon abs2 = mirror(abs0);
	
	//each slot has its own first bubble: Lindhard function and hole-energy range
	auto cf = [&](int k, double rhoInc) -> double {
		if (!open[k]) return 0.;
		return -pow(f_rho2/onShell,2)*pow(gA/2./fpi,2)*rhoInc*lindhard[k]*crho*kk*kk
		
		/M/M
		
		*kk*kk/pow(2*M_PI,4)/2.*sqrt(1-cosi*cosi)/(Emax[k]-Emin[k])/(2.*M_PI); //the last /2 is from trace 2Pi from averaging over cosTheta_nuc
	};
	//dim: 2-2+2  -2 +2 + 2 = 4
	double coeffA=cf(0,rhoInc_ch);          //ch 1: slot 0
	double coeffB=cf(1,rhoInc_ch);          //ch 2: slot 1
	double coeffC=cf(2,rhoInc_0[abs0]);     //ch 3, pi0 absorbed on abs0: slot 0
	double coeffD=cf(2,rhoInc_0[abs2]);     //ch 3, pi0 absorbed on abs2: slot 2
	
	//slot 0 is W_1 + W_3 (pi0 absorbed on a proton): at Z=N the pi0 splits in two equal halves
	auto slot0 = [&](int i, int j) -> std::complex<double> {
		if (coeffA==2.*coeffC) return coeffA*(Wmunu_1(i,j)+Wmunu_3(i,j)/2.);
		return coeffA*Wmunu_1(i,j)+coeffC*Wmunu_3(i,j);
	};

	A00[0]=slot0(0,0);
	A11[0]=slot0(1,1);
	A03[0]=slot0(0,3);
	A12[0]=slot0(1,2);
	A33[0]=slot0(3,3);
	
	A00[1]=coeffB*Wmunu_2(0,0);
	A11[1]=coeffB*Wmunu_2(1,1);
	A03[1]=coeffB*Wmunu_2(0,3);
	A12[1]=coeffB*Wmunu_2(1,2);
	A33[1]=coeffB*Wmunu_2(3,3);
	
	A00[2]=coeffD*Wmunu_3(0,0);
	A11[2]=coeffD*Wmunu_3(1,1);
	A03[2]=coeffD*Wmunu_3(0,3);
	A12[2]=coeffD*Wmunu_3(1,2);
	A33[2]=coeffD*Wmunu_3(3,3);
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

