//
//  hadron_tensor_2p2h.cc
//  Pions
//
//  Created by Joanna Sobczyk on 20/01/16.
//  Copyright (c) 2016 Joanna Sobczyk. All rights reserved.
//

/*
 CHANGES:
 - Delta CA5 form factor
 
 
 */

#include "hadron_tensor_2p2h.h"


bool Global = false;
bool FERMI=false; //separation energy
int NR_PART=12;
CHANNEL CH_2p2h=CC;




/**********************
 
 
 FULL VERSION WITH SPLITTING INTO isospin
 
 
 **********************/

void hadron_2p2h_part1_Epi_thetapi_phipi_kpi_phip_Ep_isospin(double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, const FermiGas& fg, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33) {
	
	const double rho = fg.rho;
	
	//Lindhard functions of the pion that is absorbed (second bubble), for its charge
	PionBubbles pb;
	pion_bubbles(fg, Epi, kk, pb, true);
	
	
	//INTEGRATION PARAMETERS
	double cosi=cos(thetapi);
	
	double cos_theta_nuc=0; //USUNAC
	double sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	D4V<double> q(w,0,0,qq);
	D4V<double> kPi(Epi, kk*sqrt(1-cosi*cosi)*cos(phi_pi), kk*sqrt(1-cosi*cosi)*sin(phi_pi), kk*cosi);
	
	
	if (!pion_bubbles_open(pb, pow(10,-20))) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	double q_minus_k = sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2));
	
	//the three isospin slots have their own bubble (hole and particle species, hence their own Fermi momenta)
	BubbleSpecies vb[3];
	double kF_hole[3], lindhard[3], Emin[3], Emax[3];
	bool open[3];
	int first_open=-1;
	for (int k=0; k<3; k++) {
		vb[k] = vertex_bubble(k);
		kF_hole[k] = kF_of(fg, vb[k].hole);
		lindhard[k] = DUR(q(0)-kPi(0), sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2)), kF_hole[k], kF_of(fg, vb[k].part));
		bubble_limits(q(0)-kPi(0), q_minus_k, kF_hole[k], kF_of(fg, vb[k].part), Emin[k], Emax[k]);
		open[k] = !(std::abs(lindhard[k])<pow(10,-20)) && !(Ep<Emin[k]) && !(Ep>Emax[k]);
		if (open[k] && first_open<0) first_open=k;
	}
	
	if (first_open<0) {
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
		std::cerr<<"Error DIR! CosTh="<<cos_theta_nuc<<" w="<<w<<" q="<<qq<<" Epi="<<Epi<<" kpi="<<kk<<" phi_pi="<<phi_pi<<" theta_pi="<<thetapi<<" Ep="<<Ep <<" Lind="<<lindhard[first_open]<<" Emin,max="<<Emin[first_open]<<", "<<Emax[first_open]<<"q-k="<<q_minus_k<<" w-Epi="<<w-Epi<<"\n";
		return;
	}
	sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	channel = 0; //we calculate everything with CG coefficients = 1
	D4V<DM> vector_contact, axial_contact, pion_pole_cur, pion_in_flight_cur;
	DM qslash=q*dirac;
	
	pion_pole(q,kPi,qslash, pion_pole_cur);
	contact_term(q,kPi, vector_contact, axial_contact);
	pion_in_flight(q,kPi, pion_in_flight_cur);
	
	//the Pauli blocking of the intermediate nucleon of the direct and of the crossed nucleon pole depends on its species
	bool nucl_direct=false, nucl_crossed=false;
	
	//currents of the W vertex on a hole; they depend on the species of the hole through its removal-energy shift only
	struct HoleCurrents {
		DM pslash, pslash_med, pprimeslash;
		D4V<DM> delta_cur_vec, delta_cur_ax, delta_cur_cross_vec, delta_cur_cross_ax, nucl_cur_vec, nucl_cur_ax, nucl_cur_s, nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s;
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
		rotate(p,q-kPi);
		D4V<double> p_med(aver_en_1st-shift,mom_nuc*sin_theta_nuc*cos(phi_p),mom_nuc*sin_theta_nuc*sin(phi_p),mom_nuc*cos_theta_nuc); //nucleon
		rotate(p_med,q-kPi);
		
		D4V<double> pprime=p+q-kPi;
		H.pslash=p*dirac;
		H.pslash_med=p_med*dirac;
		H.pprimeslash=pprime*dirac+M_N*identi;
		
		/**
		 * Delta current in medium and crossed. This part is the one with resonance
		 */
		delta_current_medium(rho, p_med, q, kPi, H.pslash_med, qslash, H.delta_cur_vec, H.delta_cur_ax);
		delta_current_crossed_medium(rho, p_med, q, kPi, H.pslash_med, qslash, H.delta_cur_cross_vec, H.delta_cur_cross_ax);
		
		D4V<double> p_inter=p+q;
		double p_inter_mom = sqrt(p_inter(1)*p_inter(1)+p_inter(2)*p_inter(2)+p_inter(3)*p_inter(3));
		nucl_direct  = p_inter_mom>kF_of(fg, intermediate_direct());
		nucl_crossed = p_inter_mom>kF_of(fg, intermediate_crossed());
		
		if (nucl_crossed)
			nucleon_cross_current_medium(kPi, p_med, q, H.pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, H.nucl_cur_cross_vec, H.nucl_cur_cross_ax, H.nucl_cur_cross_s,shift);
		if (nucl_direct)
			nucleon_current_medium(kPi, p_med, q, H.pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, H.nucl_cur_vec, H.nucl_cur_ax, H.nucl_cur_s,shift);
	}
	
	
	double f_pi = (lambda_pi*lambda_pi-mpi*mpi)/(lambda_pi*lambda_pi-kPi*kPi);
	double onShell = kPi*kPi-mpi*mpi;
	
	//CHARGE CURRENT
	//calculate 3 different isospin channels:
	//  n --> p + pi0  channel = 3
	//  n --> n + pi+  channel = 2
	//  p --> p + pi+  channel = 1
	
	//constant by which we multiply Delta, DeltaCrossed, Nucleon, NucleonCrossed currents
	
	double eff_long = 1;
    
    double gprime1 = gprime_l(kPi*kPi,kk*kk);
    eff_long=1.+ gprime1*onShell/((kk*kk)*f_pi*f_pi);
 
	
	
	D4V<DM> j1,j2, np1, np2, del1, del2, pf1, pf2;
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
		
		j1 =(CGDP_v[ch]/CGDP_v[0])*H.delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*H.delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*H.delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*H.delta_cur_cross_ax; // Delta added
		if (nucl_direct && nucl_crossed)
            j1+=(CGNP_v[ch]/CGNP_v[0])*H.nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*H.nucl_cur_ax +(CGNPC_v[ch]/CGNPC_v[0])*H.nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*H.nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*H.nucl_cur_cross_s+(CGNP_s[ch]/CGNP_s[0])*H.nucl_cur_s;
		else if (nucl_direct)
            j1+=(CGNP_v[ch]/CGNP_v[0])*H.nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*H.nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*H.nucl_cur_s;
		else if (nucl_crossed)
            j1+=(CGNPC_v[ch]/CGNPC_v[0])*H.nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*H.nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*H.nucl_cur_cross_s;
		// Nucleon part added
		j1*=eff_long;
        
        
		j1+=pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]) +(CGCT_v[ch]/CGCT_v[0])*vector_contact + axial_contact*(CGCT_a[ch]/CGCT_a[0]) +  pion_pole_cur*(CGPP_a[ch]/CGPP_a[0]);
		
		//j1=pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]);
		
		j2 = j1;
		
		j1 &=H.pprimeslash;
		j1 *=(H.pslash+M_N*identi);
		
		
		j2.hermit();
		j2 &=gamma0;
		j2 *=gamma0;
		
		D4T<DM> Amunu(j1,j2); // Here you construct the full Amunu
	
		 pf1=(CGPF_v[ch]/CGPF_v[0])*pion_in_flight_cur;
		 pf2=pf1;
		 pf1 &=H.pprimeslash;
		 pf1 *=(H.pslash+M_N*identi);
		 pf2.hermit();
		 pf2 &=gamma0;
		 pf2 *=gamma0;
		 //Amunu.Add(pf1*(-0.5),pf2);
         Amunu.Add(pf1*(-1.),pf2); //I remove 0.5 here and another 0.5 for 2 separate bubbles - in total -1

		
		//now we subtract NP-NP
		
		if (CGNP_v[ch]!=0 && CGNP_a[ch]!=0  && nucl_direct) {
            np1= (CGNP_v[ch]/CGNP_v[0])*H.nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*H.nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*H.nucl_cur_s;
			np1*=eff_long;
			np2=np1;
			np1 &=H.pprimeslash;
			np1 *=(H.pslash+M_N*identi);
			np2.hermit();
			np2 &=gamma0;
			np2 *=gamma0;
			
			Amunu.Add(np1*(-1.),np2);
		}
		
		
		
		D4T<std::complex<double>> Bmunu(Amunu.CT(Amunu));
		(*Wmunu[k])+=Bmunu;
	}
	
	
	
	
	double pi_fact = pow(f_pi,2)*pow(gA/2./fpi,2)*kk*kk;
    gprime1 = gprime_l(kPi*kPi,kk*kk);
	double factg = pow(gA/2./fpi,2)*gprime1;
	
	//dressing of the absorbed pion: ONLY NUCLEON IN THE NOMINATOR, nucleon and Delta in the RPA denominator
	auto pion_inc = [&](const std::complex<double>& Nabs, const std::complex<double>& Ntot, const std::complex<double>& Dtot) -> double {
		std::complex<double> utot=Ntot+Dtot;
		std::complex<double> pionImSelf = pi_fact*utot/(1.-factg*utot);
		std::complex<double> pionImSelf_nuc = pi_fact*Nabs;
		return pionImSelf_nuc.imag()/pow(std::abs(1.-factg*utot),2)/(pow(kPi*kPi-mpi*mpi-pionImSelf.real(),2)+pow(pionImSelf.imag(),2));
	};
	//charged pion (slots 0 and 1) and neutral pion absorbed on a proton or on a neutron. The pi0 of slot 0 is
	//absorbed on the nucleon of its final state N2', the pi0 of slot 2 on the other one.
	double pionInc_ch = pion_inc(pb.N_ch, pb.N_ch, pb.D_ch);
	double pionInc_0[2] = {pion_inc(pb.N_0abs[PROT], pb.N_0, pb.D_0), pion_inc(pb.N_0abs[NEUT], pb.N_0, pb.D_0)};
	Nucleon abs0 = vb[0].part;
	Nucleon abs2 = mirror(abs0);
	
	//each slot has its own first bubble: Lindhard function and hole-energy range
	auto cf = [&](int k, double pionInc) -> double {
		if (!open[k]) return 0.;
		return pow(f_pi,2)*pionInc*lindhard[k]
		
		/M/M
		
		*kk*kk/pow(2*M_PI,4)/2.*sqrt(1-cosi*cosi)/(Emax[k]-Emin[k])/(2.*M_PI); // angle /2 from Amunu trace
	};
	//dim: 2-2+2  -2 +2 + 2 = 4
	double coeffA=cf(0,pionInc_ch);          //ch 1: slot 0
	double coeffB=cf(1,pionInc_ch);          //ch 2: slot 1
	double coeffC=cf(2,pionInc_0[abs0]);     //ch 3, pi0 absorbed on abs0: slot 0
	double coeffD=cf(2,pionInc_0[abs2]);     //ch 3, pi0 absorbed on abs2: slot 2
	
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
 
 
 LINDHARD FUNCTIONS
 
 
 **********************/



//Nonrelativistic Lindhard Functions
//sum over isospin and spin of nucleon
void ulind(double QZR,double Q, double XKF,std::complex<double> &CUFUN, std::complex<double> &cudel){
    
    
    double XMN,XMD,WRES,FNS,FDS,QA,QZ,QZA,RUN,AM,AP,TERF,TERS,YMN, DCOCI;
    double B,RO,FAC,A,TDIR,TCROS,SQS,QFR,GAMH;
    std::complex<double> CYI,CCROS,CDIR,CAC,CAPC,cunuc;
    
    
    XMN=M/mpi;//6.7348*MeV;   //nucleon mass (in units of mpi)  mn/mpi
    XMD=Md/mpi;//8.8315*MeV;   //delta mass        " "
    WRES=XMD-XMN;
    FNS=1.;
    FDS=4.52;
    CYI=std::complex<double>(0.,1.);
    
    QA=Q/XKF;
    QZ=std::abs(QZR);
    QZA=std::abs(QZ*XMN/(XKF*XKF));
    
    AM=QZA/QA-QA/2.;
    AP=QZA/QA+QA/2.;
    
//    std::cerr<<QZR<<" "<<Q<<" "<<XKF<<" "<<AM<<" "<<AP<<"\n";
    
    if(Q==0.)
        RUN=FNS*XMN*XKF/(M_PI*M_PI)*2./3.*pow(QA/QZA,2);
    
    else {
        AM=QZA/QA-QA/2.;
        AP=QZA/QA+QA/2.;
        if (std::abs(1.-AM)< 0.00001) AM=1.00001;
        if (std::abs(1.-AP)<0.00001) AP=1.00001;
        TERF=std::abs((1.+AM)/(1.-AM))+1.e-15;
        TERS=std::abs((1.+AP)/(1.-AP))+1.e-15;
  //      std::cerr<<"real "<<TERF<<" "<<TERS<<" "<<AM<<" "<<AP<<"\n";
        if (QZA-0.00001<0.) {
            RUN=FNS*XMN*XKF/(M_PI*M_PI)*(-1.+(1.-AM*AM)/(2.*QA)*std::log(TERF)-(1.-AP*AP)/(2.*QA)*std::log(TERS));
        }
        else {
            DCOCI=std::abs(QA/QZA);
            
            if (DCOCI-0.1<0.) {
                RUN=FNS*XMN*XKF/(M_PI*M_PI)*2./3.*pow(QA/QZA,2);
            }
            else {
                RUN=FNS*XMN*XKF/(M_PI*M_PI)*(-1.+(1.-AM*AM)/(2. *QA)*std::log(TERF)-(1.-AP*AP)/(2.*QA)*std::log(TERS));
            }
        }
        
    }
    
    
    YMN=0.;
    
    if((QA*QA/2.+QA)>=QZA && QZA>=(QA*QA/2.-QA)  && QA>=2.){
        YMN=-2.*XMN*XKF*FNS/(4.*M_PI*QA)*(1.-AM*AM);}
    
    else if (QA<2. && (QA+QA*QA/2.)>=QZA && QZA >=(QA-QA*QA/2.)){
        YMN=-2.*XMN*XKF*FNS/(4.*M_PI*QA)*(1.-AM*AM);}
    
    else if (QA<2. && 0.<QZA && QZA<(QA-QA*QA/2.)){
        YMN=-2.*XMN*XKF*FNS/(4.*M_PI*QA)*2.*QZA;}
    
 //   std::cerr<<"ymn="<<YMN<<"\n";
    cunuc=std::complex<double>(RUN,YMN);
    
    B=Q/XMD;
    RO=pow(XKF,3)*2./(3.*M_PI*M_PI);
    
    if(Q==0.){
        FAC=0.;
    }
    else {
        FAC=pow((4./3.*XKF/(2.*M_PI)),2)/pow(B,3)*FDS;}
    
    if (std::abs(QZ)-1.<=0.) {
        A=(QZ-WRES-Q*Q/(2.*XMD))/XKF;
        AP=(-QZ-WRES-Q*Q/(2.*XMD))/XKF;
        
        if(std::abs(B/A)-0.1<=0.){
            TDIR=4./9.*FDS*RO/(A*XKF);}
        
        else {
            TERF=std::abs((A+B)/(A-B))+1.e-15;
            TDIR=FAC*(B*A+(B*B-A*A)/2.*std::log(TERF));}
        
        
        if(std::abs(B/AP)-0.1<=0.){
            TCROS=4./9.*FDS*RO/(AP*XKF);}
        else {
            TERS=std::abs((AP-B)/(AP+B))+1.e-15;
            TCROS=FAC*(B*AP-(B*B-AP*AP)/2.*std::log(TERS));}
        
        cudel=TDIR+TCROS;
    }
    else {
        SQS=XMN+std::abs(QZ);
        QFR=sqrt((QZ*QZ-1.)/(1.+QZ/XMN));
        GAMH=1./(3.*4.*M_PI)*FDS*XMN/SQS*pow(QFR,3);
        
        if (QZ>0.) {
            CAC=(QZ-WRES+CYI*GAMH-Q*Q/(2.*XMD))/XKF;
            CAPC=(-QZ-WRES-Q*Q/(2.*XMD))/XKF;
        }
        else if(QZ<=0.) {
            CAC=(QZ-WRES-Q*Q/(2.*XMD))/XKF;
            CAPC=(-QZ-WRES+CYI*GAMH-Q*Q/(2.*XMD))/XKF;
        }
        
        if(std::abs(B/CAC)-0.1<=0) {
            CDIR=4./9.*FDS*RO/(CAC*XKF);}
        else {
            CDIR=FAC*(B*CAC+(B*B-CAC*CAC)/2.*std::log((CAC+B)/(CAC-B)));
        }
        
        if(std::abs(B/CAPC)-0.1<=0) {
            CCROS=4./9.*FDS*RO/(CAPC*XKF);}
        else {
            CCROS=FAC*(B*CAPC-(B*B-CAPC*CAPC)/2.*std::log((CAPC-B)/(CAPC+B)));
        }
        cudel=CDIR+CCROS;
    }
    CUFUN=cunuc;
}


//Relativistic Im Lindhard
//not summed over isospin
double DUR (double q0, double dq, double pf1, double pf2){
    
    double CTE=-0.5*M*M/M_PI/dq;
    double q2=q0*q0-dq*dq;
    if (q2>=0) {
        //          std::cerr<<"DUR=0 dq="<<dq*197.33<<" q0="<<q0*197.33<<"\n";
        return 0.;
    }
    
    double ef1=sqrt(M*M+pf1*pf1);
    double ef2=sqrt(M*M+pf2*pf2);
    if ((ef2-q0)>ef1) {
        return 0.;
    }
    
    double a=0.5*(-q0+dq*sqrt(1.-4.*M*M/q2));
    //   std::cerr<<"a="<<a<<" ef1="<<ef1<<" ef2="<<ef2<<" q0, dq="<<q0<<" "<<dq<<"\n";
    
    if (a>ef1) {
        //             std::cerr<<"DUR=0 a>ef1 dq="<<dq*197.33<<" q0="<<q0*197.33<<" pf="<<pf1*197.33<<"\n";
        return 0.;
    }
    double pom = std::max(M,ef2-q0);
    double einf=std::max(pom,a);

    //     std::cerr<<"DUR "<<"pf="<<pf1*197<<" ef1="<<ef1<<" einf="<<einf<<" CTe="<<CTE<<"\n";
    return (ef1-einf)*CTE;
}



//Direct (pzeta) and crossed (pzetap) Delta-hole functions of delta_lind for the Fermi momentum k_fermi
static void delta_zetas(std::complex<double> q_zero, double q_mod, double k_fermi, std::complex<double>& pzeta, std::complex<double>& pzetap) {
	std::complex<double> z, zp, i=std::complex<double>(0,1);
	double fdel_f = fstar;
	double rq_zero, s, sp,srot, srotp, qcm, qcmp;
	double md=6.2433, m=4.7592, mpi=0.7045;
	double wr = md-m;
	double gamma = 0.;
	double gammap = 0.;
	if(pow(q_zero.imag(),2)<1e-36){
		rq_zero = q_zero.real();
		s = m*m+pow(rq_zero,2)-pow(q_mod,2)+  2.*rq_zero *sqrt(m*m+3./5.*pow(k_fermi,2));
		if(s>pow(m+mpi,2)){
			srot = sqrt(s);
			qcm = sqrt(s*s+pow(mpi,4)+pow(m,4)-2.*(s*mpi*mpi+s*m*m+ pow(mpi*m,2))) /(2.*srot);
			gamma = 1./3. * 1./(4.*M_PI) * pow(fdel_f,2)*pow(qcm,3)/srot*(m+sqrt(m*m+qcm*qcm))/pow(mpi,2);
		}
		sp = m*m+pow(rq_zero,2)-pow(q_mod,2)- 2.*rq_zero *sqrt(m*m+3./5.*pow(k_fermi,2));
		if(sp>pow(m+mpi,2)){
			srotp = sqrt(sp);
			qcmp=sqrt(sp*sp+pow(mpi,4)+pow(m,4)-2.*(sp*pow(mpi,2)+sp*m*m+pow(mpi*m,2)))/(2.*srotp);
			gammap = 1./3. * 1./(4.*M_PI) * pow(fdel_f,2)* pow(qcmp,3)/srotp*(m+sqrt(m*m+qcmp*qcmp))/pow(mpi,2);
		}
	}
	
	z=md/(q_mod*k_fermi)*(q_zero-pow(q_mod,2)/(2.*md) -wr +i*gamma/2.);
	
	zp=md/(q_mod*k_fermi)*(-q_zero-pow(q_mod,2)/(2.*md)-wr +i*gammap/2.);
	
	if(sqrt(pow(z.real(),2)+pow(z.imag(),2))>50.){
		pzeta =  2./(3.*z) +2./(15.*pow(z,3));
	}
	else if(sqrt(pow(z.real(),2)+pow(z.imag(),2))<1e-2){
		pzeta =  2.*z -2./3.*pow(z,3) -i*M_PI/2.*(1. - pow(z,2));
	}
	else
	pzeta =  z +  (1.-z*z) * log((z+1.)/(z-1.))/2.;
	
	
	if(sqrt(pow(zp.real(),2)+pow(zp.imag(),2))>50.)
	pzetap =  2./(3.*zp) +2./(15.*pow(zp,3));
	else if(sqrt(pow(zp.real(),2)+pow(zp.imag(),2))<1e-2)
	pzetap =  2.*zp -2./3.*pow(zp,3) -i*M_PI/2.*(1. - pow(zp,2));
	else
	pzetap =  zp + (1.-pow(zp,2)) * log((zp+1.)/(zp-1.))/2.;
	
	
	
}


std::complex<double> delta_lind (std::complex<double> q_zero,double q_mod, double rho,double k_fermi) {
	/*
	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 c   complex Lindhard function for symmetric nuclear matter:
	 c                    from Appendix of
	 c                    E.Oset et al Phys. Rept. 188:79, 1990
	 c                    formula A.4
	 c
	 c   input variables:
	 c     q_zero [fm^-1] : Energy
	 c     q_mod  [fm^-1] : Momentum
	 c     rho    [fm^3]  : Nuclear density
	 c     k_fermi[fm^-1] : Fermi momentum
	 c
	 c   All variables are real*8
	 c
	 c   output variable:
	 c     delta_lind [fm^-2]
	 c
	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 c            ATTENTION!!!
	 c Only works properly for real q_zero,
	 c if q_zero has an imaginary part calculates the L. function
	 c assuming Gamma= 0.
	 c Therefore this subroutine provides two different functions
	 c depending on whether q_zero is real or not!!!!!!!!!!!
	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 
	 implicit none
	 real*8 q_mod,rho,k_fermi,m,rq_zero,gamma,s,srot,wr
	 real*8 fdel_f,mpi,gammap,sp,srotp,qcm,qcmp,md,pi
	 complex*16 delta_lind,q_zero,z,zp,i,pzeta,pzetap
	 
	 real*8 fp0in,fp0ex,f,fstar,pilambda,cr0
	 real*8 rholambda,gp,xmpi,xmrho
	 
	 common/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
	 f       cr0,rholambda,gp,xmpi,xmrho
	 
	 c
	 c m = 939/197.3,md = 1232./197.3, mpi = 139./197.3
	 c
	 data md,m,i/6.2433d0,4.7592d0,(0.d0,1.d0)/
	 data mpi,pi/0.7045d0,3.14159265358979323846264338d0/
	 c
	 c  f*f = 2.13 --> f*^2/4pi = .36
	 c  NOTE: fdel = fdel_f because f =1.0
	 c
	 */
	
	std::complex<double> pzeta, pzetap;
	delta_zetas(q_zero, q_mod, k_fermi, pzeta, pzetap);
	double fdel_f = fstar;
	double md=6.2433;
	
	return 2./3. * rho * md/(q_mod*k_fermi) * (pzeta +pzetap) * pow(fdel_f,2);
	
	
}


/**********************
 
 
 PION DRESSING IN ASYMMETRIC NUCLEI
 
 
 **********************/

//Crossed part of the nonrelativistic nucleon-hole Lindhard function of ulind (same normalisation: sum over spin and
//isospin, units of mpi) for the Fermi momentum XKF. ulind(XKF) = direct(XKF) + crossed(XKF). The sum
//direct(kF_hole)+crossed(kF_part) is the bubble of a hole sea kF_hole and a particle sea kF_part:
//the terms with both states occupied cancel between the direct and the crossed diagram.
std::complex<double> ulind_crossed(double QZR, double Q, double XKF) {
    double XMN=M/mpi;
    double FNS=1.;
    if (Q==0. || XKF<=0.) return std::complex<double>(0.,0.);
    
    double QA=Q/XKF;
    double QZ=std::abs(QZR);
    double QZA=std::abs(QZ*XMN/(XKF*XKF));
    double AP=QZA/QA+QA/2.;
    if (std::abs(1.-AP)<0.00001) AP=1.00001;
    
    double c=FNS*XMN*XKF/(M_PI*M_PI);
    double re;
    if (std::abs(AP)>10.) re=-c/QA*(2./(3.*AP)+2./(15.*pow(AP,3))+2./(35.*pow(AP,5)));
    else re=c*(-AP/QA-(1.-AP*AP)/(2.*QA)*std::log(std::abs((1.+AP)/(1.-AP))+1.e-15));
    
    double im=0.;
    if (QZA>0. && AP<1.) im=2.*XMN*XKF*FNS/(4.*M_PI*QA)*(1.-AP*AP);
    return std::complex<double>(re,im);
}

//Direct and crossed parts of delta_lind (it is their sum), with the same arguments
void delta_lind_pieces(std::complex<double> q_zero,double q_mod, double rho,double k_fermi, std::complex<double>& direct, std::complex<double>& crossed) {
    std::complex<double> pzeta, pzetap;
    delta_zetas(q_zero, q_mod, k_fermi, pzeta, pzetap);
    double fdel_f = fstar;
    double md=6.2433;
    direct  = 2./3. * rho * md/(q_mod*k_fermi) * pzeta * pow(fdel_f,2);
    crossed = 2./3. * rho * md/(q_mod*k_fermi) * pzetap * pow(fdel_f,2);
}

/**
 * Nucleon-hole and Delta-hole Lindhard functions of the pion absorbed in the one-body diagrams, for its charge.
 * The nucleon-hole bubble of a charged pion is 2*(g=2 Lindhard function): isospin factor 2 of tau_+- times the bubble
 * with the hole sea of the nucleon that absorbs it and the particle sea of the nucleon it becomes. The pi0 is absorbed
 * in a pp or a nn bubble with isospin factor 1. The Delta-hole function weights the species with the squared
 * Clebsch-Gordan coefficients of N pi -> Delta; the crossed term has the weights of the opposite charge.
 * For Z=N all of them reduce to the isospin-averaged cunuc and udel that were used before.
 */
void pion_bubbles(const FermiGas& fg, double Epi, double kk, PionBubbles& pb, bool withDelta) {
    std::complex<double> cunuc(0.,0.), udel(0.,0.), dummy(0.,0.);
    
    if (fg.kF_p==fg.kF_n) {
        ulind(Epi/mpi, kk/mpi, fg.kF/mpi, cunuc, dummy);
        cunuc*=pow(mpi,2);
        if (withDelta) {
            udel=delta_lind (Epi*fermi,kk*fermi, fg.rho*fermi3, fg.kF*fermi);
            udel/=(fermi*fermi);
        }
        pb.N_ch=pb.N_0=cunuc;
        pb.D_ch=pb.D_0=udel;
        pb.N_0abs[PROT]=pb.N_0abs[NEUT]=cunuc*0.5;
        return;
    }
    
    const double kF[2]={fg.kF_p, fg.kF_n};
    const double rho_s[2]={fg.rho_p, fg.rho_n};
    
    //nucleon-hole Lindhard functions and their crossed parts for each species
    std::complex<double> U[2], C[2];
    for (int sp=0; sp<2; sp++) {
        if (kF[sp]<=0.) { U[sp]=C[sp]=0.; continue; }
        ulind(Epi/mpi, kk/mpi, kF[sp]/mpi, U[sp], dummy);
        U[sp]*=pow(mpi,2);
        C[sp]=ulind_crossed(Epi/mpi, kk/mpi, kF[sp]/mpi)*pow(mpi,2);
    }
    //charged pion of the one-body diagrams: pi+ (n -> p) for neutrinos, pi- (p -> n) for antineutrinos
    BubbleSpecies ch=bubble_species(NEUT,PROT);
    //With FERMI the species potentials are such that the Fermi energies of p and n coincide: a charge-exchange bubble
    //(hole of species h, particle of species p) is evaluated at the energy w - (V_p - V_h) = w - (kF_h^2 - kF_p^2)/2M.
    //The crossed term is the same function at -w, so it takes the same shift. Without FERMI there are no potentials.
    const double dV = FERMI ? (kF[ch.hole]*kF[ch.hole]-kF[ch.part]*kF[ch.part])/2./M : 0.;
    if (dV==0.) {
        pb.N_ch = U[ch.hole] - C[ch.hole] + C[ch.part];
    } else {
        //ulind and ulind_crossed take |w|. The bubble is A_h(w)+A_p(-w) with A_k(w) = direct, A_k(-w) = crossed (U = A(w)+A(-w)), so
        //for a negative shifted energy the roles of direct and crossed are exchanged: A_h(-|w|)+A_p(|w|) = C_h(|w|) + U_p(|w|) - C_p(|w|)
        const double Eeff=Epi-dV;
        const double e=std::abs(Eeff);
        auto Uk=[&](int sp){ std::complex<double> u(0.,0.); if (kF[sp]>0.) { ulind(e/mpi, kk/mpi, kF[sp]/mpi, u, dummy); u*=pow(mpi,2); } return u; };
        auto Ck=[&](int sp){ return kF[sp]>0. ? ulind_crossed(e/mpi, kk/mpi, kF[sp]/mpi)*pow(mpi,2) : std::complex<double>(0.,0.); };
        if (Eeff>=0.) pb.N_ch = Uk(ch.hole) - Ck(ch.hole) + Ck(ch.part);
        else          pb.N_ch = std::conj(Ck(ch.hole) + Uk(ch.part) - Ck(ch.part)); //the Im part of the |w| functions has the sign of w>0: Re is unchanged, Im flips (checked against numerical integration)
    }
    //Without the shift (FERMI off) the Fermi energies of p and n differ, and the direct+crossed form also contains
    //the emission p <-> n of the pion by the Fermi sea with the higher Fermi energy (e.g. n -> p pi- in a neutron-rich
    //nucleus). That is not an absorption: keep the Im part (absorption) negative.
    if (pb.N_ch.imag()>0.) pb.N_ch=std::complex<double>(pb.N_ch.real(),0.);
    pb.N_0abs[PROT]=U[PROT]*0.5;
    pb.N_0abs[NEUT]=U[NEUT]*0.5;
    pb.N_0=pb.N_0abs[PROT]+pb.N_0abs[NEUT];
    
    pb.D_ch=pb.D_0=std::complex<double>(0.,0.);
    if (!withDelta) return;
    
    //species pieces of the Delta-hole function, normalised so that the two species add up to delta_lind at Z=N
    std::complex<double> dir[2], cro[2];
    for (int sp=0; sp<2; sp++) {
        if (kF[sp]<=0.) { dir[sp]=cro[sp]=0.; continue; }
        delta_lind_pieces(Epi*fermi, kk*fermi, 2.*rho_s[sp]*fermi3, kF[sp]*fermi, dir[sp], cro[sp]);
        dir[sp]/=(fermi*fermi);
        cro[sp]/=(fermi*fermi);
    }
    //weights (fractions of 4/3) of the direct and crossed term of pi+, pi0, pi- on (p, n)
    const double wdir[3][2]={{3./4.,1./4.},{1./2.,1./2.},{1./4.,3./4.}};
    const double wcro[3][2]={{1./4.,3./4.},{1./2.,1./2.},{3./4.,1./4.}};
    int ic = CH_2p2h==CCanti ? 2 : 0;
    for (int sp=0; sp<2; sp++) {
        pb.D_ch += wdir[ic][sp]*dir[sp] + wcro[ic][sp]*cro[sp];
        pb.D_0  += wdir[1][sp]*dir[sp]  + wcro[1][sp]*cro[sp];
    }
}
