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

void hadron_2p2h_part1_Epi_thetapi_phipi_kpi_phip_Ep_isospin(double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, double rho, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33) {
	
	double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);
	
	std::complex<double> cunuc(0.,0.);
	std::complex<double> udel(0.,0.);
	std::complex<double> utot(0.,0.);
	
	ulind(Epi/mpi, kk/mpi, k_fermi/mpi,cunuc, udel);
	cunuc*=pow(mpi,2);//*fpNN2/mpi/mpi*gprime;
	udel*=pow(mpi,2);//*fpNN2/mpi/mpi*gprime;
	
	udel=delta_lind (Epi*fermi,kk*fermi, rho*fermi3, k_fermi*fermi);
	udel/=(fermi*fermi);
	
	
	//INTEGRATION PARAMETERS
	double cosi=cos(thetapi);
	
	double mom_nuc=sqrt(3./5.)*k_fermi;
	double cos_theta_nuc=0; //USUNAC
	double sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
	
	D4V<double> p(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,0,mom_nuc*sin_theta_nuc,mom_nuc*cos_theta_nuc); //nukleon
	D4V<double> p_med(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,0,mom_nuc*sin_theta_nuc,mom_nuc*cos_theta_nuc); //nukleon
	D4V<double> q(w,0,0,qq);
	D4V<double> kPi(Epi, kk*sqrt(1-cosi*cosi)*cos(phi_pi), kk*sqrt(1-cosi*cosi)*sin(phi_pi), kk*cosi);
	
	
	if (std::abs(cunuc.imag())<pow(10,-20)) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	double lindhard = DUR(q(0)-kPi(0), sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2)),k_fermi, k_fermi);
	double q_minus_k = sqrt(pow(kPi(1),2)+pow(kPi(2),2)+pow(kPi(3)-qq,2));
	
	if (std::abs(lindhard)<pow(10,-20)) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		return;}
	
	double Emax=sqrt(M*M+k_fermi*k_fermi);
	double pom = std::max(M,Emax-(w-Epi));
	double Emin=std::max(pom, 0.5*(-(w-Epi)+q_minus_k*sqrt(1.-4.*M*M/((w-Epi)*(w-Epi)-q_minus_k*q_minus_k))));
	
	

	double aver_en_1st = Ep;
	mom_nuc = sqrt(aver_en_1st*aver_en_1st-M*M);
	cos_theta_nuc = ( pow(q(0)-kPi(0),2)-q_minus_k*q_minus_k + 2.*aver_en_1st*(q(0)-kPi(0)) )/(2.*mom_nuc*q_minus_k);
	
	if (std::abs(cos_theta_nuc)>1) {
		A00[0]=A11[0]=A12[0]=A03[0]=A33[0]=0;
		A00[1]=A11[1]=A12[1]=A03[1]=A33[1]=0;
		A00[2]=A11[2]=A12[2]=A03[2]=A33[2]=0;
		std::cerr<<"Error DIR! CosTh="<<cos_theta_nuc<<" w="<<w<<" q="<<qq<<" Epi="<<Epi<<" kpi="<<kk<<" phi_pi="<<phi_pi<<" theta_pi="<<thetapi<<" Ep="<<Ep <<" Lind="<<lindhard<<" Emin,max="<<Emin<<", "<<Emax<<"q-k="<<q_minus_k<<" w-Epi="<<w-Epi<<"\n";
		return;
	}
	sin_theta_nuc=sqrt(1.-pow(cos_theta_nuc,2));
    p(0)=aver_en_1st-k_fermi*k_fermi/2./M*FERMI;
	p(1)=mom_nuc*sin_theta_nuc*cos(phi_p);
	p(2)=mom_nuc*sin_theta_nuc*sin(phi_p);
	p(3)=mom_nuc*cos_theta_nuc; //nukleon
	rotate(p,q-kPi);
	
	
	
	
	p_med(0)=aver_en_1st-k_fermi*k_fermi/2./M*FERMI;
	p_med(1)=mom_nuc*sin_theta_nuc*cos(phi_p);
	p_med(2)=mom_nuc*sin_theta_nuc*sin(phi_p);
	p_med(3)=mom_nuc*cos_theta_nuc; //nucleon
	rotate(p_med,q-kPi);
	
	//	}
	//	std::cerr<<"CosTh="<<cos_theta_nuc<<" E="<<Ep<<"\n";
	
	//    D4V<double> kPi(x[0], 0,0,0);
	D4V<double> pprime=p+q-kPi;
	
	
	DM pslash=p*dirac;
	DM pslash_med=p_med*dirac;
	
	DM qslash=q*dirac;
	
	channel = 0; //we calculate everything with CG coefficients = 1
	D4V<DM> delta_cur_vec, delta_cur_ax, delta_cur_cross_vec, delta_cur_cross_ax, nucl_cur_vec, nucl_cur_ax,nucl_cur_s,  nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s, vector_contact, axial_contact, pion_pole_cur, pion_in_flight_cur;
	/**
	 * Delta current in medium and crossed. This part is the one with resonance
	 */
    delta_current_medium(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_vec, delta_cur_ax);
	delta_current_crossed_medium(rho, p_med, q, kPi, pslash_med, qslash, delta_cur_cross_vec, delta_cur_cross_ax);
	
	pion_pole(q,kPi,qslash, pion_pole_cur);
	contact_term(q,kPi, vector_contact, axial_contact);
	pion_in_flight(q,kPi, pion_in_flight_cur);

	
	D4V<double> p_inter=p+q;
	double p_inter_mom = sqrt(p_inter(1)*p_inter(1)+p_inter(2)*p_inter(2)+p_inter(3)*p_inter(3));
	
	if (p_inter_mom>k_fermi) {
		nucleon_cross_current_medium(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_cross_vec, nucl_cur_cross_ax, nucl_cur_cross_s,k_fermi*k_fermi/2./M*FERMI);
		nucleon_current_medium(kPi, p_med, q, pslash_med, qslash, vector_nuc, axial_nuc, scalar_nuc, nucl_cur_vec, nucl_cur_ax, nucl_cur_s,k_fermi*k_fermi/2./M*FERMI);
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
 
	
	DM pprimeslash=pprime*dirac+M_N*identi;
	
	
	
	D4V<DM> j1,j2, np1, np2, del1, del2, pf1, pf2;
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
		
		j1 =(CGDP_v[ch]/CGDP_v[0])*delta_cur_vec + (CGDP_a[ch]/CGDP_a[0])*delta_cur_ax +(CGDPC_v[ch]/CGDPC_v[0])*delta_cur_cross_vec + (CGDPC_a[ch]/CGDPC_a[0])*delta_cur_cross_ax; // Delta added
		if (p_inter_mom>k_fermi)
            j1+=(CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax +(CGNPC_v[ch]/CGNPC_v[0])*nucl_cur_cross_vec+(CGNPC_a[ch]/CGNPC_a[0])*nucl_cur_cross_ax+(CGNPC_s[ch]/CGNPC_s[0])*nucl_cur_cross_s+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s;
		// Nucleon part added
		j1*=eff_long;
        
        
		j1+=pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]) +(CGCT_v[ch]/CGCT_v[0])*vector_contact + axial_contact*(CGCT_a[ch]/CGCT_a[0]) +  pion_pole_cur*(CGPP_a[ch]/CGPP_a[0]);
		
		//j1=pion_in_flight_cur*(CGPF_v[ch]/CGPF_v[0]);
		
		j2 = j1;
		
		j1 &=pprimeslash;
		j1 *=(pslash+M_N*identi);
		
		
		j2.hermit();
		j2 &=gamma0;
		j2 *=gamma0;
		
		D4T<DM> Amunu(j1,j2); // Here you construct the full Amunu
	
		 pf1=(CGPF_v[ch]/CGPF_v[0])*pion_in_flight_cur;
		 pf2=pf1;
		 pf1 &=pprimeslash;
		 pf1 *=(pslash+M_N*identi);
		 pf2.hermit();
		 pf2 &=gamma0;
		 pf2 *=gamma0;
		 //Amunu.Add(pf1*(-0.5),pf2);
         Amunu.Add(pf1*(-1.),pf2); //I remove 0.5 here and another 0.5 for 2 separate bubbles - in total -1

		
		//now we subtract NP-NP
		
		if (CGNP_v[ch]!=0 && CGNP_a[ch]!=0  && p_inter_mom>k_fermi) {
            np1= (CGNP_v[ch]/CGNP_v[0])*nucl_cur_vec + (CGNP_a[ch]/CGNP_a[0])*nucl_cur_ax+(CGNP_s[ch]/CGNP_s[0])*nucl_cur_s;
			np1*=eff_long;
			np2=np1;
			np1 &=pprimeslash;
			np1 *=(pslash+M_N*identi);
			np2.hermit();
			np2 &=gamma0;
			np2 *=gamma0;
			
			Amunu.Add(np1*(-1.),np2);
		}
		
		
		
		D4T<std::complex<double>> Bmunu(Amunu.CT(Amunu));
		if (ch==min_ch) Wmunu_1+=Bmunu;
		else if (ch==min_ch+1) Wmunu_2+=Bmunu;
		else  Wmunu_3+=Bmunu;
	}
	
	
	
	
	utot=cunuc+udel;
	
	
	double pi_fact = pow(f_pi,2)*pow(gA/2./fpi,2)*kk*kk;
    gprime1 = gprime_l(kPi*kPi,kk*kk);
	double factg = pow(gA/2./fpi,2)*gprime1;
	
	
	std::complex<double> pionImSelf = pi_fact*utot/(1.-factg*utot);
	std::complex<double> pionImSelf_nuc = pi_fact*cunuc;
	
	double pionInc = pionImSelf_nuc.imag()/pow(std::abs(1.-factg*utot),2)/(pow(kPi*kPi-mpi*mpi-pionImSelf.real(),2)+pow(pionImSelf.imag(),2)); //ONLY NUCLEON IN THE NOMINATOR

	
	
	double coeff =pow(f_pi,2)*pionInc*lindhard
	
	/M/M
	
	*kk*kk/pow(2*M_PI,4)/2.*sqrt(1-cosi*cosi)/(Emax-Emin)/(2.*M_PI); // angle /2 from Amunu trace
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
	
	std::complex<double> z, zp, pzeta, pzetap, i=std::complex<double>(0,1);
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
	
	
	
	return 2./3. * rho * md/(q_mod*k_fermi) * (pzeta +pzetap) * pow(fdel_f,2);
	
	
}
