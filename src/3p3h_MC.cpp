//
//  3p3h.cpp
//  Pions
//
//  Created by Joanna Sobczyk on 12/02/2019.
//  Copyright © 2019 Joanna Sobczyk. All rights reserved.
//

#include "3p3h_MC.h"


/***************
 *
 *
 *  3p3h - MONTE CARLO
 *
 *
 ****************/


void test_events(double Ein, nucleus_str* nuc, double& max_weight, double &total) {
	
	int nr_events=3000000;
	max_weight = 0;
	double rho, r, w, theta_lepton, qq, p_nuc, phi_nuc, theta_nuc, weight, k_fermi, w_max, theta_lepton_max;
	std::complex<double> A00, A11, A12, A03, A33;
	
	if (Global) {
		rho = (getRho(0., nuc, PROTON)+getRho(0., nuc, NEUTRON))/fermi3;
		k_fermi = pow(3*M_PI*M_PI/2.*rho,1./3.);
	}
	
	bool result_zero=false;
	
	double suma = 0;
	
	for (int i = 0; i<nr_events; i++) {
		
        double corr_in_medium = nuc->Q_2p2h;
        w = frandom()*(Ein-lepton_mass-corr_in_medium);
		theta_lepton = frandom()*M_PI;
		
		if (!Global) {
			r = frandom()*nuc->rmax;
			rho = (getRho(r, nuc, PROTON)+getRho(r, nuc, NEUTRON))/fermi3;
			k_fermi = pow(3*M_PI*M_PI/2.*rho,1./3.);
		}
		p_nuc = frandom()*k_fermi;
		theta_nuc = frandom()*M_PI;
		phi_nuc=frandom()*2.*M_PI;
		
		double eout = Ein-w;
		double kout = sqrt(eout*eout-lepton_mass*lepton_mass);
		qq = sqrt(Ein*Ein+kout*kout-2*Ein*kout*cos(theta_lepton));
		
		if (qq<w) result_zero=true;
		//if (qq>1.2*GeV) result_zero=true;
		init_deltaFF(-w*w+qq*qq,w);
		
		if (!result_zero) {
		hadron_delta_p_phi_theta(theta_nuc, phi_nuc, p_nuc, w-nuc->Q_2p2h, qq, rho, true, false, false, A00, A11, A12, A03, A33);
//		std::cerr<<"w, q "<<w<<" "<<qq<<" "<<p_nuc<<"\n";
//		std::cerr<<"ALL a: "<<A00<<" "<<A11<<" "<<A12<<" "<<A03<<" "<<A33<<"\n";
		double sint2 = pow(sin(theta_lepton/2.),2);
		double cost2 = pow(cos(theta_lepton/2.),2);
		
		double W1= A11.real()/2.;
		double W2 =0.5*(A00.real()+A11.real()+ w*w/qq/qq*(A33.real()-A11.real())-2*w/qq*A03.real() );
		double W3 =A12.imag()/qq;
		double W4 =(A33.real()-A11.real())/2./qq/qq;
		double W5 =(A03.real()-w/qq*(A33.real()-A11.real()))/qq;
		
		
		double as = lepton_mass*lepton_mass/eout/(eout+kout);
		
		double sum1 = 2*W1*sint2 + W2*cost2 - W3*(Ein+eout)*sint2;
		double sum2 = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) + W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
		
		
		
		double coef = kout*eout*GF*GF/M_PI/M_PI*sin(theta_lepton);
		
		double factor_rho = NR_PART/rho/2./M_PI;
		if (!Global) factor_rho = 2*r*r*pow(fermi,3);//4*M_PI*r*r/2./M_PI;
		
			if (!Global) weight = coef*(sum1+as*sum2)*factor_rho*2*M_PI*M_PI*k_fermi*M_PI*(Ein-lepton_mass-corr_in_medium)*nuc->rmax*2*M_PI;
			else weight = coef*(sum1+as*sum2)*factor_rho*2*M_PI*M_PI*k_fermi*M_PI*(Ein-lepton_mass-corr_in_medium)*2*M_PI;
		if (weight>max_weight) {
			max_weight = weight;
			w_max = w;
			theta_lepton_max=theta_lepton;
//				std::cerr<<"max values: w="<<w_max<<" Theta="<<theta_lepton_max<<" weight="<<max_weight<< "\n";
		}
			suma+=weight;
		}
		result_zero=false;
	}
	
	total = suma/nr_events*pow(10,38)/cm2;

	std::cerr<<"max values: w="<<w_max<<" Theta="<<theta_lepton_max<<" weight="<<max_weight<<" total="<<total<<"\n";

}


 void saveEvents_3p3h(double Ein, std::string name, int nr, nucleus_str* nuc) {
 
 event_3p3h* e = new event_3p3h;
 const char * name_c = name.c_str();
 TFile *ff = new TFile (name_c, "recreate");
 TTree* tf = new TTree ("treeout", "Tree of events");
 tf->Branch ("e", "event_3p3h", &e);
 delete e;
 
 
 int i = 0;
	 int nr_events=nr*1000000;
 int counter = 1;
	 double max_weight = 0;
	 double total_xsec;
	 double rho, r, w, theta_lepton, qq, p_nuc, phi_nuc, theta_nuc, weight, k_fermi;
	 std::complex<double> A00, A11, A12, A03, A33;
	 bool result_zero=false;
	 if (Global) {
		 rho = (getRho(0., nuc, PROTON)+getRho(0., nuc, NEUTRON))/fermi3;
		 k_fermi = pow(3*M_PI*M_PI/2.*rho,1./3.);
	 }
	 
	 std::ifstream f;
	 std::ostringstream przedNazwa1;
	 przedNazwa1<<"3p3h_weight_"<<nuc->name<<"_"<<Ein<<".txt";
	 f.open(przedNazwa1.str());
	 if ( f.fail() ) {
		 std::cerr<<"No file encountered. Generating test events...\n";
		 test_events(Ein, nuc, max_weight,total_xsec);
		 std::ofstream f2;
		 f2.open(przedNazwa1.str());
		 f2<<max_weight<<" "<<total_xsec;
		 f2.close();
	 }
	 else {
		 f>>max_weight>>total_xsec;
		 f.close();
	 }
 
 while (i<nr_events) {
	 
	 w = frandom()*(Ein-lepton_mass);
	 theta_lepton = frandom()*M_PI;
	 
	 if (!Global) {
		 r = frandom()*nuc->rmax;
		 rho = (getRho(r, nuc, PROTON)+getRho(r, nuc, NEUTRON))/fermi3;
		 k_fermi = pow(3*M_PI*M_PI/2.*rho,1./3.);
	 }
	 p_nuc = frandom()*k_fermi;
	 theta_nuc = frandom()*M_PI;
	 phi_nuc=frandom()*2.*M_PI;
	 
	 double eout = Ein-w;
	 double kout = sqrt(eout*eout-lepton_mass*lepton_mass);
	 qq = sqrt(Ein*Ein+kout*kout-2*Ein*kout*cos(theta_lepton));
	 
	 if (qq<w) result_zero=true;
	 //if (qq>1.2*GeV) result_zero=true;
	 
	 init_deltaFF(-w*w+qq*qq,w);
	 
	 if (!result_zero) {
	 hadron_delta_p_phi_theta(theta_nuc, phi_nuc, p_nuc, w-nuc->Q_2p2h, qq, rho, true, false, false, A00, A11, A12, A03, A33);

	 
	 double sint2 = pow(sin(theta_lepton/2.),2);
	 double cost2 = pow(cos(theta_lepton/2.),2);
	 
	 double W1= A11.real()/2.;
	 double W2 =0.5*(A00.real()+A11.real()+ w*w/qq/qq*(A33.real()-A11.real())-2*w/qq*A03.real() );
	 double W3 =A12.imag()/qq;
	 double W4 =(A33.real()-A11.real())/2./qq/qq;
	 double W5 =(A03.real()-w/qq*(A33.real()-A11.real()))/qq;
	 
	 
	 double as = lepton_mass*lepton_mass/eout/(eout+kout);
	 
	 double sum1 = 2*W1*sint2 + W2*cost2 - W3*(Ein+eout)*sint2;
	 double sum2 = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) + W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
	 
	 
	 
	 double coef = kout*eout*GF*GF/M_PI/M_PI*sin(theta_lepton);
	 
	 double factor_rho = NR_PART/rho/2./M_PI;
	 if (!Global) factor_rho = 2*r*r*pow(fermi,3);
	 
		 if (!Global) weight = coef*(sum1+as*sum2)*factor_rho*2*M_PI*M_PI*k_fermi*M_PI*(Ein-lepton_mass)*nuc->rmax*2*M_PI;
		 else weight = coef*(sum1+as*sum2)*factor_rho*2*M_PI*M_PI*k_fermi*M_PI*(Ein-lepton_mass)*2*M_PI;
		 
	 
	 if (weight> max_weight) {
		 std::cerr<<"max weight changed! "<<weight<<" previous: "<<max_weight<<"\n";
		 max_weight = weight;
	 }
	 
	 if (weight> max_weight*frandom()) { //we accept this event
		 
		 e= new event_3p3h;
		 e->q0 = w;
		 e->theta_lepton = theta_lepton;
		 e->r=r;
		 
		 event_3p3h_gaussian(k_fermi, w, qq, e);
		 
		 tf->Fill();
		 
		 delete e;
	 }
	 }
	result_zero=false;
	 
 if (i*10./nr_events == counter) {std::cerr<<"* "; counter++;}
 i++;
 
 }
	 
 
 TVectorD constants(4);
 constants[0]=nr_events;
 constants[1]=Ein;
 constants[2]=lepton_mass;
 constants[3]=total_xsec;
 
 constants.Write("constants");
 
 ff->Write ();
 ff->Close ();
 delete tf;
 delete ff;
 }
 
 
 
 void event_3p3h_gaussian(double k_fermi, double w, double qq,  event_3p3h* e)
 {
	 
	 bool possible_kin = false;
	 double W;
	 while (!possible_kin) {
	 
	 vec k1 = rand_from_ball(k_fermi);
	 vec k2 = rand_from_ball(k_fermi);
	 vec k3 = rand_from_ball(k_fermi);
	 vec q(0,0,qq);
	 vect total(k1+k2+k3+q, w+sqrt(k1*k1+M*M)+sqrt(k2*k2+M*M)+sqrt(k3*k3+M*M) );
		 W = sqrt(total*total);
		 if (total*total > 9.*(M*M+k_fermi*k_fermi)) possible_kin = true;
	 }

	 vect out1, out2, out3;
	 bool pauli_blocked = true;
	 while (pauli_blocked) {
		 threebody(W, M,M,M,out1,out2,out3);
		 if (out1.length()>=k_fermi && out2.length()>=k_fermi && out3.length()>=k_fermi) pauli_blocked=false;
	 }
	 
	 e->out1 = out1;
	 e->out2 = out2;
	 e->out3 = out3;

	 
 }


void threebody (double W, double m1, double m2, double m3, vect &p1, vect &p2, vect &p3)
{
	if (W < m1+m2+m3)
	{
		cout<<"impossible kinematics in 3-body decay"<<endl;
	}
	
	double W_2 = W*W;
	double m3_2 = m3*m3;
	double m1_2 = m1*m1;
	double m2_2 = m2*m2;
	
	double m23_2_los, m23_2_min, m23_2_max, E3;
	
	do
	{
		double m12_2_los = (m1+m2)*(m1+m2) + frandom()* ( (W-m3)*(W-m3) - (m1+m2)*(m1+m2) );
		double m23_test = (m2+m3)*(m2+m3) + frandom()* ( (W-m1)*(W-m1) - (m2+m3)*(m2+m3) );
		m23_2_los = m23_test;
		
		E3 = (W_2 + m3_2 - m12_2_los)/2.0/W;
		
		double E2star = (m12_2_los - m1_2 + m2_2)/2.0/sqrt(m12_2_los);
		double E3star = (W_2 - m12_2_los - m3_2)/2.0/sqrt(m12_2_los);
		
		m23_2_min = (E2star+E3star)*(E2star+E3star)
		- ( sqrt(E2star*E2star - m2_2) + sqrt (E3star*E3star - m3_2) )*
		( sqrt(E2star*E2star - m2_2) + sqrt (E3star*E3star - m3_2) );
		
		m23_2_max = (E2star+E3star)*(E2star+E3star)
		- ( sqrt(E2star*E2star - m2_2) - sqrt (E3star*E3star - m3_2) )*
		( sqrt(E2star*E2star - m2_2) - sqrt (E3star*E3star - m3_2) );
		
	}
	while ( (m23_2_los < m23_2_min) || (m23_2_los > m23_2_max) );
	
	
	double E1 = (W_2 + m1_2 - m23_2_los)/2.0/W;
	double E2 = W - E1 - E3;
	
	double mom1 = sqrt (E1*E1 - m1_2);
	double mom2 = sqrt (E2*E2 - m2_2);
	double mom3 = sqrt (E3*E3 - m3_2);
	
	vec axis = rand_dir ();		 // new -- Y axis
	vec proba = rand_dir ();
	
	vec kier1 = vecprod (axis, proba);
	kier1 = kier1/kier1.length();//unit vector -- new Z axis!
	
	// unit vetor -- new X axis
	vec kier2 = vecprod (axis, kier1);
	
	double kosalpha = (mom1*mom1 + mom2*mom2 - mom3*mom3)/2.0/mom1/mom2;
	double sinalpha = sqrt(1-kosalpha*kosalpha);
	double kosbeta = ( mom1 - mom2*kosalpha )/mom3;
	double sinbeta = sqrt(1 - kosbeta*kosbeta);
	
	vec pp1 = kier1*mom1;
	vec pp2 = mom2*(-kosalpha*kier1 + sinalpha*kier2);
	vec pp3 = mom3*(-kosbeta*kier1 - sinbeta*kier2);
	
	p1 = vect(pp1,E1);
	p2 = vect(pp2,E2);
	p3 = vect(pp3,E3);
	
}


void hadron_delta_p_phi_theta(double theta, double phi, double mom_nuc, double w, double qq, double rho, bool ca3, bool ca2, bool cq, std::complex<double>& A00, std::complex<double>& A11, std::complex<double>& A12, std::complex<double>& A03, std::complex<double>& A33){

	
	double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);
	double costheta_nuc = cos(theta);
	
	D4V<double> q(w,0,0,qq);
	
	D4V<double> p(sqrt(M*M+mom_nuc*mom_nuc),mom_nuc*sqrt(1.-costheta_nuc*costheta_nuc)*cos(phi),mom_nuc*sqrt(1.-costheta_nuc*costheta_nuc)*sin(phi),mom_nuc*costheta_nuc); //nukleon
	D4V<double> p_sl(sqrt(M*M+mom_nuc*mom_nuc)-k_fermi*k_fermi/2./M*FERMI,mom_nuc*sqrt(1.-costheta_nuc*costheta_nuc)*cos(phi),mom_nuc*sqrt(1.-costheta_nuc*costheta_nuc)*sin(phi),mom_nuc*costheta_nuc); //nukleon
	
	//	std::cerr<<"M q2 pq "<<sqrt(p*p)<<" "<<q*q<<" "<<p*q<<"\n";
	
	DM pslash=p*dirac;
	DM qslash=q*dirac;
	
	//Delta in the medium
	double p_delta2 = (p_sl+q)*(p_sl+q);
	if (p_delta2<0) {
		p_delta2=0;
	}
	double denom =im_delta_denominator(p_delta2, rho, ca3, ca2, cq);
	
	D4T<DM> projD, del_corr;
	projection_delta(p+q, projD); //P_{\alfa\beta}
	

    delta_32(p+q, projD);
    projD*=denom;

	
	D4T<DM> verD_a, verD_v;
	vertex_delta(p, q, verD_a, verD_v); //Gamma_V^{\beta\nu}, Gamma_A^{\beta\nu}
	verD_v+=verD_a;//Gamma^{\beta\nu} = Gamma_V^{\beta\nu}+Gamma_A^{\beta\nu}
	
	projD*=verD_v; //P_{\alfa\beta} Gamma^{\beta\nu}
	
	verD_v.hermit(); // (Gamma^{\beta\nu})^dagger
	
	verD_v*=gamma0; // (Gamma^{\beta\nu})^dagger \gamma0
	verD_v&=gamma0;  // \gamma (Gamma^{\beta\nu})^dagger \gamma0
	
	projD *=(pslash+M_N*identi);  //P_{\alfa\beta} Gamma^{\beta\nu}  (\pslash+M)
	
	verD_v*=projD;
	
	D4T<std::complex<double>> Amunu(verD_v.CT(verD_v));
	
	//	std::cerr<<"Amunu: "<<Amunu(0,0)<<" "<<Amunu(1,1)<<"\n";
	if (denom!=denom) std::cerr<<"p_d^2="<<(p_sl+q)*(p_sl+q)<<"\n";
	
	
	double C_n = 1;
	double C_p = sqrt(3.);
	
	double coeff=-2.*mom_nuc*mom_nuc*cabbibo2*(C_n*C_n+C_p*C_p)/pow(2.*M_PI,3)/p(0)
	*sqrt(1-costheta_nuc*costheta_nuc)/2.; // /2 because of the trace
	

	
	A00 = Amunu(0,0)*coeff;
	A11 = Amunu(1,1)*coeff;
	A12 = Amunu(2,1)*coeff;
	A33 = Amunu(3,3)*coeff;
	A03 = Amunu(0,3)*coeff;
	
}


double total_xsec_3p3h(double Ein, int nr, nucleus_str* nuc) {



int i = 0;
int nr_events=nr*100000;
int counter = 1;
    double total_xsec=0;
    double rho, r, w, theta_lepton, qq, p_nuc, phi_nuc, theta_nuc, weight, k_fermi;
    std::complex<double> A00, A11, A12, A03, A33;
    bool result_zero=false;
    if (Global) {
        rho = (getRho(0., nuc, PROTON)+getRho(0., nuc, NEUTRON))/fermi3;
        k_fermi = pow(3*M_PI*M_PI/2.*rho,1./3.);
    }
    

while (i<nr_events) {
    
    
    w = frandom()*(Ein-lepton_mass);
    if (w-nuc->Q_2p2h<0) result_zero=true;
    
    theta_lepton = frandom()*M_PI;
    
    if (!Global) {
        r = frandom()*nuc->rmax;
        rho = (getRho(r, nuc, PROTON)+getRho(r, nuc, NEUTRON))/fermi3;
        k_fermi = pow(3*M_PI*M_PI/2.*rho,1./3.);
    }
    p_nuc = frandom()*k_fermi;
    theta_nuc = frandom()*M_PI;
    phi_nuc=frandom()*2.*M_PI;
    
    double eout = Ein-w;
    double kout = sqrt(eout*eout-lepton_mass*lepton_mass);
    qq = sqrt(Ein*Ein+kout*kout-2*Ein*kout*cos(theta_lepton));
    
    if (qq<w) result_zero=true;
   // if (qq>1.2*GeV) result_zero=true;
    
    init_deltaFF(-w*w+qq*qq,w);
    
    if (!result_zero) {
    hadron_delta_p_phi_theta(theta_nuc, phi_nuc, p_nuc, w-nuc->Q_2p2h, qq, rho, true, false, false, A00, A11, A12, A03, A33);

    
    double sint2 = pow(sin(theta_lepton/2.),2);
    double cost2 = pow(cos(theta_lepton/2.),2);
    
    double W1= A11.real()/2.;
    double W2 =0.5*(A00.real()+A11.real()+ w*w/qq/qq*(A33.real()-A11.real())-2*w/qq*A03.real() );
    double W3 =A12.imag()/qq;
    double W4 =(A33.real()-A11.real())/2./qq/qq;
    double W5 =(A03.real()-w/qq*(A33.real()-A11.real()))/qq;
    
    
    double as = lepton_mass*lepton_mass/eout/(eout+kout);
    
    double sum1 = 2*W1*sint2 + W2*cost2 - W3*(Ein+eout)*sint2;
    double sum2 = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) + W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
    
    
    double sum1_anti = 2*W1*sint2 + W2*cost2 + W3*(Ein+eout)*sint2;
    double sum2_anti = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) - W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
            
    double sum_ =sum1+as*sum2;
            
    if (CH_2p2h==CCanti) sum_ =sum1_anti+as*sum2_anti;
   
    double coef = kout*eout*GF*GF/M_PI/M_PI*sin(theta_lepton);
    
    double factor_rho = NR_PART/rho/2./M_PI;
    if (!Global) factor_rho = 2*r*r*pow(fermi,3);//4*M_PI*r*r/2./M_PI;
    
    if (!Global) weight = coef*sum_*factor_rho*2*M_PI*M_PI*k_fermi*M_PI*(Ein-lepton_mass)*nuc->rmax*2*M_PI;
    else weight = coef*sum_*factor_rho*2*M_PI*M_PI*k_fermi*M_PI*(Ein-lepton_mass)*2*M_PI;
        
    total_xsec += weight;
        
    }
   result_zero=false;
    
if (i*10./nr_events == counter) {std::cerr<<"* "; counter++;}
i++;

}
     total_xsec=total_xsec/nr_events*pow(10,38)/cm2;
    return total_xsec;
}
