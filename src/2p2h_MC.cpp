
#include "2p2h_MC.h"


/**********************
 
 
 SAVE EVENTS - MONTE CARLO
 
 **********************/

void saveEvents_2p2h(double Ein, std::string name, int nr, nucleus_str* nuc) {
	
	event* e = new event;
	const char * name_c = name.c_str();
	TFile *ff = new TFile (name_c, "recreate");
	TTree* tf = new TTree ("treeout", "Tree of events");
	tf->Branch ("e", "event", &e);
	delete e;
	int aaaa=0;
	/**
	 * In the case of global Fermi gas, you have the same density for all the nucleus and hence the same Fermi momentum
	 * The two lines below only apply in the GFG case. Otherwise, we regenerate them every event in the LFG case
	 */
	double rho = getRho(0, nuc,true)/fermi3+getRho(0, nuc,false)/fermi3; //for GFG this will be constant
	double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);
	
	std::complex<double> cunuc(0.,0.);  // NUCLEAR?
	std::complex<double> udel(0.,0.);  // DELTA?
	
	int i = 0;
	int nr_events=nr*10000; //number of produced events
	int counter = 1;
	double Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, theta_lepton,Eout_lepton ,Ep2,r=0;
	bool result_zero = false;
	bool result_zero_lind = false; //if the Lindhard function is 0 ?
	
	/**
	 * We impose a cut for the initial neutrino energy. You need to have enough energy to produce the lepton mass + energy to remove the nucleus from the fermi sea
	 */
	if (Ein<lepton_mass+nuc->Q_2p2h) {
		std::cerr<<"Initial neutrino energy too small! "<<Ein<<"<"<<lepton_mass+nuc->Q_2p2h<<"\n";
		return;
	}
	
	double total_xsec=0;
    
	//here we produce events
	while (i<nr_events) {
		
		e= new event; // The event we will save the info in a tree
		
		if (!Global) {
            /**
             * This is the Radial dependency of LFG. This is our 1st variable to integrate over
             */
			r = nuc->rmax*frandom(); //generate a point in nucleus with radius r
			rho = getRho(r, nuc,true)/fermi3+getRho(r, nuc,false)/fermi3; //get the density at this point as the sum of the proton and neutron densities
			k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.); //get the Fermi momentum at this point
		}
		
        double corr_in_medium = nuc->Q_2p2h; //energy needed to remove the nucleus from the Fermi sea
        /**
         * We generate the energy of the outgoing lepton. The energy will exist between two limits: lepton mass and the energy transfered to the nucleon by the W boson
         * This is our 2nd variable to integrate over
         */
        Eout_lepton = lepton_mass+(Ein-lepton_mass-corr_in_medium)*frandom(); // corrected energy of the outgoing lepton

        w = Ein-Eout_lepton-corr_in_medium;
        /**
         * We generate the angle of the outgoing lepton. This is our 3rd variable to integrate over
         */
		theta_lepton = M_PI*frandom();  //generate the outgoing lepton angle
		double kout = sqrt(Eout_lepton*Eout_lepton-lepton_mass*lepton_mass);
		qq = sqrt(Ein*Ein + pow(kout,2) - 2*kout*Ein*cos(theta_lepton));
		
		if (qq<w) result_zero=true;
		if (qq>1.2*GeV) result_zero=true; //we impose a cut for the transfered momentum
		/**
		 * We generate the phi angle of the 1st nucleon. This is our 4th variable to integrate over. This Phi is with respect to the q-pi vector
		 */
		phi_p = 2*M_PI*frandom(); //generate the phi angle of the 1st initial nucleons
		/**
		 * Generate the Four momentum of the pion. This is our 5th, 6th, 7th and 8th variables to integrate over
		 */
		kk = (qq+3*k_fermi)*frandom();
		phi_pi = 2*M_PI*frandom();
		thetapi=M_PI*frandom();
		Epi = w*frandom();
		
		/**
		 * This is the Lindhard function. This calculates the screening of the pion in the nucleus?
		 */
		ulind(Epi/mpi, kk/mpi, k_fermi/mpi,cunuc, udel);
		cunuc*=pow(mpi,2);
		
		if (std::abs(cunuc.imag())<pow(10,-20)) result_zero=true; //if the pion cannot produce 1p1h (second bubble)
		
		if (!(qq<w) && (qq>1.2*GeV) && !(std::abs(cunuc.imag())<pow(10,-20))) aaaa++;

        /**
         * We have to now calculate the kinematic variables to apply some cuts on what events are processed. We only calculate valid kimematic events to save time
         */
		double cosi=cos(thetapi);
		double kPi1=kk*sqrt(1-cosi*cosi)*cos(phi_pi);
		double kPi2 = kk*sqrt(1-cosi*cosi)*sin(phi_pi);
		double kPi3 = kk*cosi;
		
		double lindhard = DUR(w-Epi, sqrt(pow(kPi1,2)+pow(kPi2,2)+pow(kPi3-qq,2)),k_fermi, k_fermi);
		double q_minus_k=sqrt(pow(kPi1,2)+pow(kPi2,2)+pow(kPi3-qq,2));
		if (std::abs(lindhard)<pow(10,-20)) result_zero_lind = true;
		
		double Emax=sqrt(M*M+k_fermi*k_fermi);
		double pom = std::max(M,Emax-(w-Epi));
        double Emin=std::max(pom, 0.5*(-(w-Epi)+q_minus_k*sqrt(1.-4.*M*M/((w-Epi)*(w-Epi)-q_minus_k*q_minus_k))));
		
        if (Emin>Emax) result_zero= true;
        
        if (Emin>Emax and result_zero_lind==false) std::cerr<<"Wrong energy "<<Emax<<" "<<Emin<<"\n";
        double Emax_2 = sqrt(M*M+k_fermi*k_fermi);
        double pom_2 = std::max(M,Emax_2-Epi);
        double Emin_2 = std::max(pom_2, 0.5*(-(Epi)+kk*sqrt(1.-4.*M*M/(Epi*Epi-kk*kk))));

        /**
         * We generate the energy of the 2nd nucleon. This is our 9th variable to integrate over
         */
        
        Ep2 =  Emin_2+(Emax_2-Emin_2)*frandom(); //generate energy of the 2nd nucleon
        if (Emin_2>Emax_2) result_zero = true;

        /**
         * We generate the energy of the 1st nucleon. This is our 10th variable to integrate over
         */
		Ep =Emin+(Emax-Emin)*frandom(); //generate energy of the 1st nucleon
		

        if (Emin>Emax) result_zero= true;
        if (Emin_2>Emax_2) result_zero = true;

        /**
         * This is the ROOT tree in which we save all the variables. Some are added here while the rest are done inside the event_2p2h_rho_gaussian function
         */
		
		e->q0 = Ein-Eout_lepton;
        e->qz = qq;
		e->theta_lepton = theta_lepton;
		e->rho = rho;
		e->r=r;
		e->E_pi = Epi;
		e->theta_pi = thetapi;
		e->phi_pi=phi_pi;
		e->p_pi = kk;
		e->E_nuc =Ep;
		e->phi_nuc = phi_p;
		e->E_nuc2 = Ep2;
        e->fermi_momentum = k_fermi;
        e->lindhard_Value = lindhard;
        e->isCC = (CH_2p2h == 0);
        e->energy_neutrino = Ein;

        /**
         * This is the main function that calculates the cross section of the 2p2h event. It will also calculate the weights of the different channels. This contains 3 diagram classes
         * Some info is added to the tree that has additonal info on what the event is internally, What diagram and so on. this is needed to be able to understand the topology of the event
         */

		
		if (!result_zero && !result_zero_lind) event_2p2h_rho_gaussian(r, Ein, theta_lepton, Ep2, Ep, phi_p, kk, phi_pi, thetapi, Epi, w, qq, rho, 1, Emax_2-Emin_2, e, nuc,corr_in_medium);
		else {
			e->weight_pp=0;
			e->weight_np=0;
			e->weight_pn=0;
			e->weight_anti_nn=0;
			e->weight_anti_pn=0;
			e->weight_anti_np=0;
		}

        /**
         * Calculate the integral volume of this event. This will contain the Jacobian of the transformation of variables and the range of integration
         */
		double coeff = w*sin(theta_lepton)*(Emax-Emin)*(qq+3*k_fermi)*(Ein-lepton_mass-corr_in_medium)*Pi*2*Pi*2*Pi*Pi;
//                | E_pi|     theta_mu    |    EN1    |  P_pi        |           E_mu                 |theta_pi| phiN1 | Phi_pi | theta_mu(its pi*sin(theta_lepton))


		if (!Global) coeff*=nuc->rmax; //for LFG we integrate over r

        e->int_vol_E_pi = w;
        e->int_vol_theta_lepton = sin(theta_lepton)*Pi;
        e->int_vol_E_nuc = (Emax-Emin);
        e->int_vol_p_pi = (qq+3*k_fermi);
        e->int_vol_phi_nuc = 2*Pi;
        e->int_vol_E_nuc2 = Emax_2-Emin_2;
        e->int_vol_E_lepton = Ein-lepton_mass-corr_in_medium;
        e->int_vol_theta_pi = Pi;
        e->int_vol_phi_pi = 2*Pi;
        e->int_vol_r = nuc->rmax;
		
		
		e->weight_pp *=coeff ;//range of integration
		e->weight_np *=coeff;
		e->weight_pn*=coeff;
		

        
        if ((e->weight_pp+e->weight_np+e->weight_pn)!=0) {
            tf->Fill(); //we do not include 0
            i++;

        }
		
		delete e;
		result_zero=false;
		result_zero_lind=false;
		
		if (i*10./nr_events == counter) {std::cerr<<"* "; counter++;}

		
		total_xsec+=(e->weight_pp+e->weight_np+e->weight_pn);

		
	}

    total_xsec*=1./nr_events*pow(10,38)/cm2*2*Pi;
	
	std::cerr<<"TOTAL X-sec: "<<total_xsec<<"\n";
	
	TVectorD constants(5);
	constants[0]=nr_events;
	constants[1]=Ein;
	constants[2]=lepton_mass;
	constants[4]=nuc->Q_2p2h;
	constants[3]=pow(10,38)/cm2*GeV;
		
	constants.Write("constants");
	
	ff->Write ();
	ff->Close ();
	delete ff;
}


void event_2p2h_rho_gaussian(double r, double Ein, double theta_lepton, double Ep2, double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double w, double qq, double rho, double cos_theta_nuc, double E_range, event* e, nucleus_str* nuc, double corr_in_medium)
{
	
	init_deltaFF(-w*w+qq*qq,w);
	
	D4V<double> q(w, 0,0,qq); //boson W
	D4V<DM> axial_nuc, vector_nuc, scalar_nuc;
	Wnn(q, vector_nuc,axial_nuc,scalar_nuc);

    /**
     * This is the hadron tensor being defined here.
     * The A00, A11, A12, A03, A33 are the 5 components of the first diagram class
     * The B00, B11, B12, B03, B33 are the 5 components of the second diagram class
     * The C00, C11, C12, C03, C33 are the 5 components of the third diagram class
     */
	std::complex<double> A00[3], B00[3], C00[3];
	std::complex<double> A11[3], B11[3], C11[3];
	std::complex<double> A12[3], B12[3], C12[3];
	std::complex<double> A03[3], B03[3], C03[3];
	std::complex<double> A33[3], B33[3], C33[3];
	for (int i=0; i<3; i++) {
		A00[i]=A11[i]=A12[i]=A03[i]=A33[i]=std::complex<double>(0,0);
		B00[i]=B11[i]=B12[i]=B03[i]=B33[i]=std::complex<double>(0,0);
		C00[i]=C11[i]=C12[i]=C03[i]=C33[i]=std::complex<double>(0,0);
		
	}
	//w-=nuc->Q;
	/**
	 * Rho Transversal component
	 */
    hadron_2p2h_rho_Epi_thetapi_phipi_kpi_phip_Ep_isospin(Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, A00, A11, A12, A03, A33);
    /**
     * Pion Longitudinal component
     */
    hadron_2p2h_part1_Epi_thetapi_phipi_kpi_phip_Ep_isospin(Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, B00, B11, B12, B03, B33);
    /**
     * This is the Interference term. Its small but has the influence on negative weights. This needs to be studied in particular
     */
    hadron_2p2h_part2_Epi_thetapi_phipi_kpi_E1_E2_isospin(Ep2, Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, C00, C11, C12, C03, C33);
    for (int i=0; i<3; i++) {
        C00[i]*=E_range/(2*M_PI); //because in this case we don't integrate over phi_nuc but instead over Ep_2
        C11[i]*=E_range/(2*M_PI);
        C12[i]*=E_range/(2*M_PI);
        C03[i]*=E_range/(2*M_PI);
        C33[i]*=E_range/(2*M_PI);
    }
	
	if (A00[0]!=A00[0] || A11[0]!=A11[0] || A12[0]!=A12[0] || A03[0]!=A03[0] || A33[0]!=A33[0]) {
		std::cerr<<"ERROR "<<Ep<<" phi="<<phi_p<<" kpi="<<kk<<" phipi="<<phi_pi<<" thpi="<<thetapi<<" Epi="<<Epi<<" q="<<qq<<" w="<<w<<" rho="<<rho<<" "<<A00[0]<<" "<<A11[0]<<" "<<A12[0]<<" "<<A03[0]<<" "<<A33[0]<<"\n";
	}

//    for(int ite=0; ite<3; ite++){
//        e->rho_transversal_00 [ite] = A00[ite];
//        e->rho_transversal_11 [ite] = A11[ite];
//        e->rho_transversal_12 [ite] = A12[ite];
//        e->rho_transversal_03 [ite] = A03[ite];
//        e->rho_transversal_33 [ite] = A33[ite];
//
//        e->pi_longitudinal_00 [ite] = B00[ite];
//        e->pi_longitudinal_11 [ite] = B11[ite];
//        e->pi_longitudinal_12 [ite] = B12[ite];
//        e->pi_longitudinal_03 [ite] = B03[ite];
//        e->pi_longitudinal_33 [ite] = B33[ite];
//
//        e->both_interference_00 [ite] = C00[ite];
//        e->both_interference_11 [ite] = C11[ite];
//        e->both_interference_12 [ite] = C12[ite];
//        e->both_interference_03 [ite] = C03[ite];
//        e->both_interference_33 [ite] = C33[ite];
//    }



    /**
     * We add the 3 diagram classes together
     */
    
    
	for (int i=0; i<3; i++) {
		A00[i]+=B00[i]+C00[i];
		A11[i]+=B11[i]+C11[i];
		A12[i]+=B12[i]+C12[i];
		A03[i]+=B03[i]+C03[i];
		A33[i]+=B33[i]+C33[i];
	}
	
	w+=corr_in_medium;
	
	double eout = Ein-w;
	double kout = sqrt(eout*eout-lepton_mass*lepton_mass);
	
	double sint2 = pow(sin(theta_lepton/2.),2);
	double cost2 = pow(cos(theta_lepton/2.),2);
	

	
	for (int i=0; i<3; i++) {

        /**
         * These are the Coefficients of the hadron tensor.
         * W1 -> Transverse Response
         * W2 -> Longitudinal Response + Transverse Response
         * W3 -> Parity violating Response
         * W4 -> Interference term (higher order contribution)
         * W5 -> Axial Vector Response
         */
		
		double W1= A11[i].real()/2.;
		double W2 =0.5*(A00[i].real()+A11[i].real()+ w*w/qq/qq*(A33[i].real()-A11[i].real())-2*w/qq*A03[i].real() );
		double W3 =A12[i].imag()/qq;
		double W4 =(A33[i].real()-A11[i].real())/2./qq/qq;
		double W5 =(A03[i].real()-w/qq*(A33[i].real()-A11[i].real()))/qq;
		
		e->transverse_response = W1;
        e->mixed_response = W2;
        e->parity_violating_response = W3;
        e->interference_response = W4;
        e->axial_vector_response = W5;
		
		double as = lepton_mass*lepton_mass/eout/(eout+kout);
		
		double sum1 = 2*W1*sint2 + W2*cost2 - W3*(Ein+eout)*sint2;
		double sum2 = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) + W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
		
		double sum1_anti = 2*W1*sint2 + W2*cost2 + W3*(Ein+eout)*sint2;
		double sum2_anti = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) - W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
		
		double coef = kout*eout*GF*GF/M_PI/M_PI;
		
		double factor_rho = NR_PART/rho/2./M_PI;
		if (!Global) factor_rho = 2*r*r*pow(fermi,3);
		
		
        if (CH_2p2h == CCanti) {
            
            switch(i){
                case 0:
                    e->weight_pp = 0;
                    e->weight_anti_nn = coef*(sum1_anti+as*sum2_anti)*factor_rho;
                    break;
                case 1:
                    e->weight_np = 0;
                    e->weight_anti_pn = coef*(sum1_anti+as*sum2_anti)*factor_rho;
                    break;
                case 2:
                    e->weight_pn = 0;
                    e->weight_anti_np = coef*(sum1_anti+as*sum2_anti)*factor_rho;
                    break;
            }
            
        }
        else
		switch(i){
			case 0:
				e->weight_pp = coef*(sum1+as*sum2)*factor_rho;
				e->weight_anti_nn = 0;
				break;
			case 1:
				e->weight_np = coef*(sum1+as*sum2)*factor_rho;
				e->weight_anti_pn = 0;
				break;
			case 2:
				e->weight_pn = coef*(sum1+as*sum2)*factor_rho;
				e->weight_anti_np = 0;
				break;
		}
		
	}
	
	
}


/**********************
 
 
 TOTAL CROSS SECTION - MONTE CARLO INTEGRATION
 
 **********************/



void total_xsec_2p2h(double Ein, int nr, nucleus_str* nuc, double& total_xsec) {
	
	
	total_xsec = 0;

	double weight=0;
	double rho = getRho(0, nuc,true)/fermi3+getRho(0, nuc,false)/fermi3;
	double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);
	
	std::complex<double> cunuc(0.,0.);
	std::complex<double> udel(0.,0.);
	
	int i = 0;
	int nr_events=nr*1000000;
	int counter = 1;
	double Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, theta_lepton,Eout_lepton ,Ep2,r=0;
	bool result_zero = false;
	
	if (Ein<lepton_mass+nuc->Q_2p2h) {
		std::cerr<<"Initial neutrino energy too small! "<<Ein<<"<"<<lepton_mass+nuc->Q_2p2h<<"\n";
		return;
	}
	
	while (i<nr_events) {
		
		if (!Global) {
			r = nuc->rmax*frandom();
			rho = getRho(r, nuc,true)/fermi3+getRho(r, nuc,false)/fermi3;
			k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);
		}
		
        double corr_in_medium = nuc->Q_2p2h;
        
        Eout_lepton = lepton_mass+corr_in_medium+(Ein-lepton_mass-corr_in_medium)*frandom();
        w = Ein-Eout_lepton-corr_in_medium;
        
        
		theta_lepton = M_PI*frandom();
		double kout = sqrt(Eout_lepton*Eout_lepton-lepton_mass*lepton_mass);
		qq = sqrt(Ein*Ein + pow(kout,2) - 2*kout*Ein*cos(theta_lepton));
		
		if (qq<w) result_zero=true;
		if (qq>1.2*GeV) result_zero=true;
		
		phi_p = 2*M_PI*frandom();
		kk = (qq+3*k_fermi)*frandom();
		phi_pi = 2*M_PI*frandom();
		thetapi=M_PI*frandom();
		Epi = w*frandom();
		
		
        ulind(Epi/mpi, kk/mpi, k_fermi/mpi,cunuc, udel);
		cunuc*=pow(mpi,2);
		
		if (std::abs(cunuc.imag())<pow(10,-20)) result_zero=true;
		
		double Emax_2 = sqrt(M*M+k_fermi*k_fermi);
		double pom_2 = std::max(M,Emax_2-Epi);
		double Emin_2 = std::max(pom_2, 0.5*(-(Epi)+kk*sqrt(1.-4.*M*M/(Epi*Epi-kk*kk))));
		
		Ep2 =  Emin_2+(Emax_2-Emin_2)*frandom();
		
		double cosi=cos(thetapi);
		double kPi1=kk*sqrt(1-cosi*cosi)*cos(phi_pi);
		double kPi2 = kk*sqrt(1-cosi*cosi)*sin(phi_pi);
		double kPi3 = kk*cosi;
		
		double lindhard = DUR(w-Epi, sqrt(pow(kPi1,2)+pow(kPi2,2)+pow(kPi3-qq,2)),k_fermi, k_fermi);
		double q_minus_k=sqrt(pow(kPi1,2)+pow(kPi2,2)+pow(kPi3-qq,2));
		if (std::abs(lindhard)<pow(10,-20)) result_zero = true;
		
		double Emax=sqrt(M*M+k_fermi*k_fermi);
		double pom = std::max(M,Emax-(w-Epi));
		double Emin=std::max(pom, 0.5*(-(w-Epi)+q_minus_k*sqrt(1.-4.*M*M/((w-Epi)*(w-Epi)-q_minus_k*q_minus_k))));
		
		Ep =Emin+(Emax-Emin)*frandom();
        
        if (Emin>Emax) result_zero= true;
        if (Emin_2>Emax_2) result_zero = true;
        
		weight=0;
        
		if (!result_zero)  MCpoint_2p2h_gaussian(r, Ein, theta_lepton, Ep2, Ep, phi_p, kk, phi_pi, thetapi, Epi, w, qq, rho, 1, Emax_2-Emin_2, weight,nuc,corr_in_medium);

				
		double coeff = w*sin(theta_lepton)*(Emax-Emin)*(qq+3*k_fermi)*(Ein-lepton_mass-corr_in_medium)*Pi*2*Pi*2*Pi*Pi;
		
		if (!Global) coeff*=nuc->rmax;
		
		
		weight*=coeff;
		
		total_xsec+=weight;
        

		
		result_zero=false;
		if (i*10./nr_events == counter) {std::cerr<<"* "; counter++;}
		i++;
		
	}
	
	total_xsec*=1./nr_events*pow(10,38)/cm2*2*Pi;


	
	
}

    


void MCpoint_2p2h_gaussian(double r, double Ein, double theta_lepton, double Ep2, double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double w, double qq, double rho, double cos_theta_nuc, double E_range, double& weight, nucleus_str* nuc, double corr_in_medium)
{
	
	init_deltaFF(-w*w+qq*qq,w);
	
	D4V<double> q(w, 0,0,qq); //boson W
	D4V<DM> axial_nuc, vector_nuc, scalar_nuc;
	Wnn(q, vector_nuc,axial_nuc,scalar_nuc);
	
	std::complex<double> A00[3], B00[3], C00[3];
	std::complex<double> A11[3], B11[3], C11[3];
	std::complex<double> A12[3], B12[3], C12[3];
	std::complex<double> A03[3], B03[3], C03[3];
	std::complex<double> A33[3], B33[3], C33[3];
	for (int i=0; i<3; i++) {
		A00[i]=A11[i]=A12[i]=A03[i]=A33[i]=std::complex<double>(0,0);
		B00[i]=B11[i]=B12[i]=B03[i]=B33[i]=std::complex<double>(0,0);
		C00[i]=C11[i]=C12[i]=C03[i]=C33[i]=std::complex<double>(0,0);
		
	}
	
    hadron_2p2h_rho_Epi_thetapi_phipi_kpi_phip_Ep_isospin(Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, A00, A11, A12, A03, A33);
    hadron_2p2h_part1_Epi_thetapi_phipi_kpi_phip_Ep_isospin(Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, B00, B11, B12, B03, B33);
    hadron_2p2h_part2_Epi_thetapi_phipi_kpi_E1_E2_isospin(Ep2, Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, rho, axial_nuc, vector_nuc, scalar_nuc, C00, C11, C12, C03, C33);
    for (int i=0; i<3; i++) {
        C00[i]*=E_range/(2*M_PI); //because in this case we don't integrate over phi_nuc but instead over Ep_2
        C11[i]*=E_range/(2*M_PI);
        C12[i]*=E_range/(2*M_PI);
        C03[i]*=E_range/(2*M_PI);
        C33[i]*=E_range/(2*M_PI);
    }
	
	if (A00[0]!=A00[0] || A11[0]!=A11[0] || A12[0]!=A12[0] || A03[0]!=A03[0] || A33[0]!=A33[0]) {
		std::cerr<<"ERROR "<<Ep<<" phi="<<phi_p<<" kpi="<<kk<<" phipi="<<phi_pi<<" thpi="<<thetapi<<" Epi="<<Epi<<" q="<<qq<<" w="<<w<<" rho="<<rho<<" "<<A00[0]<<" "<<A11[0]<<" "<<A12[0]<<" "<<A03[0]<<" "<<A33[0]<<"\n";
	}
	
	for (int i=0; i<3; i++) {
		A00[i]+=B00[i]+C00[i];
		A11[i]+=B11[i]+C11[i];
		A12[i]+=B12[i]+C12[i];
		A03[i]+=B03[i]+C03[i];
		A33[i]+=B33[i]+C33[i];
	}
	
	w+=corr_in_medium;
	
	double eout = Ein-w;
	double kout = sqrt(eout*eout-lepton_mass*lepton_mass);
	
	double sint2 = pow(sin(theta_lepton/2.),2);
	double cost2 = pow(cos(theta_lepton/2.),2);
	
	for (int i=0; i<3; i++) {
		
		double W1= A11[i].real()/2.;
		double W2 =0.5*(A00[i].real()+A11[i].real()+ w*w/qq/qq*(A33[i].real()-A11[i].real())-2*w/qq*A03[i].real() );
		double W3 =A12[i].imag()/qq;
		double W4 =(A33[i].real()-A11[i].real())/2./qq/qq;
		double W5 =(A03[i].real()-w/qq*(A33[i].real()-A11[i].real()))/qq;
		
		
		
		double as = lepton_mass*lepton_mass/eout/(eout+kout);
		
		double sum1 = 2*W1*sint2 + W2*cost2 - W3*(Ein+eout)*sint2;
		double sum2 = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) + W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
		
		double sum1_anti = 2*W1*sint2 + W2*cost2 + W3*(Ein+eout)*sint2;
		double sum2_anti = W1*cos(theta_lepton)-W2/2.*cos(theta_lepton) - W3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + W4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - W5*(eout+kout)/2.;
		
		double coef = kout*eout*GF*GF/M_PI/M_PI;
		
		double factor_rho = NR_PART/rho/2./M_PI;
		if (!Global) factor_rho = 2*r*r*pow(fermi,3);
        
        if (CH_2p2h!=CCanti)
            weight += coef*(sum1+as*sum2)*factor_rho;
        
        if (CH_2p2h==CCanti)
            weight += coef*(sum1_anti+as*sum2_anti)*factor_rho;
		
	}
	
	
}

void saveEvents_2p2h_variedEnergy(double Enu_low, double Enu_high, std::string name, int nr, nucleus_str* nuc){

    event* e = new event;
    const char * name_c = name.c_str();
    TFile *ff = new TFile (name_c, "recreate");
    TTree* tf = new TTree ("treeout", "Tree of events");
    tf->Branch ("e", "event", &e);
    delete e;
    int aaaa=0;

    double rho = getRho(0, nuc,true)/fermi3+getRho(0, nuc,false)/fermi3; //for GFG this will be constant
    double k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);

    std::complex<double> cunuc(0.,0.);
    std::complex<double> udel(0.,0.);

    int i = 0;
    int nr_events=nr*100000; //number of produced events
    int counter = 1;
    double Ep, phi_p, kk, phi_pi, thetapi, Epi, qq, w, theta_lepton,Eout_lepton ,Ep2,r=0;
    bool result_zero = false;
    bool result_zero_lind = false;
    double Ein =0;

    if (Enu_low<lepton_mass+nuc->Q_2p2h) {
        std::cerr<<"Initial neutrino energy too small! "<<Enu_low<<"<"<<lepton_mass+nuc->Q_2p2h<<"\n";
        Enu_low = lepton_mass+nuc->Q_2p2h;
    }

    double total_xsec=0;

    //here we produce events
    while (i<nr_events) {
        if(Enu_low == Enu_high){ // Here we have a fixed energy
            Ein = Enu_low;
        }
        else { // Here generate random Ein between E_low and E_high
            Ein = Enu_high * frandom();   // Generate a random energy for neutrino
            if (Ein < Enu_low) continue; // Only select energies above the minimum value
        }


//        if (Ein<lepton_mass+nuc->Q_2p2h) {
//            std::cerr<<"Initial neutrino energy too small! "<<Ein<<"<"<<lepton_mass+nuc->Q_2p2h<<"\n";
//            Enu_low = lepton_mass+nuc->Q_2p2h;
//        }

        e= new event;

        if (!Global) {
            r = nuc->rmax*frandom(); //generate a point in nucleus
            rho = getRho(r, nuc,true)/fermi3+getRho(r, nuc,false)/fermi3;
            k_fermi = pow(3*M_PI*M_PI*rho/2., 1./3.);
        }

        double corr_in_medium = nuc->Q_2p2h;

        Eout_lepton = lepton_mass+(Ein-lepton_mass-corr_in_medium)*frandom(); // corrected energy of the outgoing lepton

        w = Ein-Eout_lepton-corr_in_medium;

        theta_lepton = M_PI*frandom();  //generate the outgoing lepton angle
        double kout = sqrt(Eout_lepton*Eout_lepton-lepton_mass*lepton_mass);
        qq = sqrt(Ein*Ein + pow(kout,2) - 2*kout*Ein*cos(theta_lepton));

        if (qq<w) result_zero=true;
        //if (qq>1.2*GeV) result_zero=true; //we impose a cut for the transfered momentum

        phi_p = 2*M_PI*frandom(); //generate the phi angle of the 1st initial nucleons
        //generate the 4-momentum of "pion" :
        kk = (qq+3*k_fermi)*frandom();
        phi_pi = 2*M_PI*frandom();
        thetapi=M_PI*frandom();
        Epi = w*frandom();


        ulind(Epi/mpi, kk/mpi, k_fermi/mpi,cunuc, udel);
        cunuc*=pow(mpi,2);

        if (std::abs(cunuc.imag())<pow(10,-20)) result_zero=true; //if the pion cannot produce 1p1h (second bubble)

        if (!(qq<w) && (qq>1.2*GeV) && !(std::abs(cunuc.imag())<pow(10,-20))) aaaa++;

        double cosi=cos(thetapi);
        double kPi1=kk*sqrt(1-cosi*cosi)*cos(phi_pi);
        double kPi2 = kk*sqrt(1-cosi*cosi)*sin(phi_pi);
        double kPi3 = kk*cosi;

        double lindhard = DUR(w-Epi, sqrt(pow(kPi1,2)+pow(kPi2,2)+pow(kPi3-qq,2)),k_fermi, k_fermi);
        double q_minus_k=sqrt(pow(kPi1,2)+pow(kPi2,2)+pow(kPi3-qq,2));
        if (std::abs(lindhard)<pow(10,-20)) result_zero_lind = true;

        double Emax=sqrt(M*M+k_fermi*k_fermi);
        double pom = std::max(M,Emax-(w-Epi));
        double Emin=std::max(pom, 0.5*(-(w-Epi)+q_minus_k*sqrt(1.-4.*M*M/((w-Epi)*(w-Epi)-q_minus_k*q_minus_k))));

        if (Emin>Emax) result_zero= true;

        if (Emin>Emax and result_zero_lind==false) std::cerr<<"Wrong energy "<<Emax<<" "<<Emin<<"\n";
        double Emax_2 = sqrt(M*M+k_fermi*k_fermi);
        double pom_2 = std::max(M,Emax_2-Epi);
        double Emin_2 = std::max(pom_2, 0.5*(-(Epi)+kk*sqrt(1.-4.*M*M/(Epi*Epi-kk*kk))));

        Ep2 =  Emin_2+(Emax_2-Emin_2)*frandom(); //generate energy of the 2nd nucleon
        if (Emin_2>Emax_2) result_zero = true;

        Ep =Emin+(Emax-Emin)*frandom(); //generate energy of the 1st nucleon


        if (Emin>Emax) result_zero= true;
        if (Emin_2>Emax_2) result_zero = true;

        e->q0 = Ein-Eout_lepton;
        e->qz = qq;
        e->theta_lepton = theta_lepton;
        e->rho = rho;
        e->r=r;
        e->E_pi = Epi;
        e->theta_pi = thetapi;
        e->phi_pi=phi_pi;
        e->p_pi = kk;
        e->E_nuc =Ep;
        e->phi_nuc = phi_p;
        e->E_nuc2 = Ep2;
        e->energy_neutrino = Ein;

        if (!result_zero && !result_zero_lind) event_2p2h_rho_gaussian(r, Ein, theta_lepton, Ep2, Ep, phi_p, kk, phi_pi, thetapi, Epi, w, qq, rho, 1, Emax_2-Emin_2, e, nuc,corr_in_medium);
        else {
            e->weight_pp=0;
            e->weight_np=0;
            e->weight_pn=0;
            e->weight_anti_nn=0;
            e->weight_anti_pn=0;
            e->weight_anti_np=0;
        }

        double coeff = w*sin(theta_lepton)*(Emax-Emin)*(qq+3*k_fermi)*(Ein-lepton_mass-corr_in_medium)*Pi*2*Pi*2*Pi*Pi;

        if (!Global) coeff*=nuc->rmax;


        e->weight_pp *=coeff ;//range of integration
        e->weight_np *=coeff;
        e->weight_pn*=coeff;



        if ((e->weight_pp+e->weight_np+e->weight_pn)!=0) {
            tf->Fill(); //we do not include 0

            if (i*10./nr_events == counter) {std::cerr<<"* "; counter++;}
            i++;

            total_xsec+=(e->weight_pp+e->weight_np+e->weight_pn);


        }

        delete e;
        result_zero=false;
        result_zero_lind=false;

//        if (i*10./nr_events == counter) {std::cerr<<"* "; counter++;}
//        i++;
//
//        total_xsec+=(e->weight_pp+e->weight_np+e->weight_pn);


    }

    total_xsec*=1./nr_events*pow(10,38)/cm2*2*Pi;

    std::cerr<<"TOTAL X-sec: "<<total_xsec<<"\n";

    TVectorD constants(5);
    constants[0]=nr_events;
    constants[1]=Ein;
    constants[2]=lepton_mass;
    constants[4]=nuc->Q_2p2h;
    constants[3]=pow(10,38)/cm2*GeV;

    constants.Write("constants");

    ff->Write ();
    ff->Close ();
    delete ff;

}