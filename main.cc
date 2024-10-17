#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <sys/types.h>
#include "generatormt.h"

#include "1p1h1pion.h"
#include "hadron_tensor_2p2h.h"
#include "hadron_tensor_2p2h_2b.h"
#include "hadron_tensor_2p2h_rho.h"
#include "3p3h_MC.h"
#include "2p2h_MC.h"


void saveEvents_2p2h(std::string suffix);
void saveEvents_3p3h(std::string suffix);
void total_2p2h(double Emin, double Emax, double nr_d, std::string name);
void total_3p3h(double Emin, double Emax, double nr_d, std::string suffix);

int main(int argc, char *argv[])

{
	
	M = 940*MeV;
	mpi = 139.57*MeV;
	
	mpip=mpi;  //pi+- mass
	mpi0=mpi;  //pi0 mass
	mpi2=mpi*mpi;
	mpi4=mpi2*mpi2;
	DP_const=delusive*fstar/mpip*sqrt3;
	CDP_const=delusive*fstar/mpip/sqrt3;

	lepton_mass=mmu;
	
	int mode_2p2h;
	std::cout<<"CHOOSE MODE\n";
	std::cout<<"(0) Save 2p2h events to a ROOT file \n";
    std::cout<<"(1) Save 3p3h events to a ROOT file \n";
	std::cout<<"(2) Calculate total xsection for 2p2h\n";
    std::cout<<"(3) Calculate total xsection for 3p3h\n";

	std::cin>>mode_2p2h;
    
    CH_2p2h=CC;
    int mode_el = 0;
    std::cout<<"Which process: CC(0) CCanti(1)\n";
    std::cin>>mode_el;
    if (mode_el==1) CH_2p2h=CCanti;
    
    std::string suffix = "_ReDelta_0";
    ReDelta=0.;
    int rede = 0;
    std::cout<<"Delta self-energy: 0 MeV(0), -50 rho/rho0 MeV(1), 30 rho/rho0 MeV(2)\n";
    std::cin>>rede;
    switch(rede) {
        case 0:
            ReDelta=0; suffix = "_ReDelta_0"; break;
        case 1:
            ReDelta=-50; suffix = "_ReDelta_m50";  break;
        case 2:
            ReDelta=30; suffix = "_ReDelta_30";  break;
    }
    
    
	double emin_delta, emax_delta;
	int nr_d;
	switch (mode_2p2h) {
		case 0:
			saveEvents_2p2h(suffix);
			break;
        case 1:
            saveEvents_3p3h(suffix);
            break;
		case 2:
			std::cout<<"Emin [Mev]= ";
			std::cin>>emin_delta;
			std::cout<<"Emax= ";
			std::cin>>emax_delta;
			std::cout<<"NR of steps= ";
			std::cin>>nr_d;
			total_2p2h(emin_delta, emax_delta,nr_d,suffix);
			break;
        case 3:
            std::cout<<"Emin [Mev]= ";
            std::cin>>emin_delta;
            std::cout<<"Emax= ";
            std::cin>>emax_delta;
            std::cout<<"NR of steps= ";
            std::cin>>nr_d;
            total_3p3h(emin_delta, emax_delta,nr_d,suffix);
            break;
	}
	
	return 0;
	
}



void saveEvents_3p3h(std::string suffix) {
	
	frandom_init(0);
	
	
	std::string gfg_;
	std::cout<<"Global Fermi Gas? 0-no 1-yes\n";
	std::cin>>Global;
	
	int oxygen;
	nucleus_str nuc;
	std::cout<<"Nucleus: C12 (0), O16 (1)\n";
	std::cin>>oxygen;
	
	
    if (oxygen==1 && !Global) {
		NR_PART=16;
		nuc = O16;
	}
	else if (oxygen==1 && Global) {
		NR_PART=16;
		nuc = GFG_16;
	}
	else if (oxygen==0 && !Global){
		NR_PART=12;
		nuc=C12;
	}
	else if (oxygen==0 && Global){
		NR_PART=12;
		nuc=GFG_12;
	}
	
	
	
	double en;
	double nr=1;
	std::cout<<"E= [MeV]\n";
	std::cin>>en;
	en*=MeV;
	
	std::cout<<"How many points generate? x(*10^6)\n";
	std::cin>>nr;
	
	std::ostringstream name;
	
	name<<"events_3p3h_CC="<<CH_2p2h<<"_"<<nuc.name<<"_E="<<en<<"_"<<suffix<<".root";
	saveEvents_3p3h(en, name.str(), nr, &nuc);

}

void saveEvents_2p2h(std::string suffix){
	
	frandom_init(0);
	std::string gfg_;
	int oxygen;
	nucleus_str nuc;
	
	CH_2p2h=CC;
	FERMI=true;
	NR_PART=12;
	Global = false;
	gfg_="LFG";
	nuc=C12;
	
	std::cout<<"Global Fermi Gas? 0-no 1-yes\n";
    std::cin>>Global;
		
    std::cout<<"Nucleus:  C12 (0), O16 (1)\n";
    std::cin>>oxygen;
		
		
    if (oxygen==1 && !Global) {
			NR_PART=16;
			nuc = O16;
		}
		else if (oxygen==1 && Global) {
			NR_PART=16;
			nuc = GFG_16;
		}
		else if (oxygen==0 && !Global){
			NR_PART=12;
			nuc=C12;
		}
		else if (oxygen==0 && Global){
			NR_PART=12;
			nuc=GFG_12;
		}

    bool efixed = false;

    std::cout<<"Fixed energy? 0-no 1-yes\n";
    std::cin>>efixed;

    double nr=1;


   if(efixed) {
       double en;
       std::cout<<"Energy of incoming neutrino [MeV]\n";
       std::cin>>en;
       en*=MeV;


       std::cout<<"How many points generate? x(*3*10^6)\n";
       std::cin>>nr;


       std::ostringstream name_rho;
       name_rho<<"events_CC="<<CH_2p2h<<"_"<<nuc.name<<"_E="<<en<<"_"<<suffix<<".root";

       saveEvents_2p2h(en, name_rho.str(), nr, &nuc);
    }
    else {
        double enl;
        double enh;

        std::cout<<"Energy low of incoming neutrino [MeV]\n";
        std::cin>>enl;
        enl*=MeV;

        std::cout<<"Energy high of incoming neutrino [MeV]\n";
        std::cin>>enh;
        enh*=MeV;

       std::cout<<"How many points generate? x(*3*10^6)\n";
       std::cin>>nr;

       saveEvents_2p2h_variedEnergy(enl, enh, suffix, nr, &nuc);
    }
}

void total_2p2h(double Emin, double Emax, double nr_d, std::string suffix) {
	
	std::ofstream f;
	
	frandom_init(0);
	int oxygen;
	nucleus_str nuc;

	NR_PART=12;
	nuc=C12;
	
		
    std::cout<<"Global Fermi Gas? 0-no 1-yes\n";
    std::cin>>Global;
		
		
    std::cout<<"Nucleus:    C12 (0), O16 (1)\n";
    std::cin>>oxygen;
		
		
    if (oxygen==1 && !Global) {
        NR_PART=16;
        nuc = O16;
    }
    else if (oxygen==1 && Global) {
        NR_PART=16;
        nuc = GFG_16;
    }
    else if (oxygen==0 && !Global){
        NR_PART=12;
        nuc=C12;
    }
    else if (oxygen==0 && Global){
        NR_PART=12;
        nuc=GFG_12;
    }
		
	
	
	double nr = 1.;
    /*
	std::cout<<"How many points generate? x(*10^6)\n";
	std::cin>>nr;
	*/
	
	std::ostringstream name;
	
	name<<"total_2p2h_CC="<<CH_2p2h<<"_"<<nuc.name<<"_"<<suffix<<".txt";
	
	f.open(name.str());
	double res;
	for (int i=0; i<=nr_d; i++) {
		double E = Emin + (Emax-Emin)*i/nr_d;
        if (E>1000) nr=5;
		total_xsec_2p2h(E, nr, &nuc, res);
		f<<E<<" "<<res<<"\n";
		f.flush();
		
	}
	f.close();
	
	
}




void total_3p3h( double Emin, double Emax, double nr_d, std::string suffix) {
    
    std::ofstream f;
    
    frandom_init(0);
    int oxygen;
    nucleus_str nuc;
    
    NR_PART=12;
    nuc=C12;

        
        
    std::cout<<"Global Fermi Gas? 0-no 1-yes\n";
    std::cin>>Global;
    
        
    std::cout<<"Nucleus:    C12 (0), O16 (1)\n";
    std::cin>>oxygen;
        
        if (oxygen==1 && !Global) {
            NR_PART=16;
            nuc = O16;
        }
        else if (oxygen==1 && Global) {
            NR_PART=16;
            nuc = GFG_16;
        }
        else if (oxygen==0 && !Global){
            NR_PART=12;
            nuc=C12;
        }
        else if (oxygen==0 && Global){
            NR_PART=12;
            nuc=GFG_12;
        }
        
    
    
    double nr = 1.;
    
    std::ostringstream name;
    
    name<<"total_3p3h_CC="<<CH_2p2h<<"_"<<nuc.name<<"_"<<suffix<<".txt";
    
    f.open(name.str());
    double res;
    for (int i=0; i<=nr_d; i++) {
        
        double E = Emin + (Emax-Emin)*i/nr_d;
        if (E>1000) nr=2;
        if (E>1500) nr=5;
        res = total_xsec_3p3h(E, nr, &nuc);
        f<<E<<" "<<res<<"\n";
        f.flush();
        
    }
    f.close();
    
    
}






