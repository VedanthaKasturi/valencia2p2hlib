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
void saveEvents_2p2h(bool Global, int oxygen, bool efixed, double en_low, double en_high, double nr, std::string suffix);
void saveEvents_3p3h(std::string suffix);
void saveEvents_3p3h(bool Global, int oxygen, double en, double nr, std::string suffix);
void total_2p2h(double Emin, double Emax, double nr_d, std::string name);
void total_2p2h(double Emin, double Emax, double nr_d, bool Global, int oxygen, std::string suffix);
void total_3p3h(double Emin, double Emax, double nr_d, std::string suffix);
void total_3p3h(double Emin, double Emax, double nr_d, bool Global, int oxygen, std::string suffix);
void generate_table_2p2h_MC(std::string suffix);

int main(int argc, char *argv[]) {
    if (argc < 10) {
        std::cerr << "Usage: ./program <mode_2p2h> <mode_el> <ReDelta> <emin> <emax> <nr_d> <Global> <oxygen> <efixed> <en_low> <en_high> <nr>\n";
        std::cerr << "mode_2p2h: 0 - save 2p2h events, 1 - save 3p3h events, 2 - total 2p2h, 3 - total 3p3h, 4 - Inclusive 2p2h, 5 - save hadronic tensor (MC method)\n";
        std::cerr << "mode_el: 0 - CC, 1 - CCanti\n";
        std::cerr << "ReDelta: 0 - 0 MeV, 1 - -50 rho/rho0 MeV, 2 - 30 rho/rho0 MeV\n";
        std::cerr << "emin: minimum energy of incoming neutrino\n";
        std::cerr << "emax: maximum energy of incoming neutrino\n";
        std::cerr << "nr_d: number of points to generate\n";
        std::cerr << "Global: 0 - no, 1 - yes\n";
        std::cerr << "oxygen: 0 - C12, 1 - O16\n";
        std::cerr << "efixed: 0 - no, 1 - yes\n";
        std::cerr << "nr: number of points to generate\n";
        std::cerr << "en_low: low energy of incoming neutrino\n";
        std::cerr << "en_high: high energy of incoming neutrino\n";
        return 1;
    }

    // Print out the options selected

    std::cout<<"Selected options:\n";
    std::cout<<"mode_2p2h: "<<argv[1]<<"\n";
    std::cout<<"mode_el: "<<argv[2]<<"\n";
    std::cout<<"ReDelta: "<<argv[3]<<"\n";
    std::cout<<"emin: "<<argv[4]<<"\n";
    std::cout<<"emax: "<<argv[5]<<"\n";
    std::cout<<"nr_d: "<<argv[6]<<"\n";
    std::cout<<"Global: "<<argv[7]<<"\n";
    std::cout<<"oxygen: "<<argv[8]<<"\n";
    std::cout<<"efixed: "<<argv[9]<<"\n";
    std::cout<<"en_low: "<<argv[10]<<"\n";
    std::cout<<"en_high: "<<argv[11]<<"\n";
    std::cout<<"nr: "<<argv[12]<<"\n";


    M = 940*MeV;
    mpi = 139.57*MeV;

    mpip=mpi;  //pi+- mass
    mpi0=mpi;  //pi0 mass
    mpi2=mpi*mpi;
    mpi4=mpi2*mpi2;
    DP_const=delusive*fstar/mpip*sqrt3;
    CDP_const=delusive*fstar/mpip/sqrt3;

    lepton_mass=mmu;

    // Parse command-line arguments
    int mode_2p2h = std::stoi(argv[1]);
    int mode_el = std::stoi(argv[2]);
    int rede = std::stoi(argv[3]);
    double emin_delta = std::stod(argv[4]);
    double emax_delta = std::stod(argv[5]);
    int nr_d = std::stoi(argv[6]);
    bool Global = std::stoi(argv[7]);
    int oxygen = std::stoi(argv[8]);
    bool efixed = std::stoi(argv[9]);

    double en_low = 0;
    double en_high = 0;
    double nr = 1;
    if (argc >= 12) {
        en_low = std::stod(argv[10]);
        if (argc >= 13) {
            en_high = std::stod(argv[11]);
        }
        nr = std::stod(argv[12]);
    }

    CH_2p2h = CC;
    if (mode_el == 1) CH_2p2h = CCanti;

    std::string suffix = "_ReDelta_0";
    switch (rede) {
        case 0:
            ReDelta = 0;
            suffix = "_ReDelta_0";
            break;
        case 1:
            ReDelta = -50;
            suffix = "_ReDelta_m50";
            break;
        case 2:
            ReDelta = 30;
            suffix = "_ReDelta_30";
            break;
    }

    switch (mode_2p2h) {
        case 0:
            saveEvents_2p2h(Global, oxygen, efixed, en_low, en_high, nr, suffix);
            break;
        case 1:
            saveEvents_3p3h(Global, oxygen, en_low, nr, suffix);
            break;
        case 2:
            total_2p2h(emin_delta, emax_delta, nr_d, Global, oxygen, suffix);
            break;
        case 3:
            total_3p3h(emin_delta, emax_delta, nr_d, Global, oxygen, suffix);
            break;
        case 4:
            total_2p2h(emin_delta, emax_delta, nr_d, Global, oxygen, suffix);
            break;
        case 5:
            generate_table_2p2h_MC(suffix);
            break;
    }

    return 0;
}


void saveEvents_3p3h(bool Global, int oxygen, double en, double nr, std::string suffix) {
    frandom_init(0);

    nucleus_str nuc;
    if (oxygen == 1 && !Global) {
        NR_PART = 16;
        nuc = O16;
    } else if (oxygen == 1 && Global) {
        NR_PART = 16;
        nuc = GFG_16;
    } else if (oxygen == 0 && !Global) {
        NR_PART = 12;
        nuc = C12;
    } else if (oxygen == 0 && Global) {
        NR_PART = 12;
        nuc = GFG_12;
    }

    en *= MeV;

    std::ostringstream name;
    name << "events_3p3h_CC=" << CH_2p2h << "_" << nuc.name << "_E=" << en << "_" << suffix << ".root";

    saveEvents_3p3h(en, name.str(), nr, &nuc);
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

void saveEvents_2p2h(bool Global, int oxygen, bool efixed, double en_low, double en_high, double nr, std::string suffix) {
    frandom_init(0);
    nucleus_str nuc;
    std::string gfg_;
//    CH_2p2h=CC;
    FERMI=true;
    NR_PART=12;
    Global = false;
    gfg_="LFG";
    nuc=C12;

    if (oxygen == 1 && !Global) {
        NR_PART = 16;
        nuc = O16;
    } else if (oxygen == 1 && Global) {
        NR_PART = 16;
        nuc = GFG_16;
    } else if (oxygen == 0 && !Global) {
        NR_PART = 12;
        nuc = C12;
    } else if (oxygen == 0 && Global) {
        NR_PART = 12;
        nuc = GFG_12;
    }

    if (efixed) {
        en_low *= MeV;

        std::ostringstream name_rho;
        name_rho << "events_CC=" << CH_2p2h << "_" << nuc.name << "_E=" << en_low<<"_"<<suffix << ".root";

        saveEvents_2p2h(en_low, name_rho.str(), nr, &nuc);
    } else {
        en_low *= MeV;
        en_high *= MeV;

        std::ostringstream name_rho;
        name_rho << "events_CC=" << CH_2p2h << "_" << nuc.name << "_E=" << en_low << "_" << en_high<<suffix << ".root";

        std::cout<<"Generating events for energy range: "<<en_low/MeV<<" - "<<en_high/MeV<<" MeV\n";
        std::cout<<"Generating "<<nr<<"*10^5 events\n";
        std::cout<<"nucleus: "<<nuc.name<<"\n";
        std::cout<<"Global Fermi Gas: "<<Global<<"\n";
        std::cout<<"filename: "<<name_rho.str()<<"\n";
        std::cout<<"Fixed energy: "<<efixed<<"\n";
        std::cout<<"CC: "<<CH_2p2h<<"\n";


        saveEvents_2p2h_variedEnergy(en_low, en_high, name_rho.str(), nr, &nuc);
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


void total_2p2h(double Emin, double Emax, double nr_d, bool Global, int oxygen, std::string suffix) {
    std::ofstream f;
    frandom_init(0);

    nucleus_str nuc;
    if (oxygen == 1 && !Global) {
        NR_PART = 16;
        nuc = O16;
    } else if (oxygen == 1 && Global) {
        NR_PART = 16;
        nuc = GFG_16;
    } else if (oxygen == 0 && !Global) {
        NR_PART = 12;
        nuc = C12;
    } else if (oxygen == 0 && Global) {
        NR_PART = 12;
        nuc = GFG_12;
    }

    std::ostringstream name;
    name << "total_2p2h_CC=" << CH_2p2h << "_" << nuc.name << "_" << suffix << ".txt";

    f.open(name.str());
    double res;
    for (int i = 0; i <= nr_d; i++) {
        double E = Emin + (Emax - Emin) * i / nr_d;
        if (E > 1000) nr_d = 5;
        total_xsec_2p2h(E, nr_d, &nuc, res);
        f << E << " " << res << "\n";
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

void total_3p3h(double Emin, double Emax, double nr_d, bool Global, int oxygen, std::string suffix) {
    std::ofstream f;
    frandom_init(0);

    nucleus_str nuc;
    if (oxygen == 1 && !Global) {
        NR_PART = 16;
        nuc = O16;
    } else if (oxygen == 1 && Global) {
        NR_PART = 16;
        nuc = GFG_16;
    } else if (oxygen == 0 && !Global) {
        NR_PART = 12;
        nuc = C12;
    } else if (oxygen == 0 && Global) {
        NR_PART = 12;
        nuc = GFG_12;
    }

    std::ostringstream name;
    name << "total_3p3h_CC=" << CH_2p2h << "_" << nuc.name << "_" << suffix << ".txt";

    f.open(name.str());
    double res;
    for (int i = 0; i <= nr_d; i++) {
        double E = Emin + (Emax - Emin) * i / nr_d;
        if (E > 1000) nr_d = 2;
        if (E > 1500) nr_d = 5;
        res = total_xsec_3p3h(E, nr_d, &nuc);
        f << E << " " << res << "\n";
        f.flush();
    }
    f.close();
}





