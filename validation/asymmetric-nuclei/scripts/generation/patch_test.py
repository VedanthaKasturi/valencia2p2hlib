import sys, re, pathlib
root = pathlib.Path(sys.argv[1]); which = sys.argv[2]
def sub(path, old, new, count=1):
    p = root/path; s = p.read_text()
    n = s.count(old)
    assert n == count, f"{path}: expected {count} matches, got {n} for {old[:60]!r}"
    p.write_text(s.replace(old, new))
# 1) seed from env MC_SEED
sub("src/generatormt.cpp", "case 0: init_genrand(time(NULL));break;",
    'case 0: { const char* s=getenv("MC_SEED"); init_genrand(s? strtoul(s,0,10) : time(NULL)); } break;')
s=(root/"src/generatormt.cpp").read_text()
if "<cstdlib>" not in s: (root/"src/generatormt.cpp").write_text("#include <cstdlib>\n"+s)
# 2) Q override from env Q2P2H_OVERRIDE (MeV)
if which == "new":
    sub("main.cc", "    if (CH_2p2h == CCanti) sel.Q_2p2h = sel.Q_2p2h_anti;\n",
        "    if (CH_2p2h == CCanti) sel.Q_2p2h = sel.Q_2p2h_anti;\n    if (const char* qo=getenv(\"Q2P2H_OVERRIDE\")) sel.Q_2p2h = atof(qo);\n")
else:
    sub("main.cc", "    return *n;\n}", "    nucleus_str sel=*n;\n    if (const char* qo=getenv(\"Q2P2H_OVERRIDE\")) sel.Q_2p2h = atof(qo);\n    std::cout << \"Q used: \" << sel.Q_2p2h << \"\\n\";\n    return sel;\n}")
s=(root/"main.cc").read_text()
if "<cstdlib>" not in s: (root/"main.cc").write_text("#include <cstdlib>\n"+s)
# 3) new only: per-piece weights
if which == "new":
    sub("lib/event.h", "    double kF_n;           // neutron Fermi momentum at r\n",
        "    double kF_n;           // neutron Fermi momentum at r\n    double w_rho[3]={0,0,0}; // TEST: rho-transverse piece of the slot weights\n    double w_pi[3]={0,0,0};  // TEST: pion-longitudinal piece\n    double w_int[3]={0,0,0}; // TEST: interference piece\n")
    # keep copies of the pieces before they are summed
    sub("src/2p2h_MC.cpp", """	for (int i=0; i<3; i++) {
		A00[i]+=B00[i]+C00[i];""", """	std::complex<double> P00[3][3],P11[3][3],P12[3][3],P03[3][3],P33[3][3]; //TEST: [piece][slot]
	for (int i=0; i<3; i++) {
		P00[0][i]=A00[i];P11[0][i]=A11[i];P12[0][i]=A12[i];P03[0][i]=A03[i];P33[0][i]=A33[i];
		P00[1][i]=B00[i];P11[1][i]=B11[i];P12[1][i]=B12[i];P03[1][i]=B03[i];P33[1][i]=B33[i];
		P00[2][i]=C00[i];P11[2][i]=C11[i];P12[2][i]=C12[i];P03[2][i]=C03[i];P33[2][i]=C33[i];
	}
	for (int i=0; i<3; i++) {
		A00[i]+=B00[i]+C00[i];""", count=2)
    # compute piece weights at the end of event_2p2h_rho_gaussian: insert before the switch on CCanti inside the slot loop (first occurrence = event function)
    s=(root/"src/2p2h_MC.cpp").read_text()
    anchor = "        if (CH_2p2h == CCanti) {\n            \n            switch(i){"
    assert s.count(anchor)==1, s.count(anchor)
    piece = """        for (int pc=0; pc<3; pc++) { //TEST: same contraction for each piece of the hadron tensor
            double pW1=P11[pc][i].real()/2.;
            double pW2=0.5*(P00[pc][i].real()+P11[pc][i].real()+ w*w/qq/qq*(P33[pc][i].real()-P11[pc][i].real())-2*w/qq*P03[pc][i].real());
            double pW3=P12[pc][i].imag()/qq;
            double pW4=(P33[pc][i].real()-P11[pc][i].real())/2./qq/qq;
            double pW5=(P03[pc][i].real()-w/qq*(P33[pc][i].real()-P11[pc][i].real()))/qq;
            double sg = (CH_2p2h==CCanti) ? -1. : 1.;
            double ps1 = 2*pW1*sint2 + pW2*cost2 - sg*pW3*(Ein+eout)*sint2;
            double ps2 = pW1*cos(theta_lepton)-pW2/2.*cos(theta_lepton) + sg*pW3/2.*(eout+kout-(Ein+eout)*cos(theta_lepton)) + pW4/2.*(lepton_mass*lepton_mass*cos(theta_lepton)+2*eout*(eout+kout)*(1-cos(theta_lepton)*cos(theta_lepton))) - pW5*(eout+kout)/2.;
            double pw = coef*(ps1+as*ps2)*factor_rho;
            if (pc==0) e->w_rho[i]=pw; else if (pc==1) e->w_pi[i]=pw; else e->w_int[i]=pw;
        }
"""
    (root/"src/2p2h_MC.cpp").write_text(s.replace(anchor, piece+anchor))
    sub("src/2p2h_MC.cpp", "e->weight_anti_np *= coeff;", "e->weight_anti_np *= coeff; for(int k=0;k<3;k++){e->w_rho[k]*=coeff; e->w_pi[k]*=coeff; e->w_int[k]*=coeff;}", count=3)
print("patched", which)
