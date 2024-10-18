#Usage: ./program <mode_2p2h> <mode_el> <ReDelta> <emin> <emax> <nr_d> <Global> <oxygen> <efixed> <en_low> <en_high> <nr>

#mode_2p2h: 0 - save 2p2h events, 1 - save 3p3h events, 2 - total 2p2h, 3 - total 3p3h
#mode_el: 0 - CC, 1 - CCanti
#ReDelta: 0 - 0 MeV, 1 - -50 rho/rho0 MeV, 2 - 30 rho/rho0 MeV
#emin: minimum energy of incoming neutrino
#emax: maximum energy of incoming neutrino
#nr_d: number of points to generate
#Global: 0 - no, 1 - yes
#oxygen: 0 - C12, 1 - O16
#efixed: 0 - no, 1 - yes
#nr: number of points to generate
#en_low: low energy of incoming neutrino
#en_high: high energy of incoming neutrino

./build/program 0 0 0 0 2000 200 0 0 0 0 2000 100