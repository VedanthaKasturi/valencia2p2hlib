//
// Created by Vedantha Srinivas Kasturi on 01.07.2025.
//

#ifndef VALENCIA2P2HLIB_TABLE_MC_H
#define VALENCIA2P2HLIB_TABLE_MC_H

#include <stdio.h>
#include "hadron_tensor_2p2h.h"
#include "hadron_tensor_2p2h_2b.h"
#include "hadron_tensor_2p2h_rho.h"
#include "Calka.h"
#include "utils.h"
#include "3p3h_MC.h"

extern int isospin_channel;

void hadron_tensor_2p2h_part1_part2_rho_MC ( double w, double qq, nucleus_str *nuc, double* w_tensor, double*  w_tensor_2, double* w_tensor_3);
void hadron_tensor_3p3h_MC ( double w, double qq, nucleus_str *nuc, double* w_tensor );

#endif //VALENCIA2P2HLIB_TABLE_MC_H
