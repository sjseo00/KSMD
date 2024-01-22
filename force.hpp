#ifndef FORCE_H
#define FORCE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <CL/sycl.hpp>
#include <math.h>

#include "util.hpp"

void BondKernel(queue &q, buffer<Atom> &atoms, buffer<Bond> &bonds, EnvInfo env_info, SimInfo sim_info);
void AngleKernel(queue &q, buffer<Atom> &atoms, buffer<Angle> &angles, int numAtoms, int numAngles);
void DihedralKernel(queue &q, buffer<Atom> &atoms, buffer<Dihedral> &dihedrals, int numAtoms, int numDihedrals);
void LennardJonesKernel(queue &q, buffer<Atom> &atoms, EnvInfo env_info, SimInfo sim_info);
void EwaldKernel(queue &q, buffer<Atom> &atoms, int numAtoms, EnvInfo info);


#endif
