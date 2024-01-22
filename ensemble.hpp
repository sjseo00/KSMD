
#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <CL/sycl.hpp>

#include "util.hpp"


float ComputePressure(queue &q, buffer<Atom> &atoms, int numAtoms, float volume);
float ComputeKineticEnergy(queue &q, buffer<Atom, 1> &atoms, int numAtoms);
void UpdateThermostat(queue &q, buffer<Atom> &atoms, int numAtoms, float dt, float T0, float tauT, Thermostat& thermostat);

    
#endif
