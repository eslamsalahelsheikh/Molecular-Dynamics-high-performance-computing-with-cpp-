#ifndef __VERLET_H
#define __VERLET_H
#include <iostream>
#include <math.h>

#include "atoms.h"

void verlet_step1(Atoms &atoms, double timestep, double mass);
void verlet_step2(Atoms &atoms, double timestep, double mass);
void compute_force(Atoms &atoms, float sigma, float epsilon) ;
#endif // __VERLET_H