#ifndef __VERLET_H
#define __VERLET_H
# include <iostream>
# include "atoms.h"

void verlet_step1(Atoms &atoms, double timestep, double mass);
void verlet_step2(Atoms &atoms, double timestep, double mass);

#endif // __VERLET_H