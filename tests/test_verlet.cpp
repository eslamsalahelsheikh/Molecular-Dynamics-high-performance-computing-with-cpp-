#include <gtest/gtest.h>

#include "../include/Molecular-Dynamics/atoms.h"
#include "../include/Molecular-Dynamics/verlet.h"


TEST(VerletTest, BasicAssertions) {
    int nbAtoms = 100;
    Positions_t  positions(3,nbAtoms);
    Velocities_t  velocities(3, nbAtoms);
    Forces_t forces(3,nbAtoms);
    double timestep = 0.1;
    double mass = 1.0;
    positions.setRandom();
    velocities.setRandom();
    Atoms atoms{positions, velocities};
    double total_steps = 1000;
    for (int i = 0; i < total_steps; ++i) {
        verlet_step1(atoms, timestep, mass);
        verlet_step2(atoms, timestep, mass);
    }
    for (int i = 0; i < nbAtoms; ++i) {
        EXPECT_NEAR(atoms.positions(0,i), positions(0,i) + velocities(0,i) * timestep * total_steps, 1e-10);
        EXPECT_NEAR(atoms.positions(1,i), positions(1,i) + velocities(1,i) * timestep * total_steps, 1e-10);
        EXPECT_NEAR(atoms.positions(2,i), positions(2,i) + velocities(2,i) * timestep * total_steps, 1e-10);
        EXPECT_NEAR(atoms.velocities(0,i), velocities(0,i), 1e-10);
        EXPECT_NEAR(atoms.velocities(1,i), velocities(1,i), 1e-10);
        EXPECT_NEAR(atoms.velocities(2,i), velocities(2,i), 1e-10);
    }
}
