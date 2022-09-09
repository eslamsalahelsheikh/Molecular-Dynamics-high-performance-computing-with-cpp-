#ifndef MOLECULARDYNAMICS_SIMULATION_H
#define MOLECULARDYNAMICS_SIMULATION_H

#include "simulation_data.h"
#include "../../../include/Molecular-Dynamics/domain.h"
#include "mpi.h"

class Simulation : public SimulationData, public Energy {

public:
    // Constructors
    Simulation();
    Simulation(Atoms &new_atoms);
    // Destructor
    ~Simulation();
    // Methods
    void initial_loop(Domain domain);
    void relaxation_loop(int iteration);
    void add_heat();
private:
    // initialize members
    Atoms atoms_;
    NeighborList neighbor_list_;
    std::ofstream traj_;
};




#endif //MOLECULARDYNAMICS_SIMULATION_H
