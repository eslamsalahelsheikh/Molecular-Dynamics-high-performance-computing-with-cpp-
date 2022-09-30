#ifndef MOLECULARDYNAMICS_SIMULATION_H
#define MOLECULARDYNAMICS_SIMULATION_H

#include "simulation_data.h"
#include "../../../include/Molecular-Dynamics/domain.h"
#include "mpi.h"

class Simulation : public SimulationData, public Energy {

public:
    // Constructors

    /**
     * @brief Default constructor
     */
    Simulation();

    /**
     * @brief Construct a new Simulation object
     * @param new_atoms Atoms class object
     */
    Simulation(Atoms &new_atoms);
    // Destructor

    /**
     * @brief Destroy the Simulation object
     */
    ~Simulation();

    // Methods

    /**
     * @brief Run simulation until equilibrium is reached
     */
    void initial_loop();

    /**
     * @brief Run simulation until equilibrium is reached using domain decomposition and MPI parallelization
     * @param domain Domain class object
     */
    void initial_loop(Domain domain);

    /**
     * @brief After adding heat, run the simulation for some time and compute the average temperature during this time
     * @param iteration The number of the current experiment
     */
    void relaxation_loop(int iteration);

    /**
     * @brief Adding heat to the system
     */
    void add_heat();

private:
    // Private attributes
    Atoms atoms_;
    NeighborList neighbor_list_;
    std::ofstream traj_;
};

#endif //MOLECULARDYNAMICS_SIMULATION_H
