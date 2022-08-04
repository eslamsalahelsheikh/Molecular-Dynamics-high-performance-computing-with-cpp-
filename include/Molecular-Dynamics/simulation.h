#ifndef MOLECULARDYNAMICS_SIMULATION_H
#define MOLECULARDYNAMICS_SIMULATION_H

#include "../../../include/Molecular-Dynamics/energy.h"
#include "../../../include/Molecular-Dynamics/verlet.h"
#include "../../../include/Molecular-Dynamics/xyz.h"
#include "../../../include/Molecular-Dynamics/utils.h"
#include <filesystem>

class SimulationData {
protected:
    double mass;
    double sigma;
    double epsilon;
    double total_steps;
    double time_step;
    double cutoff_radius;
    double relaxation_time_multiplier;
    double relaxation_time;     // relaxation time
    double desired_temperature; // desired temperature
    SimulationData() {
        // Initialize atoms parameters
        mass = 196.9665; // atomic mass of Gold (https://www.nuclear-power.com/gold-atomic-number-mass-density/)
        sigma = 1.0;
        epsilon = 1.0;
//        total_steps = 5000;
        time_step = 1.0; // time step in fs
        cutoff_radius = 10.0;    // cutoff radius for LJ potential
        relaxation_time_multiplier = 10.0; // relaxation time in fs
    }
};

class Simulation : public SimulationData, public Energy {

public:
    // Constructor
    Simulation(Atoms &new_atoms);
    // Destructor
    ~Simulation();
    // Methods
    void initial_loop(double total_steps);
    double relaxation_loop(double total_steps);
    void add_heat();
private:
    // initialize members
    Atoms atoms_;
    NeighborList neighbor_list_;
    std::ofstream traj_;
};




#endif //MOLECULARDYNAMICS_SIMULATION_H
