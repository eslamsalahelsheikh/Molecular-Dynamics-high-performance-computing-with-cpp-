#ifndef MOLECULARDYNAMICS_SIMULATION_DATA_H
#define MOLECULARDYNAMICS_SIMULATION_DATA_H

#include "energy.h"
#include "verlet.h"
#include "xyz.h"
#include "utils.h"

class SimulationData {
public:
    /**
     * @brief Construct a new Simulation Data object
     */
    SimulationData();

    // Simulation parameters
    double expermint_num; // number of experiments
    bool continue_old_experiment; // continue old experiment or not
    std::string old_experiment_file; // name of the old experiment file
    std::ofstream energy_file; // file to write energies to
    std::string directory; // directory to write files to
    std::string cluster_file; // name of the cluster file
    int layer_numbers; // number of layers in the cluster
    double atomic_distance; // distance between atoms in the cluster
    double cutoff_radius;   // cutoff radius for EAM potential in Angstrom
    Eigen::Array<double, 3, 1> domain_length; // domain length in Angstrom
    Eigen::Array<int, 3, 1> domain_grid; // number of domains in each direction
    Eigen::Array<int, 3, 1> domain_periodicity; // periodicity of the domain in each direction

protected:
    // Protected attributes
    double mass;    // mass of Gold
    double sigma;   // distance at which the potential is a minimum
    double epsilon; // depth of the potential well
    double total_steps; // total number of steps for simulation
    double time_step;   // time step in fs
    double relaxation_time_multiplier;  // relaxation time = relaxation time multiplier * time_step in fs
    double relaxation_time_multiplier_final_value;  // after the system arrives at desired temp (stop thermostat and relax system)
    double relaxation_time;     // relaxation time
    double desired_temperature; // desired temperature
    double relaxation_steps; // steps for relaxation after adding heat
    double add_energy; // energy added in each experiment
    double stop_thermostate_after_steps; // stop thermostat after this number of steps

private:
    // Private methods
    /**
     * @brief Create directories and files to store data
     */
    void create_directories_and_files();

    // Private attributes
    std::string milestone; // milestone number
};

#endif //MOLECULARDYNAMICS_SIMULATION_DATA_H
