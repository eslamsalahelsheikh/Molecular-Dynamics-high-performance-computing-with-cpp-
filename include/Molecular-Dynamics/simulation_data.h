#ifndef MOLECULARDYNAMICS_SIMULATION_DATA_H
#define MOLECULARDYNAMICS_SIMULATION_DATA_H

#include "energy.h"
#include "verlet.h"
#include "xyz.h"
#include "utils.h"

class SimulationData {
public:
    SimulationData();
    double expermint_num; // number of experiments
    bool continue_old_experiment; // continue old experiment or not
    std::string old_experiment_file; // name of the old experiment file
    std::ofstream energy_file; // file to write energies to
    std::string directory; // directory to write files to
    std::string cluster_file; // name of the cluster file
protected:
    double mass;
    double sigma;
    double epsilon;
    double total_steps;
    double time_step;
    double cutoff_radius;
    double relaxation_time_multiplier;
    double relaxation_time_multiplier_final_value;
    double relaxation_time;     // relaxation time
    double desired_temperature; // desired temperature
    double relaxation_steps; // steps for relaxation after adding heat
    double add_energy; // energy added in each experiment
    double stop_thermostate_after_steps; // stop thermostat after this number of steps
private:
    void create_directories_and_files();
    std::string cluster_name; // name of the cluster

};

#endif //MOLECULARDYNAMICS_SIMULATION_DATA_H
