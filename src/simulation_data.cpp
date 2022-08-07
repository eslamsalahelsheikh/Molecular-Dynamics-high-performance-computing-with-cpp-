#include "../../../include/Molecular-Dynamics/simulation_data.h"

SimulationData::SimulationData() {
    // Initialize all simulation parameters
    cluster_name = "cluster_923";

    mass = 196.9665* 103.6; // atomic mass of Gold (https://www.nuclear-power.com/gold-atomic-number-mass-density/)
    total_steps = 5000;
    time_step = 10; // time step in fs
    cutoff_radius = 5.0;    // cutoff radius for EAM potential
    relaxation_time_multiplier = 10; // relaxation time = relaxation time multiplier * time_step in fs
    stop_thermostate_after_steps = 2000; // stop thermostat after this number of steps
    relaxation_time_multiplier_final_value = 1e10; // after the system arrives at desired temp (stop thermostat and relax system)
    desired_temperature = 500; // desired temperature (only in the start) in K

    // Relaxation experiment parameters
    relaxation_steps = 500; // number of relaxation steps
    expermint_num = 500;    // number of experiments
    add_energy = 0.01;    // energy added in each experiment

    // choose whether to continue old experiment or not
    continue_old_experiment = false;
    old_experiment_file = directory + "traj_923_114_relax.xyz"; // name of the old experiment file (ONLY if continuing old experiment)

    create_directories_and_files();
}

void SimulationData::create_directories_and_files() {
    // getting number of atoms from the name of the cluster
    std::size_t pos = cluster_name.find("_");      // position of "_" in str
    std::string number_of_atoms = cluster_name.substr (pos+1);     // get from "live" to the end
    // creating directory for the experiment
    directory = "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/"+number_of_atoms+"/";
    std::filesystem::create_directory(directory);
    // creating energy file
    energy_file = std::ofstream(directory +  number_of_atoms + "_energies.csv");
    // getting cluster file
    cluster_file = "/home/eslam/Desktop/Molecular-Dynamics/clusters/"+cluster_name+".xyz";
}