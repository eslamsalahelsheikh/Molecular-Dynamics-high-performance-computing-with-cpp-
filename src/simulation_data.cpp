#include "../../../include/Molecular-Dynamics/simulation_data.h"

SimulationData::SimulationData() {
    // Initialize all simulation parameters
    cluster_file = "/home/eslam/Desktop/Molecular-Dynamics/milestones/08/cluster_923.xyz";
    milestone = "08";   // milestone number
    layer_numbers = 12; // number of layers
    atomic_distance = 2.885; // atomic distance
    mass = 196.9665 * 103.6; // atomic mass of Gold (https://www.nuclear-power.com/gold-atomic-number-mass-density/)
    total_steps = 5000; // total number of steps for simulation
    time_step = 1.0; // time step in fs
    cutoff_radius = 15.0;    // cutoff radius for EAM potential in Angstrom
    relaxation_time_multiplier = 10; // relaxation time = relaxation time multiplier * time_step in fs
    stop_thermostate_after_steps = 2000; // stop thermostat after this number of steps
    relaxation_time_multiplier_final_value = 1e10; // after the system arrives at desired temp (stop thermostat and relax system)
    desired_temperature = 1.0; // goal temperature (only in the start) in K

    // Relaxation experiment parameters
    relaxation_steps = 2000; // number of relaxation steps for simulation
    expermint_num = 40;    // number of experiments (Adding heat and relaxing)
    add_energy = 0.01;    // energy/atom to be added in each experiment

    // MPI parameters
    domain_length = {30.0, 30.0, 30.0}; // domain length in Angstrom
    domain_grid = {2, 2, 2}; // number of domains in each direction
    domain_periodicity = {0, 0, 0}; // periodicity of the domain in each direction

    // choose whether to continue old experiment or not
    continue_old_experiment = false;
    old_experiment_file = "traj_1415_4990_initial.xyz"; // name of the old experiment file (ONLY if continuing old experiment)

    create_directories_and_files(); // create directories and files
}

void SimulationData::create_directories_and_files() {
    // getting layer numbers and number of ranks if using domain decomposition
    std::string number_of_layers = std::to_string(layer_numbers);   // number of layers
    std::string number_of_cores = std::to_string(domain_grid.prod());   // number of ranks
    // creating directory for the experiment
    if (milestone == "07")
        directory =
                "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/" + number_of_layers +
                "/";
    else if (milestone == "08") directory = "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_08/923/";
    // creating energy file
    std::filesystem::create_directory(directory); // create directory
    energy_file = std::ofstream(directory + number_of_cores + "_energies.csv");
    // getting cluster file
    old_experiment_file = directory + old_experiment_file;
}