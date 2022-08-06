#include "../../../include/Molecular-Dynamics/simulation.h"


int main() {
    // Reading initial positions and velocities from xyz file
    Atoms atoms;
    SimulationData data;
    if (data.continue_old_experiment) {
        auto [names, positions, Velocities_t]{read_xyz_with_velocities(data.old_experiment_file)};
        atoms = Atoms{positions, Velocities_t};
    } else {
        auto [names, positions]{read_xyz("/home/eslam/Desktop/Molecular-Dynamics/clusters/"+data.cluster_name+".xyz")};
        atoms = Atoms{positions};
    }
    // initializing atoms with poses
    std::cout << "number of atoms: " << atoms.nb_atoms() << std::endl;
    // Initialize simulation for given atoms
    Simulation simulation(atoms);

    // Run initial simulation
    if (!data.continue_old_experiment) simulation.initial_loop();

    std::vector<double> average_temps{};
    std::vector<double> potential_energies{};
    std::vector<double> total_energies{};
    //    TODO:: use relative paths
    std::string directory = "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/"+ std::to_string(atoms.nb_atoms());
    std::filesystem::create_directory(directory);

    simulation.energy_update(atoms);
    for (int i = 0; i < data.expermint_num; i++) {
        std::cout << "--------------------starting relaxation loop number: "<< i << "---------------------------" << std::endl;
        // Deposit heat
        simulation.add_heat();
        // Run relaxation loop
        double average_temp = simulation.relaxation_loop(i);
        average_temps.push_back(average_temp);
        potential_energies.push_back(simulation.get_potential_energy());
        total_energies.push_back(simulation.get_total_energy());
        // adding energies to a file
        export_data(directory, atoms.nb_atoms(), simulation.get_total_energy(), average_temp, simulation.get_potential_energy());
    }
//
//    export_data(directory+"/total_energies_"+std::to_string(atoms.nb_atoms())+".txt", total_energies);
//    export_data(directory+"/average_temps_"+std::to_string(atoms.nb_atoms())+".txt", average_temps);
//    export_data(directory+"/potential_energies_"+std::to_string(atoms.nb_atoms())+".txt", potential_energies);
    return 0;
}