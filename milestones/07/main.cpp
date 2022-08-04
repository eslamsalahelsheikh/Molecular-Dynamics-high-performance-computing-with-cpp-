#include "../../../include/Molecular-Dynamics/simulation.h"


int main() {
    // Reading initial positions and velocities from xyz file
    auto [names, positions]{read_xyz("/home/eslam/Desktop/Molecular-Dynamics/clusters/cluster_561.xyz")};
    // initializing atoms with poses
    Atoms atoms{positions};
    std::cout << "number of atoms: " << atoms.nb_atoms() << std::endl;
    // Initialize simulation for given atoms
    Simulation simulation(atoms);
    // Run simulation
    simulation.initial_loop(1000);

    std::vector<double> average_temps{};
    std::vector<double> potential_energies{};
    std::vector<double> total_energies{};

    int num_experiments = 30;
    for (int i = 0; i < num_experiments; i++) {
        std::cout << "--------------------starting relaxation loop number: "<< i << "---------------------------" << std::endl;
        // Deposit heat
        simulation.add_heat();
        // Run relaxation loop
        int relax_steps = 200;
        double average_temp = simulation.relaxation_loop(relax_steps);
        average_temps.push_back(average_temp);
        potential_energies.push_back(simulation.get_potential_energy());
        total_energies.push_back(simulation.get_total_energy());
    }
//    TODO:: use relative paths
    std::string directory = "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/"+ std::to_string(atoms.nb_atoms());
    std::filesystem::create_directory(directory);
    export_data(directory+"/total_energies_"+std::to_string(atoms.nb_atoms())+".txt", total_energies);
    export_data(directory+"/average_temps_"+std::to_string(atoms.nb_atoms())+".txt", average_temps);
    export_data(directory+"/potential_energies_"+std::to_string(atoms.nb_atoms())+".txt", potential_energies);
    return 0;
}