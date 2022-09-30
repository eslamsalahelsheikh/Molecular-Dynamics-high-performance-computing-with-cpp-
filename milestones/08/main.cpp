#include "../../../include/Molecular-Dynamics/simulation.h"
#include "../../../include/Molecular-Dynamics/ih.h"

int main(int argc, char *argv[]) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    // Initialize atoms class and getting simulation data
    Atoms atoms;
    SimulationData data;
    // Reading initial positions from cluster xyz file
    auto [names, positions]{read_xyz("/home/eslam/Desktop/Molecular-Dynamics/milestones/08/cluster_923.xyz")};
    atoms = Atoms{positions};
    // initialize Domain class for domain decomposition
    Domain domain(MPI_COMM_WORLD, data.domain_length, data.domain_grid, data.domain_periodicity);
    // Decompose atoms into domains
    domain.enable(atoms);
    // Exchange atoms between domains
    domain.exchange_atoms(atoms);
    // Update ghost atoms
    domain.update_ghosts(atoms, data.cutoff_radius * 2);
    // Initialize simulation class for given atoms
    Simulation simulation(atoms);

    // Run simulation
    simulation.initial_loop(domain);
    // Disable domain decomposition
    domain.disable(atoms);
    // Finalize MPI environment
    MPI_Finalize();
    return 0;
}