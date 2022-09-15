#include "../../../include/Molecular-Dynamics/simulation.h"
#include "../../../include/Molecular-Dynamics/ih.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    Atoms atoms;
    SimulationData data;
    if (data.continue_old_experiment) {
        auto [names, positions, Velocities_t]{read_xyz_with_velocities(data.old_experiment_file)};
        atoms = Atoms{positions, Velocities_t};
    } else {
//        Positions_t positions = generate_cluster(data.layer_numbers, data.atomic_distance);
        auto [names, positions]{read_xyz("/home/eslam/Desktop/Molecular-Dynamics/milestones/08/cluster_923.xyz")};

        atoms = Atoms{positions};
    }

    Domain domain(MPI_COMM_WORLD, data.domain_length, data.domain_grid, data.domain_periodicity);
    domain.enable(atoms);
    // Initialize simulation for given atoms
    Simulation simulation(atoms);

    // Run initial simulation
    if (!data.continue_old_experiment) simulation.initial_loop(domain);
    else simulation.energy_update(atoms);

    domain.disable(atoms);
    MPI_Finalize();
    return 0;
}