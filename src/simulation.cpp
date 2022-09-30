
#include "../../../include/Molecular-Dynamics/simulation.h"

// constructor
Simulation::Simulation() : Energy() {}

Simulation::Simulation(Atoms &new_atoms) : atoms_{new_atoms}, Energy(atoms_, epsilon, sigma, mass) {
    energy_update(atoms_, epsilon, sigma); // Update energy class
    neighbor_list_ = NeighborList(cutoff_radius); // NeighborList object for cutoff
    std::cout << "initialized Simulation parameters " << std::endl;
}

Simulation::~Simulation() {}

// Main simulation loop
void Simulation::initial_loop() {
    bool equilibrum = false;
    for (int i = 0; i < total_steps; ++i) {
        std::cout << "step: " << i << std::endl;
        if (equilibrum) std::cout << " equlibruim reached" << std::endl;
        std::cout << "initial steps: " << i << "  current_temp: " << get_temperature() << "  current_potential: "
                  << get_potential_energy() << "  total_energy: " << get_total_energy() << std::endl;
        verlet_step1(atoms_, time_step, mass);  // update positions
        neighbor_list_.update(atoms_);  // update neighbor list
        update_gupta(atoms_, neighbor_list_, cutoff_radius);    // update gupta potential
        verlet_step2(atoms_, time_step, mass);  // update velocities
        // thermal bathing
        // Apply the Berendsen thermostat for stop_thermostate_after_steps and then increase coupling constant to big value
        if (i == stop_thermostate_after_steps) {
            equilibrum = true;
            relaxation_time_multiplier = relaxation_time_multiplier_final_value;    // this should be big enough to reduce thermostat effect
        }
        relaxation_time = relaxation_time_multiplier * time_step; // relaxation time
        // apply thermostat
        berendsen_thermostat(atoms_, desired_temperature, time_step, relaxation_time);
        if (i % 10 == 0) export_xyz_initial(directory, i, atoms_); // export to xyz file
    }
}

void Simulation::initial_loop(Domain domain) {
    // Simulation loop
    for (int i = 0; i < total_steps; ++i) {
        std::cout << "step: " << i << std::endl;
        double old_total_energy = get_total_energy();   // store old total energy
        verlet_step1(atoms_, time_step, mass);  // update positions
        domain.exchange_atoms(atoms_);  // exchange atoms between domains after updating positions
        domain.update_ghosts(atoms_, cutoff_radius * 2); // update ghost atoms before calculating forces
        neighbor_list_.update(atoms_);  // update neighbor list
        update_gupta(atoms_, neighbor_list_, cutoff_radius, domain);    // update gupta potential
        verlet_step2(atoms_, time_step, mass);  // update velocities
        // thermal bathing
        if (i == stop_thermostate_after_steps)
            relaxation_time_multiplier = relaxation_time_multiplier_final_value;    // this should be big enough to reduce thermostat effect
        relaxation_time = relaxation_time_multiplier * time_step; // relaxation time
        berendsen_thermostat(atoms_, desired_temperature, time_step, relaxation_time);

        // Check if total energy is conserved
//        if (abs(get_total_energy() - old_total_energy)<0.01)
//        {
//            domain.disable(atoms_);
//            std::cout << "quilibrium reached, exiting now!!" << std::endl;
//            MPI_Finalize();
//            exit(1);
//        }
        if (i % 10 == 0) {
            domain.disable(atoms_); // disable domain decomposition before exporting to xyz file
            export_xyz_initial(directory, i, atoms_);   // export to xyz file
            export_data(i, energy_file, get_total_energy(), get_potential_energy()); // export energy data to file
            domain.enable(atoms_);  // enable domain decomposition
            domain.exchange_atoms(atoms_);  // exchange atoms between domains
            domain.update_ghosts(atoms_, cutoff_radius * 2);  // update ghost atoms
        }
    }
}

void Simulation::relaxation_loop(int iteration) {
    // Relaxation loop experiment
    std::cout << "--------------------starting relaxation loop number: " << iteration << "---------------------------"
              << std::endl;
    double total_temp = 0.0;
    for (int i = 0; i < relaxation_steps; ++i) {
        std::cout << "relaxation steps: " << i << "  current_temp: " << get_temperature() << "  current_potential: "
                  << get_potential_energy() << "  total_energy: " << get_total_energy() << std::endl;
        total_temp += get_temperature();    // Sum up temperature over all steps
        verlet_step1(atoms_, time_step, mass);  // update positions
        neighbor_list_.update(atoms_);  // update neighbor list
        update_gupta(atoms_, neighbor_list_, cutoff_radius);    // update gupta potential
        verlet_step2(atoms_, time_step, mass);  // update velocities
        if (i % 10 == 0) {
            export_xyz_relax(directory, iteration * relaxation_steps + i, atoms_);
        } // write xyz file every 10 steps
    }
    // adding energies to a file
    double average_temperature = total_temp / relaxation_steps; // calculate average temperature
    export_data(iteration, energy_file, get_total_energy(), average_temperature, get_potential_energy(),
                continue_old_experiment);   // export energy data to file
}

void Simulation::add_heat() {
    // Add heat to the system
    deposit_heat(atoms_, add_energy);
}

