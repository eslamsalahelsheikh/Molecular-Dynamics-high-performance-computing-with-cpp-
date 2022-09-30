#include "../include/Molecular-Dynamics/energy.h"
#include <iostream>

Energy::Energy() {}; // default constructor

Energy::Energy(Atoms &atoms, double epsilon, double sigma, double mass) {
    kinetic_energy_ = 0.0;
    potential_energy_ = 0.0;
    total_energy_ = 0.0;
    temperature_ = 0.0;
    mass_ = mass;
}

double Energy::kinetic_energy(Atoms &atoms) {
    double KE = 0.0;
    // loop over all atoms to calculate the kinetic energy
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        KE += 0.5 * mass_ * pow(atoms.velocities.col(i).matrix().norm(), 2);
    }
    return KE;
}

double Energy::kinetic_energy(Atoms &atoms, double local_atoms_num) {
    double KE = 0.0;
    // loop over the local atoms to calculate the kinetic energy
    for (int i = 0; i < local_atoms_num; ++i) {
        KE += 0.5 * mass_ * pow(atoms.velocities.col(i).matrix().norm(), 2);
    }
    return KE;
}

double Energy::potential_energy(Atoms &atoms, double epsilon, double sigma) {
    return lj_direct_summation(atoms, epsilon, sigma); // calculate the potential energy using the direct summation
}

double Energy::potential_energy_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma) {
    return lj_direct_summation_neighbors(atoms, neighbor_list, epsilon,
                                         sigma); // calculate the potential energy using the direct summation with neighbor list
}

double Energy::total_energy() {
    return potential_energy_ + kinetic_energy_; // calculate the total energy
}

double Energy::calculate_temperature(Atoms &atoms, bool use_exist_kinetic) {
    double KE = use_exist_kinetic ? kinetic_energy_ : kinetic_energy(atoms);
    const double BOLTZMANN_CONSTANT{8.617333262e-5}; // eV/K
    if (atoms.nb_atoms() > 0)
        return 2.0 / 3.0 * (KE / (BOLTZMANN_CONSTANT * atoms.nb_atoms())); // calculate the temperature
    else return 0.0;
}

void Energy::energy_update(Atoms &atoms, double epsilon, double sigma) {
    kinetic_energy_ = kinetic_energy(atoms);    // update the kinetic energy
    potential_energy_ = potential_energy(atoms, epsilon,
                                         sigma);    // update the potential energy using the direct summation
    total_energy_ = total_energy(); // update the total energy
    temperature_ = calculate_temperature(atoms, kinetic_energy_);   // update the temperature
}

void Energy::berendsen_thermostat(Atoms &atoms, double desired_temperature, double timestep, double relaxation_time) {
    if (atoms.nb_atoms() <= 0) return;
    double current_temperature = calculate_temperature(atoms,
                                                       false); // recalculate temperature after second verlet step
    double temperature_ratio = desired_temperature / current_temperature; // calculate the temperature ratio
    double scaling_factor = sqrt(
            1 + ((temperature_ratio - 1) * timestep / relaxation_time)); // calculate the scaling factor
    atoms.velocities *= scaling_factor; // scaling velocities
}

void Energy::update_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma, double mass) {
    kinetic_energy_ = kinetic_energy(atoms); // update kinetic energy
    potential_energy_ = potential_energy_neighbors(atoms, neighbor_list, epsilon,
                                                   sigma); // update potential energy using neighbor list
    total_energy_ = total_energy(); // update total energy
    temperature_ = calculate_temperature(atoms, true); // update temperature
}

void Energy::update_gupta(Atoms &atoms, NeighborList &neighbor_list, double cutoff_radius) {
    kinetic_energy_ = kinetic_energy(atoms);    // update the kinetic energy
    potential_energy_ = gupta(atoms, neighbor_list, cutoff_radius);  // update the potential energy using gupta
    total_energy_ = total_energy(); // update the total energy
    temperature_ = calculate_temperature(atoms, true);  // update the temperature
}

void Energy::update_gupta(Atoms &atoms, NeighborList &neighbor_list, double cutoff_radius, Domain &domain) {
    double per_atom_potential_energy = 0.0;
    double per_atom_kinetic_energy = 0.0;
    if (atoms.nb_atoms() > 0) {
        // calculate potential energy and kinetic energy for only the local atoms in this rank
        per_atom_potential_energy = gupta(domain.nb_local(), atoms, neighbor_list,
                                          cutoff_radius);    // update the potential energy using gupta
        per_atom_kinetic_energy = kinetic_energy(atoms, domain.nb_local());  // update the kinetic energy
    }
    // summing up potential & kinetic energy over all ranks
    potential_energy_ = MPI::allreduce(per_atom_potential_energy, MPI_SUM, MPI_COMM_WORLD);
    kinetic_energy_ = MPI::allreduce(per_atom_kinetic_energy, MPI_SUM, MPI_COMM_WORLD);
    total_energy_ = total_energy(); // update the total energy
    temperature_ = calculate_temperature(atoms, true);  // update the temperature
    std::cout << "total energy: " << total_energy_ << " temperature: " << temperature_ << " potential energy: "
              << potential_energy_ << std::endl;
}

void Energy::deposit_heat(Atoms &atoms, double added_energy) {
    if (atoms.nb_atoms() <= 0) return;
    std::cout << "before adding heat, temp is: " << get_temperature() << std::endl;
    double total_added_energy = added_energy * atoms.nb_atoms(); // adding energy per atom
    atoms.velocities *= sqrt(1. + total_added_energy / kinetic_energy_); // scaling velocities
    std::cout << "after adding heat, temp is: " << calculate_temperature(atoms, false) << std::endl;
}