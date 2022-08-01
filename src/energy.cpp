#include "../include/Molecular-Dynamics/energy.h"
#include <iostream>

Energy::Energy() {
}
Energy::Energy(Atoms &atoms, double epsilon, double sigma, double mass) {
    kinetic_energy_ = 0.0;
    potential_energy_ = 0.0;
    total_energy_ = 0.0;
    temperature_ = 0.0;
    energy_update(atoms, epsilon, sigma, mass);
}

double Energy::kinetic_energy(Atoms &atoms, double mass) {
    double KE = 0.0;
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        KE += 0.5 * mass * pow(atoms.velocities.col(i).matrix().norm(), 2);
    }
    return KE;
}

double Energy::potential_energy(Atoms &atoms, double epsilon, double sigma) {
    return lj_direct_summation(atoms, epsilon, sigma);
}
double Energy::potential_energy_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma){
    return lj_direct_summation_neighbors(atoms, neighbor_list, epsilon, sigma);
}
double Energy::total_energy() {
    return potential_energy_ + kinetic_energy_;
}
double Energy::temperature(Atoms &atoms) {
    return 2.0 / 3.0 * (kinetic_energy_ / atoms.nb_atoms()); // assuming kB=1
}

void Energy::energy_update(Atoms &atoms, double epsilon, double sigma, double mass) {
    kinetic_energy_ = kinetic_energy(atoms, mass);
    potential_energy_ = potential_energy(atoms, epsilon, sigma);
    total_energy_ = total_energy();
    temperature_ = temperature(atoms);
//    std::cout << "total energy: " << total_energy_
//              << " temperature: " << temperature_ << std::endl;
}

void Energy::berendsen_thermostat(Atoms &atoms, double desired_temperature,
                                  double timestep, double relaxation_time) {
//     std::cout << "relaxation time: " << relaxation_time << std::endl;
    double current_temperature =
        temperature(atoms); // recalculate temperature after second verlet step
    double temperature_ratio = desired_temperature / current_temperature;
    double scaling_factor =
        sqrt(1 + ((temperature_ratio - 1) * timestep / relaxation_time));
    atoms.velocities *= scaling_factor; // scaling velocities
}

void Energy::update_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma, double mass) {
    kinetic_energy_ = kinetic_energy(atoms, mass);
    potential_energy_ = potential_energy_neighbors(atoms, neighbor_list, epsilon, sigma);
    total_energy_ = total_energy();
    temperature_ = temperature(atoms);
}

void Energy::update_gupta(Atoms &atoms, NeighborList &neighbor_list,double mass,double cutoff_radius) {
    kinetic_energy_ = kinetic_energy(atoms, mass);
    potential_energy_ = gupta(atoms, neighbor_list,cutoff_radius);
    total_energy_ = total_energy();
    temperature_ = temperature(atoms);
}

void Energy::deposit_heat(Atoms &atoms, double heat) {
    atoms.velocities *= sqrt(1. + (heat*atoms.nb_atoms()) / kinetic_energy_);
}