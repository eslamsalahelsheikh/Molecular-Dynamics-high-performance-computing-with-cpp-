#ifndef MOLECULAR_DYNAMICS_ENERGY_H
#define MOLECULAR_DYNAMICS_ENERGY_H

# include "atoms.h"
# include "lj_direct_summation.h"

class Energy {
public:
    // Constructors
    Energy();
    Energy(Atoms &atoms, double epsilon, double sigma, double mass);
    // Methods
    void berendsen_thermostat(Atoms &atoms, double desired_temperature, double timestep, double relaxation_time);
    void update(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0, double mass = 1.0);
    void update_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon = 1.0, double sigma = 1.0, double mass = 1.0);
    // getters
    double get_potential_energy() { return potential_energy_;}
    double get_kinetic_energy() {return kinetic_energy_;}
    double get_total_energy() {return total_energy_;}
    double get_temperature() {return temperature_;}
    double temperature(Atoms &atoms);
private:
    // Private attributes
    double epsilon;
    double sigma;
    double kinetic_energy_;
    double potential_energy_;
    double total_energy_;
    double temperature_;
    // Private methods
    double kinetic_energy(Atoms &atoms, double mass);
    double potential_energy(Atoms &atoms, double epsilon, double sigma);
    double potential_energy_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma);
    double total_energy();

};


#endif //MOLECULAR_DYNAMICS_ENERGY_H
