#ifndef MOLECULAR_DYNAMICS_ENERGY_H
#define MOLECULAR_DYNAMICS_ENERGY_H

# include "atoms.h"
# include "lj_direct_summation.h"
# include "gupta.h"
# include "domain.h"

class Energy {
public:
    // Constructors

    /**
    * @brief Default constructor
    */
    Energy();

    /**
    * @brief Construct a new Energy object
    * @param atoms Atoms object
    * @param epsilon Depth of the potential well
    * @param sigma Distance at which the potential is a minimum
    * @param mass mass of the atoms
    */
    Energy(Atoms &atoms, double epsilon, double sigma, double mass);

    // Methods

    /**
    * @brief Apply the Berendsen thermostat to the system
    * @param atoms Atoms object
    * @param desired_temperature Goal temperature
    * @param timestep Timestep
    * @param relaxation_time Relaxation time
    */
    void berendsen_thermostat(Atoms &atoms, double desired_temperature, double timestep, double relaxation_time);

    /**
    * @brief Update kinetic, potential and total energy using the direct summation
    * @param atoms Atoms object
    * @param epsilon Depth of the potential well
    * @param sigma Distance at which the potential is a minimum
    */
    void energy_update(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

    /**
    * @brief Update kinetic, potential and total energy using the direct summation with a neighbor list
    * @param atoms Atoms object
    * @param neighbor_list NeighborList object
    * @param epsilon Depth of the potential well
    * @param sigma Distance at which the potential is a minimum
    * @param mass mass of the atoms
    */
    void update_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon = 1.0, double sigma = 1.0,
                          double mass = 1.0);

    /**
    * @brief Update kinetic, potential and total energy using Gupta potential
    * @param atoms Atoms object
    * @param neighbor_list NeighborList object
    * @param cutoff_radius cutoff radius for Gupta potential
    */
    void update_gupta(Atoms &atoms, NeighborList &neighbor_list, double cutoff_radius = 10.0);

    /**
    * @brief Update kinetic, potential and total energy using Gupta potential for just the local atoms in this rank
    * @param atoms Atoms object
    * @param neighbor_list NeighborList object
    * @param cutoff_radius cutoff radius for Gupta potential
    * @param domain Domain object
    */
    void update_gupta(Atoms &atoms, NeighborList &neighbor_list, double cutoff_radius, Domain &domain);

    /**
    * @brief Add heat to the system
    * @param atoms Atoms object
    * @param added_energy Amount of energy/atom to be added to the system
    */
    void deposit_heat(Atoms &atoms, double added_energy);

    /**
    * @brief Calculate the kinetic energy of the system
    * @param atoms Atoms object
    * @param use_exist_kinetic boolean to use the existing kinetic energy or not
    * @return the current temperature of the system
    */
    double calculate_temperature(Atoms &atoms, bool use_exist_kinetic);

    // getters
    double get_potential_energy() { return potential_energy_; }

    double get_kinetic_energy() { return kinetic_energy_; }

    double get_total_energy() { return total_energy_; }

    double get_temperature() { return temperature_; }

private:
    // Private attributes
    double kinetic_energy_;
    double potential_energy_;
    double total_energy_;
    double temperature_;
    double mass_;

    // Private methods

    /**
    * @brief Calculate the kinetic energy of the system
    * @param atoms Atoms object
    * @return the current kinetic energy of the system
    */
    double kinetic_energy(Atoms &atoms);

    /**
    *
    * @brief Calculate the kinetic energy of the system for just the local atoms in this rank
    * @param atoms Atoms object
    * @param local_atoms_num   number of local atoms in this rank
    * @return The current kinetic energy of the local atoms in this rank
    */
    double kinetic_energy(Atoms &atoms, double local_atoms_num);

    /**
    * @brief Calculate the potential energy of the system
    * @param atoms Atoms object
    * @param epsilon Depth of the potential well
    * @param sigma Distance at which the potential is a minimum
    * @return The current potential energy of the system
    */
    double potential_energy(Atoms &atoms, double epsilon, double sigma);

    /**
    * @brief Calculate the potential energy of the system using the direct summation with a neighbor list
    * @param atoms Atoms object
    * @param neighbor_list NeighborList object
    * @param epsilon Depth of the potential well
    * @param sigma Distance at which the potential is a minimum
    * @return The current potential energy of the system
    */
    double potential_energy_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma);

    /**
    * @brief Calculate the total energy of the system
    * @return The current total energy of the system
    */
    double total_energy();
};


#endif //MOLECULAR_DYNAMICS_ENERGY_H