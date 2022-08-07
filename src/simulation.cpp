
#include "../../../include/Molecular-Dynamics/simulation.h"

Simulation::Simulation(): Energy() {}
Simulation::Simulation(Atoms &new_atoms) : atoms_{new_atoms}, Energy(atoms_,epsilon,sigma,mass) {
    energy_update(atoms_, epsilon, sigma); // Update energy class
    neighbor_list_ = NeighborList(cutoff_radius); // NeighborList object for cutoff
    // Creating visualization file
    std::string directory = "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/"+ std::to_string(new_atoms.nb_atoms());
    std::filesystem::create_directory(directory);
//    traj_ = std::ofstream(directory+"/traj_"+std::to_string(new_atoms.nb_atoms())+".xyz");
    std::cout << "initialized Simulation parameters " << std::endl;
}
Simulation::~Simulation() {
    traj_.close();   // closing output file
}

// Main simulation loop
void Simulation::initial_loop() {
    bool equilibrum = false;
    for (int i = 0; i < total_steps; ++i) {
        if (equilibrum) std::cout << " equlibruim reached" << std::endl;
        std::cout << "initial steps: " << i << "  current_temp: " << get_temperature() << "  current_potential: " << get_potential_energy() << "  total_energy: " << get_total_energy() << std::endl;
        if (i % 10 == 0) {export_xyz_initial(i, atoms_);}
        neighbor_list_.update(atoms_);
        verlet_step1(atoms_, time_step, mass);
        update_gupta(atoms_, neighbor_list_, cutoff_radius);
        verlet_step2(atoms_, time_step, mass);
        // thermal bathing
        // to preserve the temperature, we assume that the temperature is constant, and then we apply the Berendsen thermostat
//        desired_temperature =
//                i == 0 ? get_temperature() : desired_temperature;
//        if (abs(get_total_energy() - old_total_energy) <= 0.01)  // reached equilibrium point, decrease coupling constant
//        if (i %100 == 0 && i != 0 && i <2000 ) { desired_temperature += 0.2 ;}
        if (i == 500)
        {
            equilibrum = true;
            relaxation_time_multiplier = relaxation_time_multiplier_final_value;    // this should be big enough to reduce thermostat effect
        }
        relaxation_time =   relaxation_time_multiplier * time_step; // relaxation time
        berendsen_thermostat(atoms_, desired_temperature, time_step,
                                    relaxation_time);
    }
}

double Simulation::relaxation_loop(int iteration) {
    double total_temp = 0.0;
    for (int i = 0; i < relaxation_steps; ++i) {
        std::cout << "relaxation steps: " << i << "  current_temp: " << get_temperature() << "  current_potential: " << get_potential_energy() << "  total_energy: " << get_total_energy() << std::endl;
        total_temp += get_temperature();
        if (i % 10 == 0) {export_xyz_relax(iteration*relaxation_steps+i, atoms_);} // write xyz file every 10 steps
        neighbor_list_.update(atoms_);
        verlet_step1(atoms_, time_step, mass);
        update_gupta(atoms_, neighbor_list_, cutoff_radius);
        verlet_step2(atoms_, time_step, mass);
    }
    return total_temp / total_steps;
}

void Simulation::add_heat(){
    deposit_heat(atoms_,add_energy);
}

