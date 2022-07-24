#include "../../../include/Molecular-Dynamics/energy.h"
#include "../../../include/Molecular-Dynamics/verlet.h"
#include "../../../include/Molecular-Dynamics/xyz.h"

int main() {
    // Initialize atoms parameters
    double mass = 1.0;
    double sigma = 1.0;
    double epsilon = 1.0;
    double time = sqrt(mass * pow(sigma, 2) / epsilon);
    double total_time = 5000 * time;
    std::cout << "total_time: " << total_time << std::endl;
    double time_step = time / 1000.0;
    std::cout << "time_step: " << time_step << std::endl;
    double relaxation_time_multiplier = 1.0;
    double relaxation_time;     // relaxation time
    double desired_temperature; // desired temperature

    // Initialize atoms on a cubic lattice
    double lattice_constant = sigma * 0.8;
    Atoms atoms = lattice(3, 3, 3, lattice_constant);
    Energy energy(atoms, epsilon, sigma, mass); // initialize energy class

    // Creating visualization file milestones/04/output
    std::cout << "checkpoint 1 " << std::endl;
    std::ofstream traj(
        "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_05/traj.xyz");

    // Main simulation loop
    for (int i = 0; i < total_time; ++i) {
        std::cout << "step: " << i << std::endl;
        if (i % 10 == 0) {
            write_xyz(traj, atoms);
        }
        //        if (abs(desired_temperature - energy.get_temperature()) <=
        //        0.01 && energy.get_temperature() != 0.0) {
        //            relaxation_time_multiplier += 0.1;
        //        }
        //        else {
        //            relaxation_time_multiplier <= 1.0 ?
        //            relaxation_time_multiplier = 1.0 :
        //            relaxation_time_multiplier -= 0.1;
        //        }
        double old_total_energy = energy.get_total_energy();
        verlet_step1(atoms, time_step, mass);
        energy.update(atoms, epsilon, sigma, mass);
        verlet_step2(atoms, time_step, mass);
        // thermal bathing
        // to preserve the temperature, we assume that the temperature is constant, and then we apply the Berendsen thermostat
        desired_temperature =
            i == 0 ? energy.get_temperature() : desired_temperature;
        if (abs(energy.get_total_energy() - old_total_energy) <= 0.01) { // reached equilibrium point, decrease coupling constant
            relaxation_time_multiplier = 10;
        }
        relaxation_time =   relaxation_time_multiplier * time_step; // relaxation time
        energy.berendsen_thermostat(atoms, desired_temperature, time_step,
                                    relaxation_time);
    }
    traj.close();
    std::cout << "Finished: Hoooooooooooooooray!!!!" << std::endl;

    return 0;
}