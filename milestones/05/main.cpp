# include "../../../include/Molecular-Dynamics/verlet.h"
# include "../../../include/Molecular-Dynamics/xyz.h"
# include "../../../include/Molecular-Dynamics/energy.h"

int main() {
    // Initialize atoms parameters
    double mass = 1.0;
    double sigma = 1.0;
    double epsilon = 1.0;
    double time = sqrt(mass * pow(sigma,2) / epsilon);
    double total_time = 100* time;
    std::cout << "total_time: " << total_time << std::endl;
    double time_step = time/1000.0;
    std::cout << "time_step: " << time_step << std::endl;
    //    double relaxation_time = 50 * time_step; // relaxation time
    double relaxation_time = 10.0 * time_step; // relaxation time
    double desired_temperature = 10.0; // desired temperature



    // Initialize atoms on a cubic lattice
    double lattice_constant = sigma * 1.0;
    Atoms atoms = lattice(10, 10, 5, lattice_constant);
    Energy energy(atoms, epsilon, sigma, mass); // initialize energy class

    // Creating visualization file milestones/04/output
    std::cout << "checkpoint 1 " <<  std::endl;
    std::ofstream traj("/home/eslam/Desktop/Molecular-Dynamics/output/traj_2.xyz");

    // Main simulation loop
    for (int i = 0; i < total_time; ++i) {
        if (i % 10 == 0 ){
            write_xyz(traj, atoms);
//            relaxation_time = (i/10.0 +1) * time_step;
        }
        verlet_step1(atoms, time_step, mass);
        energy.update(atoms, epsilon, sigma, mass);
        verlet_step2(atoms, time_step, mass);
        // thermal bathing
        energy.berendsen_thermostat(atoms,  desired_temperature,  time_step, relaxation_time);
    }
    traj.close();
    std::cout << "Finished: Hoooooooooooooooray!!!!" << std::endl;

    return 0;
}