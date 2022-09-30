#ifndef MOLECULARDYNAMICS_UTILS_H
#define MOLECULARDYNAMICS_UTILS_H

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

/**
 * @brief  Export the energy and force of atoms to a file
 * @param iteration Iteration number
 * @param energy_file Energy file
 * @param total_energy  Total energy
 * @param average_temp  Average temperature
 * @param potential_energy Potential energy
 * @param continue_old_experiment Boolean if the simulation is continued from an old experiment
 */
inline void export_data(int iteration, std::ofstream &energy_file, double total_energy, double average_temp,
                        double potential_energy, bool &continue_old_experiment) {

    if (iteration == 0 or continue_old_experiment) {
        continue_old_experiment = false;
        //    TODO:: use relative paths
        energy_file << "iteration,total_energy,average_temp,potential_energy" << std::endl;
    }
    energy_file << iteration << "," << total_energy << "," << average_temp << "," << potential_energy << std::endl;
}
/**
 * @brief  Export the energy and force of atoms to a file
 * @param iteration Iteration number
 * @param energy_file Energy file
 * @param total_energy Total energy
 * @param potential_energy Potential energy
 */
inline void export_data(int iteration, std::ofstream &energy_file, double total_energy, double potential_energy) {
    if (iteration == 0) energy_file << "iteration,total_energy,potential_energy" << std::endl;
    energy_file << iteration << "," << total_energy << "," << potential_energy << std::endl;
}
/**
 * @brief  Export xyz file for relaxation loop
 * @param directory Directory to export the file to
 * @param iteration_num Iteration number
 * @param atoms Atoms object
 */
inline void export_xyz_relax(std::string directory, int iteration_num, Atoms &atoms) {
    std::ofstream traj3(directory + "/traj_" + std::to_string(atoms.nb_atoms()) + "_" + std::to_string(iteration_num) +
                        "_relax.xyz");
    write_xyz(traj3, atoms); // write xyz file every 10 steps        verlet_step1(atoms_, time_step, mass);
}
/**
 * @brief  Export xyz file for initial simulation loop
 * @param directory Directory to export the file to
 * @param iteration_num Iteration number
 * @param atoms Atoms object
 */
inline void export_xyz_initial(std::string directory, int iteration_num, Atoms &atoms) {
    std::ofstream traj2(directory + "/traj_" + std::to_string(iteration_num) + "_initial.xyz");
    write_xyz(traj2, atoms); // write xyz file every 10 steps        verlet_step1(atoms_, time_step, mass);
}


#endif //MOLECULARDYNAMICS_UTILS_H
