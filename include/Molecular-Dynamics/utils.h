#ifndef MOLECULARDYNAMICS_UTILS_H
#define MOLECULARDYNAMICS_UTILS_H
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

//inline void export_data(std::string filename, std::vector<double> data){
//    std::ofstream file;
//    file.open(filename);
//    for (auto &i : data) {
//        file << i << std::endl;
//    }
//    file.close();
//}
inline void export_data(int iteration, std::ofstream &energy_file, double total_energy, double average_temp, double potential_energy, bool &continue_old_experiment){

    if (iteration == 0 or continue_old_experiment) {
        continue_old_experiment = false;
        //    TODO:: use relative paths
        energy_file << "iteration,total_energy,average_temp,potential_energy" << std::endl;
    }
    energy_file << iteration << "," << total_energy << "," << average_temp << "," << potential_energy << std::endl;
}
inline void export_data(int iteration, std::ofstream &energy_file, double total_energy,double potential_energy){
    energy_file << "iteration,total_energy,potential_energy" << std::endl;
    energy_file << iteration << "," << total_energy << "," << "," << potential_energy << std::endl;
}

inline void export_xyz_relax(std::string directory, int iteration_num, Atoms &atoms){
    std::ofstream traj3 (directory+"/traj_"+std::to_string(atoms.nb_atoms())+"_"+std::to_string(iteration_num)+"_relax.xyz");
    write_xyz(traj3, atoms); // write xyz file every 10 steps        verlet_step1(atoms_, time_step, mass);
}
inline void export_xyz_initial(std::string directory, int iteration_num, Atoms &atoms){
    std::ofstream traj2 (directory+"/traj_"+std::to_string(iteration_num)+"_initial.xyz");
    write_xyz(traj2, atoms); // write xyz file every 10 steps        verlet_step1(atoms_, time_step, mass);
}



#endif //MOLECULARDYNAMICS_UTILS_H
