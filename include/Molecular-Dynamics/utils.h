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
inline void export_data(std::string directory, int atoms_num, double total_energy, double average_temp, double potential_energy){
    std::string total_energy_filename = directory +"/total_energies_"+std::to_string(atoms_num)+".txt";
    std::string average_temp_filename = directory +"/average_temps_"+std::to_string(atoms_num)+".txt";
    std::string potential_energy_filename = directory +"/potential_energies_"+std::to_string(atoms_num)+".txt";
    std::ofstream file;
    file.open(total_energy_filename);
    file << total_energy << std::endl;
    file.close();

    file.open(average_temp_filename);
    file << average_temp << std::endl;
    file.close();

    file.open(potential_energy_filename);
    file << potential_energy << std::endl;
    file.close();
}


inline void export_xyz_relax(int iteration_num, Atoms &atoms){
    std::string directory = "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/"+ std::to_string(atoms.nb_atoms());
    std::filesystem::create_directory(directory);
    std::ofstream traj3 (directory+"/traj_"+std::to_string(atoms.nb_atoms())+"_"+std::to_string(iteration_num)+"_relax.xyz");
    write_xyz(traj3, atoms); // write xyz file every 10 steps        verlet_step1(atoms_, time_step, mass);
}
inline void export_xyz_initial(int iteration_num, Atoms &atoms){
    std::string directory = "/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/"+ std::to_string(atoms.nb_atoms());
    std::filesystem::create_directory(directory);
    std::ofstream traj2 (directory+"/traj_"+std::to_string(atoms.nb_atoms())+"_"+std::to_string(iteration_num)+"_initial.xyz");
    write_xyz(traj2, atoms); // write xyz file every 10 steps        verlet_step1(atoms_, time_step, mass);
}



#endif //MOLECULARDYNAMICS_UTILS_H
