#ifndef MOLECULARDYNAMICS_UTILS_H
#define MOLECULARDYNAMICS_UTILS_H
#include <vector>
#include <string>
#include <fstream>

inline void export_data(std::string filename, std::vector<double> data){
    std::ofstream file;
    file.open(filename);
    for (auto &i : data) {
        file << i << std::endl;
    }
    file.close();
}

#endif //MOLECULARDYNAMICS_UTILS_H
