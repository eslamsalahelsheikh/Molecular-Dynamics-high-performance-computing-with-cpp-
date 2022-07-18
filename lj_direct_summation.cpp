
# include "lj_direct_summation.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double total_energy = 0.0;
    for (int i = 0; i < atoms.nb_atoms(); i++) {
        for (int j = i + 1; j < atoms.nb_atoms(); j++) {
            Eigen::Array3d r{atoms.positions.col(j) - atoms.positions.col(i)};
            Eigen::Array3d r_squared = r.array().square();
            double r_squared_sum = r_squared.sum();
            double r_squared_sqrt = std::sqrt(r_squared_sum);
            double r_squared_sqrt_sigma = r_squared_sqrt / sigma;
            double r_squared_sqrt_sigma_6 = pow(r_squared_sqrt_sigma, 6);
            double r_squared_sqrt_sigma_12 = pow(r_squared_sqrt_sigma, 12);
            double V = 4.0 * epsilon * (r_squared_sqrt_sigma_12 - r_squared_sqrt_sigma_6);
            total_energy += V;
        }
    }
    return total_energy;
}