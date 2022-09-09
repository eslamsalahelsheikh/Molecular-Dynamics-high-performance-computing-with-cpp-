/*
* Copyright 2021 Lars Pastewka
*
* ### MIT license
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#ifndef YAMD_DOMAIN_H
#define YAMD_DOMAIN_H

#include <mpi.h>

#include "atoms.h"
#include "mpi_support.h"

class Domain {
public:
    Domain(const MPI_Comm &comm, const Eigen::Array3d &domain_length, const Eigen::Array3i &decomposition,
           const Eigen::Array3i &periodicity);

    ~Domain();

    /*
     * Enable domain decomposition: After a call to this method, every process
     * remains only with the atoms in its local domain.
     */
    void enable(Atoms &atoms);

    /*
     * Disable domain decomposition: After this call to this method, all
     * processes contain identical copies of the Atoms object.
     */
    void disable(Atoms &atoms);

    /*
     * Communicate atoms that have left the local domain to the neighboring
     * domains. Ghost buffers will be invalidated after a call to this method.
     */
    void exchange_atoms(Atoms &atoms);

    /*
     * Communicate atoms into the ghost buffers of neighboring cells.
     */
    void update_ghosts(Atoms &atoms, double border_width);
    
    /*
     * Set new domain length and (affinely) rescale atom positions.
     */
    void scale(Atoms &atoms, Eigen::Array3d domain_length);

    /*
     * Return MPI communicator.
     */
    MPI_Comm communicator() { return comm_; }

    /*
     * Return rank of the current process.
     */
    int rank() const { return rank_; }

    /*
     * Return size of the communicator group.
     */
    int size() const { return size_; }

    /*
     * Return whether domain decomposition is enabled.
     */
    bool is_enabled() const { return is_enabled_; }

    /*
     * Return length of the global domain.
     */
    Eigen::Array3d domain_length() const { return domain_length_; }

    /*
     * Return length of the global domain in dimension dim.
    */
    double domain_length(int dim) const { return domain_length_(dim); }

    /*
     * Coordinate of local process.
     */
    Eigen::Array3i coordinate() const { return coordinate_; }

    /*
     * Coordinate of local process in dimension dim.
     */
    int coordinate(int dim) const { return coordinate_(dim); }

    /*
     * Number of process-local atoms.
     */
    int nb_local() const { return nb_local_; }

protected:
    /*
     * Get domain index given a position
     */
    Eigen::Array3i get_coordinate(const Eigen::Array3d &position) {
        return (position * decomposition_.cast<double>() / domain_length_).floor().cast<int>();
    }

    /*
     * Get domain index given an array of position
     */
    Eigen::ArrayXi get_coordinates(const Eigen::Array3Xd &positions, int dim) {
        return (positions.row(dim) * (static_cast<double>(decomposition_(dim) / domain_length_(dim))))
            .floor().cast<int>();
    }

    /*
     * Check whether domain decomposition is enabled and raise an error is this is not the case.
     */
    void assert_enabled() const {
        if (!is_enabled_) throw std::runtime_error("Expected domain decomposition to be enabled.");
    }

    /*
     * Check whether domain decomposition is disable and raise an error is this is not the case.
     */
    void assert_disabled() const {
        if (is_enabled_) throw std::runtime_error("Expected domain decomposition to be disabled.");
    }

    /*
     * Update offsets after domain length has changed.
     */
    void _update_offsets();

    /*
     * Communicate atoms that have left the local domain to the neighboring
     * domains in Cartesian direction *dim*. Ghost buffers will be
     * invalidated after a call to this method.
     */
    Eigen::Index _exchange_atoms(Atoms &atoms, int dim);

    /*
     * Communicate atoms into the ghost buffers of neighboring cells in
     * Cartesian direction *dim*.
     */
    std::tuple<Eigen::Index, Eigen::Index>
    _update_ghosts(Atoms &atoms, double border_width, int dim,
                   Eigen::Index left_start, Eigen::Index left_len,
                   Eigen::Index right_start, Eigen::Index right_len);

    // MPI communicator
    MPI_Comm comm_;

    // Length of the cell
    Eigen::Array3d domain_length_;

    // Cartesian decomposition of grid
    const Eigen::Array3i decomposition_;

    // Periodicity
    const Eigen::Array3i periodicity_;

    // Present size and rank
    int size_, rank_;

    // Cartesian coordinate of the present rank
    Eigen::Array3i coordinate_;

    // Ranks of process to the left and to the right
    Eigen::Array3i left_, right_;

    // Number of process-local atoms
    int nb_local_;

    // Is the domain decomposition enabled?
    bool is_enabled_;

    // Offsets for periodic boundary conditions
    Eigen::Matrix3d offset_left_, offset_right_;
};


#endif //YAMD_DOMAIN_H
