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

#include <iostream>
#include <numeric>

#include "../include/Molecular-Dynamics/neighbors.h"

NeighborList::NeighborList(double interaction_range) :
        interaction_range_{interaction_range}, seed_{1}, neighbors_{1} {}

const std::tuple<const Eigen::ArrayXi &, const Eigen::ArrayXi &> NeighborList::update(const Atoms &atoms) {
    // Shorthand for atoms.positions.
    auto &&r{atoms.positions};

    // Origin stores the bottom left corner of the enclosing rectangles and lengths the three Cartesian lengths.
    Eigen::Array3d origin{3}, lengths{3};

    // This is the number of cells/grid points that fit into the enclosing rectangle. The grid is such that a sphere of
    // diameter interaction_range_ fits into each cell.
    Eigen::Array3i nb_grid_pts{3};

    // Compute box that encloses all atomic positions. Make sure that box lengths are exactly divisible by the
    // interaction range. Also compute the number of cells in each Cartesian direction.
    for (int i{0}; i < 3; ++i) {
        origin(i) = r.row(i).minCoeff();
        lengths(i) = r.row(i).maxCoeff() - origin(i);
        nb_grid_pts(i) = static_cast<int>(std::ceil(lengths(i) / interaction_range_));

        // This can only happen if all atoms sit within a plane. It is unlikely, but...
        if (nb_grid_pts(i) <= 0) {
            nb_grid_pts(i) = 1;
        }

        double padding_length{nb_grid_pts(i) * interaction_range_ - lengths(i)};
        origin(i) -= padding_length / 2;
        lengths(i) += padding_length;
    }

    // Compute cell indices. The follow array contains the cell index for each atom.
    Eigen::ArrayXi atom_to_cell{coordinate_to_index(
            ((r.row(0) - origin(0)) * nb_grid_pts(0) / lengths(0)).floor().cast<int>(),
            ((r.row(1) - origin(1)) * nb_grid_pts(1) / lengths(1)).floor().cast<int>(),
            ((r.row(2) - origin(2)) * nb_grid_pts(2) / lengths(2)).floor().cast<int>(),
            nb_grid_pts)};

    // We now sort the cell indices. This will allow us to search for the atoms that sit in neighboring cells.
    Eigen::ArrayXi sorted_atom_indices{atom_to_cell.size()};

    // Fill array `sorted_atom_indices` with consecutive numbers, i.e. 0, 1, 2, 3, 4, 5, ...
    std::iota(sorted_atom_indices.begin(), sorted_atom_indices.end(), 0);

    // Sort the array `sorted_atom_indices` by cell index. This yields an array of atom indices, but now sorted by
    // cell, i.e. all atoms within cell 0 are at the beginning of the array, followed by all atoms in cell 1 etc..
    // Example:
    //     sorted_atom_indices:                2 4 9 6 7 8 0 1 3 9
    //     atom_to_cell(sorted_atom_indices):  0 0 0 0 1 1 1 2 2 3
    //                                         ^       ^     ^   ^
    //     cell_index:                         0       1     2   3
    //     entry_index:                        0       4     7   9
    std::sort(sorted_atom_indices.begin(), sorted_atom_indices.end(),
              [&](int i, int j) { return atom_to_cell[i] < atom_to_cell[j]; });

    // We now build an array that points to the first entry within a certain cell in the `sorted_atom_indices` array.
    // We use a std::vector because we need to dynamically grow this array.
    std::vector<std::tuple<int, int>> binned_atoms{};
    int cell_index{atom_to_cell(sorted_atom_indices(0))};
    int entry_index{0};

    // This stores the index of the first entry for each cell.
    binned_atoms.push_back({cell_index, entry_index});

    // We now loop over the sorted atom indices and check when the cell index changes.
    for (int i{1}; i < sorted_atom_indices.size(); ++i) {
        if (atom_to_cell(sorted_atom_indices(i)) != cell_index) {
            cell_index = atom_to_cell(sorted_atom_indices(i));
            entry_index = i;
            binned_atoms.push_back({cell_index, entry_index});
        }
    }

    // We are now in a position to build a neighbor list in linear order. We are doing a bit of optimization here.
    // Since we have a dynamically growing list, we don't want to resize every time we add a neighbor. We are therefore
    // doubling the size when necessary and then resizing once (to a shorter array) when the list has been build.
    seed_.resize(atoms.nb_atoms() + 1);

    int n{0};
    auto interaction_range_sq{interaction_range_ * interaction_range_};
    for (int i{0}; i < atoms.nb_atoms(); ++i) {
        seed_(i) = n;

        Eigen::Array3i cell_coord{
                static_cast<int>(std::floor(nb_grid_pts(0) * (r(0, i) - origin(0)) / lengths(0))),
                static_cast<int>(std::floor(nb_grid_pts(1) * (r(1, i) - origin(1)) / lengths(1))),
                static_cast<int>(std::floor(nb_grid_pts(2) * (r(2, i) - origin(2)) / lengths(2)))};

        // Loop over neighboring cells.
        for (int x = -1; x <= 1; ++x) {
            for (int y = -1; y <= 1; ++y) {
                for (int z = -1; z <= 1; ++z) {
                    Eigen::Array3i neigh_cell_coord{cell_coord(0) + x, cell_coord(1) + y, cell_coord(2) + z};

                    if ((neigh_cell_coord >= 0).all() && (neigh_cell_coord < nb_grid_pts).all()) {
                        int cell_index{coordinate_to_index(neigh_cell_coord, nb_grid_pts)};

                        // Find first entry within the cell neighbor list.
                        auto cell{std::lower_bound(binned_atoms.begin(), binned_atoms.end(), cell_index,
                                                   [&](const auto &i, const auto &j) {
                                                       return std::get<0>(i) < j;
                                                   })};
                        if (cell != binned_atoms.end() && std::get<0>(*cell) == cell_index) {
                            for (int j{std::get<1>(*cell)};
                                 j < atom_to_cell.size() && atom_to_cell(sorted_atom_indices(j)) == cell_index; ++j) {
                                auto neighi{sorted_atom_indices(j)};

                                // Exclude the atom from being its own neighbor
                                if (neighi != i) {
                                    Eigen::Array3d distance_vector{r.col(i) - r.col(neighi)};
                                    double distance_sq{(distance_vector * distance_vector).sum()};
                                    if (distance_sq <= interaction_range_sq) {
                                        if (n >= neighbors_.size()) {
                                            neighbors_.conservativeResize(2 * neighbors_.size());
                                        }
                                        neighbors_(n) = neighi;
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    seed_(atoms.nb_atoms()) = n;
    neighbors_.conservativeResize(n);

    return {seed_, neighbors_};
}