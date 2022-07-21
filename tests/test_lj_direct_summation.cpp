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

#include "gtest/gtest.h"

#include "../include/Molecular-Dynamics/lj_direct_summation.h"

TEST(LJDirectSummationTest, Forces) {
  constexpr int nb_atoms = 10;
  constexpr double epsilon =
      0.7; // choose different to 1 to pick up missing factors
  constexpr double sigma = 0.3;
  constexpr double delta = 0.0001; // difference used for numerical (finite
                                   // difference) computation of forces
  Atoms atoms{nb_atoms};
  atoms.positions.setRandom(); // random numbers between -1 and 1

  // compute and store energy of the indisturbed configuration
  double e0{lj_direct_summation(atoms, epsilon, sigma)};
  Forces_t forces0{atoms.forces};

  // loop over all atoms and compute forces from a finite differences
  // approximation
  for (int i{0}; i < nb_atoms; ++i) {
    // loop over all Cartesian directions
    for (int j{0}; j < 3; ++j) {
      // move atom to the right
      atoms.positions(j, i) += delta;
      double eplus{lj_direct_summation(atoms, epsilon, sigma)};
      // move atom to the left
      atoms.positions(j, i) -= 2 * delta;
      double eminus{lj_direct_summation(atoms, epsilon, sigma)};
      // move atom back to original position
      atoms.positions(j, i) += delta;

      // finite-differences forces
      double fd_force{-(eplus - eminus) / (2 * delta)};

      // check whether finite-difference and analytic forces agree
      if (abs(forces0(j, i)) > 1e-10) {
        EXPECT_NEAR(abs(fd_force - forces0(j, i)) / forces0(j, i), 0, 1e-5);
      } else {
        EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
      }
    }
  }
}