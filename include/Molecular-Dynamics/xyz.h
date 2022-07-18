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

#ifndef YAMD_XYZ_H
#define YAMD_XYZ_H

#include <fstream>

#include "atoms.h"

/*
 * Type Names_t, if not defined use
 * using Names_t = std::vector<std::string>;
 */

/*
 * Read positions from an XYZ file. XYZ files are structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z
 *         where Name is some name for the atom and X Y Z the position
 */
std::tuple<Names_t, Positions_t> read_xyz(const std::string &filename);

/*
 * Read positions and velocities from an XYZ file.
 * The XYZ file is structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z VX VY VZ
 *         where Name is some name for the atom, X Y Z the position
 *         and VX, VY, VZ the velocity of the atom
 */
std::tuple<Names_t, Positions_t, Velocities_t> read_xyz_with_velocities(const std::string &filename);

/*
 * Write positions and velocities to an XYZ file.
 * The XYZ file is structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z VX VY VZ
 *         where Name is some name for the atom, X Y Z the position
 *         and VX, VY, VZ the velocity of the atom
 */
void write_xyz(std::ofstream &file, Atoms& atoms);

/*
 * Write positions and velocities to an XYZ file.
 * The XYZ file is structured a follows:
 *     line 1: Number of atoms
 *     line 2: Comment line (is ignored)
 *     following lines: Name X Y Z VX VY VZ
 *         where Name is some name for the atom, X Y Z the position
 *         and VX, VY, VZ the velocity of the atom
 */
void write_xyz(const std::string &filename, Atoms& atoms);

#endif //YAMD_XYZ_H