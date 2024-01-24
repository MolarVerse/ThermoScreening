/**
 * @file gibbs_dftb.c
 * @author Stefanie Kröll (S.Kroell@student.uibk.ac.at) and Josef M. Gallmetzer (Josef.Gallmetzer@student.uibk.ac.at)
 * @brief Thermochemistry program to calculate thermodynamically properties like Gibbs Free energy
 * @version 0.1
 * @date 2022-06-10
 *
 * @copyright Copyright (c) 2022 
 */

 /**
 * institute: 
 * Department of General, Inorganic and Theoretical Chemistry, University of Innsbruck (6020)
 */

/**
 * group: 
 * T. Hofer (t.hofer@uibk.ac.at)
 */

 /**
 * license:
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/***  Define Headers ***/

#include <stdlib.h>    // macros: EXIT_SUCCESS, EXIT_FAILURE
#include <stdio.h>     // function: printf, scanf
#include <string.h>    // function:
#include <strings.h>   // function:
#include <errno.h>     // function:
#include <math.h>      // function:
#include "constants.h" // external: Definition of physical constants

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


/*** Define Temperature and Pressure of the system ***/

const double T = 298.15; // in unit K
const double P = 101325; // in unit Pa

/*** Define Prototype Functions for the external programms symmetry_i.c, getsymnum.f90 and atom_properties.c ***/

/**
 * symmetry_i.c (https://github.com/nquesada/symmetry): Taken from the xtb-source code
 * (https://github.com/grimme-lab/xtb/tree/main/symmetry
 * Explaination:
 * Brute force symmetry analyzer.
 * This is actually C++ program, masquerading as a C one!
 *
 * (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca
 * modifications by S. Dohm and S. Ehlert
*  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * $Log: symmetry.c,v $
 * Revision 1.16  2003/04/04  13:05:03  patchkov
 * Revision 1.15  2000/01/25  16:47:17  patchkov
 * Revision 1.14  2000/01/25  16:39:08  patchkov
 * Revision 1.13  1996/05/24  12:32:08  ps
 * Revision 1.12  1996/05/23  16:10:47  ps
 * First reasonably stable version.
 */

/**
 * getsymnum.f90: Taken from the xtb-source code
 * (https://github.com/grimme-lab/xtb/tree/main/symmetry)
 * Explaination: Modification of symmetry_i.c
 * Copyright (C) 2017-2020 Stefan Grimme
 * xtb is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * xtb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU Lesser General Public License for more details.
 You should have received a copy of the GNU Lesser General Public License along with xtb.  If not, see <https://www.gnu.org/licenses/>.
*/


/**
 * atom_properties.c written by T. Hofer
 * and generous made avaiable
 * Explaination: contains the atom masses in Molar Mass. 
 * Masses taken from www.webelements.com
 */


/**
 * @brief Compute the symmetry of the system in Schoenflies notation by the external program symmetry_i.c
 * Use: int getsymnum()
 * 
 * @param natoms number of atoms
 * @param attype type of atom from data type int
 * @param coord atom coordinates in Angstrom unit
 * @param symbol Schoenfliess notation of the point group
 * @param paramar list of parameter set for calculation
 */
void schoenflies(int natoms, int *attype, double *coord, char *symbol, double *paramar);



/**
 * @brief  Compute symmetry number by using symmetry_i.c and getsymnum.f90
 * 
 * @param pgroup point group in Schoenflies notation computed by void schoenflies()
 * @param linear if linear molecule: 1 (true) else 0 (false)
 * @return int symmetry number
 */
int getsymnum(char *pgroup, int linear);




/**
 * @brief Get the atom mass in Molar Mass unit by using atom_properties.c
 * 
 * @param name atom type
 * @return int mass
 */

int periodic_number(char name[]);
/**
 * @brief 
 * 
 * @param name 
 * @return double 
 */
double init_mass(char name[]);




/*** Declaration of Symmetry Number and Spin variable ***/

int sigmaR; // symmetry number for rotational contribution (calculated)
int s;      // spin (user input)



/*** Declaration of Functions ***/

/**
 * @brief Compute symmetry number sigma
 * Use: in rotational contribution
 *
 * @param natoms number of atoms
 * @param coords atom coordinates in Angstrom unit
 * @param atom_types type of atom from data type int
 * @param linear lineare molecule ? 1 : 0
 * @return int symmetry number
 */
int initialize_sigma(int natoms, double *coords, int *atom_types, int linear)
{ // Declaration and initialization of Schoenfliess symbol
  char symbol[100] = {0};
  // Declaration of a parameter set array for the calculations setting
  double paramar[11];

  // Initialization of Schoenflies calculation settings
  /** The same settings like in xtb source-code:
   * https://github.com/grimme-lab/xtb/blob/b9dbd81ecc5158898a213ee80afcc0b39a7efee9/src/thermo.f90
   */
  paramar[0] = -1;     // Verbosity
  paramar[1] = 10;     // MaxAxisOrder
  paramar[2] = 100;    // MaxOptCycles
  paramar[3] = 1e-3;   // ToleranceSame
  paramar[4] = 5e-1;   // TolerancePrimary
  paramar[5] = 1e-3;   // ToleranceFinal
  paramar[6] = 5e-1;   // MaxOptStep
  paramar[7] = 1e-7;   // MinOptStep
  paramar[8] = 1e-7;   // GradientStep
  paramar[9] = 1e-8;   // OptChangeThreshold
  paramar[10] = 1e-10; // OptChangeThreshold

  // Compute Schoenflies point group
  schoenflies(natoms, atom_types, coords, symbol, paramar);

  // Get symmetry number
  int symnum = getsymnum(symbol, linear);

  // return symmetry number used in the rotational contribution
  return symnum;
}

/**
 * @brief Calculation of the Mass center position for the moment of inertia calculations later.
 * Use: It is used for the rotational component. 
 * Formula: r_c = {Sum_{i=0}^{numAtom} r_{i}S * mass_{i} } / M , r={x,y,z} 
 * Explaination: Sum over all products of atom coordinates and atom masses and then divided by system mass
 *
 * @param mass mass of the atom in Molar Mass unit
 * @param xAxis x-coordinate of atom geometry in Angstrom unit
 * @param yAxis y-coordinate of atom geometry in Angstrom unit
 * @param zAxis z-coordinate of atom geometry in Angstrom unit
 * @param massCenter 3D-array to return the calculated mass center in Angstrom unit for x-,y-, z-direction
 * @param numAtoms number of atoms in system
 */

void massCenterPosition(double *mass, double *xAxis, double *yAxis, double *zAxis, double massCenter[3], int numAtoms)
{ 
  /** Declaration and initialization of total system mass.
   *  Sum over all atom masses in Molar Mass unit.
   */
  double M = 0;

  // Mass center calculations: Sum over all atoms

  for (size_t i = 0; i < (size_t)numAtoms; i++) 
  {
    // Partially calculation of mass center
    massCenter[0] += mass[i] * xAxis[i]; // x-coordinate mass center component
    massCenter[1] += mass[i] * yAxis[i]; // y-coordinate mass center component
    massCenter[2] += mass[i] * zAxis[i]; // z-coordinate mass center component

    // Sum over all atom masses to calculate total system mass
    M += mass[i]; 
  }
  // Division of the mass center by the total system mass
  massCenter[0] /= M; 
  massCenter[1] /= M;
  massCenter[2] /= M;
}

/**
 * @brief Transform the atom coordinates in center of mass system (CMS). 
 * Use: It is used for the rotational component. 
 * Formula: r_{CMS} = r - r_c 
 * Explaination: Atom coordinates are subtracted by the mass center.
 *
 * @param xAxis x-coordinate of atom geometry in Angstrom unit
 * @param yAxis y-coordinate of atom geometry in Angstrom unit
 * @param zAxis z-coordinate of atom geometry in Angstrom unit
 * @param massCenter mass center in Angstrom unit computed by void massCenterPosition()
 * @param numAtoms number of atoms in system
 */

void relocate(double *xAxis, double *yAxis, double *zAxis, double massCenter[3], int numAtoms)
{
  // Sum over all atoms
  for (size_t i = 0; i < (size_t)numAtoms; i++)
  {
    xAxis[i] -= massCenter[0]; // transformed x-coordinates
    yAxis[i] -= massCenter[1]; // transformed y-coordinates
    zAxis[i] -= massCenter[2]; // transformed z-coordinates
  }
}
/**
 * @brief Computation of moment of inertia tensor diagonal elements I_{ii}
 * Use: in void inertiaTensor() for rotational contribution 
 * Formula: I_{ii} = Sum_{i=0}^{numAtoms} mass_{i}*(iAxis^2 + jAxis_i^2), i is runtime variable over all atoms
 * @param iAxis parameter which coordinate (x,y,z)
 * @param jAxis parameter which coordinate (x,y,z)
 * @param mass atom mass in Molar mass units
 * @param numAtoms number of atoms
 * @return double
 */
double inertiaFunctionII(double *iAxis, double *jAxis, double *mass, int numAtoms)
{
  // Declaration and initialization diagonal element the moment of inertia tensor
  double sum = 0;
  // Sum over all atoms
  for (size_t i = 0; i < (size_t)numAtoms; i++)
  {
    sum += mass[i] * (iAxis[i] * iAxis[i] + jAxis[i] * jAxis[i]);
  }
  // return the diagonal element I_{ii}
  return sum;
}
/**
 * @brief Computation of moment of inertia tensor non-diagonal elements I_{ij}
 * Use: in void inertiaTensor() for rotational contribution 
 * Formula: I_{ij} = -Sum_{i=0}^{numAtoms} mass_{i}* iAxis_i * jAxis_i, i is runtime variable over all atoms
 * @param iAxis parameter which coordinate (x,y,z)
 * @param jAxis parameter which coordinate (x,y,z)
 * @param mass atom mass in Molar mass units
 * @param numAtoms number of atoms
 * @return double
 */
double inertiaFunctionIJ(double *iAxis, double *jAxis, double *mass, int numAtoms)
{
  // Declaration and initialization non-diagonal element the moment of inertia tensor
  double sum = 0;
  // Sum over all atoms
  for (size_t i = 0; i < (size_t)numAtoms; i++)
  {
    sum -= mass[i] * iAxis[i] * jAxis[i];
  }
  // return non-diagonal element
  return sum;
}
/**
 * @brief Compute inertia tensor
 * Use: void eigenValues() for rotational contribution
 *
 * @param xAxis x-coordinates in Angstrom
 * @param yAxis y-coordinates in Angstrom
 * @param zAxis z-coordinates in Angstrom
 * @param mass atom mass in Molar mass unit
 * @param numAtoms number of atoms
 * @param tensor parameter of computed tensor
 */
void inertiaTensor(double *xAxis, double *yAxis, double *zAxis, double *mass, int numAtoms, double tensor[3][3])
{
  // Computation of all tensor elements by using void inertiaFunctionII() and void inertiaFunctionIJ()
  tensor[0][0] = inertiaFunctionII(yAxis, zAxis, mass, numAtoms);
  tensor[0][1] = inertiaFunctionIJ(xAxis, yAxis, mass, numAtoms);
  tensor[1][0] = tensor[0][1];
  tensor[0][2] = inertiaFunctionIJ(xAxis, zAxis, mass, numAtoms);
  tensor[2][0] = tensor[0][2];
  tensor[1][1] = inertiaFunctionII(xAxis, zAxis, mass, numAtoms);
  tensor[1][2] = inertiaFunctionIJ(yAxis, zAxis, mass, numAtoms);
  tensor[2][1] = tensor[1][2];
  tensor[2][2] = inertiaFunctionII(xAxis, yAxis, mass, numAtoms);
}

/**
 * @brief Compute eigenvalues
 * Use: void rot_contribution()
 * @param tensor inertia tensor computed by void inertiaTensor()
 * @param eigVals parameter of computed eigenvalues
 */
void eigenValues_old(double tensor[3][3], double eigVals[3])
{
  // store the tensor elements in different variables for better readability 
  /**
   *        [a b c]
   * T =    [d e f]
   *        [g h i]
   */

  double a = tensor[0][0];
  double b = tensor[0][1];
  double c = tensor[0][2];
  double d = tensor[1][0];
  double e = tensor[1][1];
  double f = tensor[1][2];
  double g = tensor[2][0];
  double h = tensor[2][1];
  double i = tensor[2][2];

  /**
   * Eigenvalue equation det(d*I-T) = 0
   * Transform into cubic equation:
   * A*x^3 + B*x^2 + C*x + D = 0
   */


  // Cubic coefficient: minus sign is pulled out
  double A = -1.0;

  // Quadratic coefficient: Sum of diagonal elements
  double B = i + e + a;

  // Linear coefficient: 
  double C = -e * i - a * i - a * e + g * c + b * d - h * f;

  // Constant Coefficient: 
  double D = a * e * i + c * d * h + b * f * g - g * c * e - d * b * i - h * f * a;

  /**
   * Solve Cubic equations: see therefore
   * https://stackoverflow.com/questions/13328676/c-solving-cubic-equations
   *
  */

  /**
   * Cubic discrimant: D = B^2*C^2-4*A*C^3-4*B^3*D+18*ABCD-27*A^2*D^2
   * D = q^3+r^2 with q = (3*C - B^2)/9, r = (-(27*D) + B * (9*C - 2*B^2))) / 54
   * 
   */

  // Transform to normal form x^3 + B*x^2 + C*x^ + D = 0

  B /= A;
  C /= A;
  D /= A;

  /** 
   * The cubic equation can be transformed into a depressed form:
   * z^3 + p*z + q = 0, with p = (3*C - B^2 / 9 and q = (-(27*D) + B*(9*C - 2*B^2)) / 54)
   * The discriminant is then D = p^3 + q^2
   */

  double p = (3.0 * C - (B * B)) / 9.0;
  double q = (-(27.0 * D) + B * (9.0 * C - 2.0 * (B * B))) / 54.0;
  double disc = p * p * p + q * q;

  double term1 = B / 3.0;

  if (disc > 0)
  {
    double s = q + sqrt(disc);
    s = (s < 0) ? -pow(-s, (1.0 / 3.0)) : pow(s, 1.0 / 3.0);
    double t = q - sqrt(disc);
    t = (t < 0) ? -pow(-t, (1.0 / 3.0)) : pow(t, 1.0 / 3.0);
    eigVals[0] = -term1 + s + t;
    term1 += (s + t) / 2.0;
    eigVals[1] = -term1;
    eigVals[2] = -term1;
    return;
  }

  double q13;

  if (disc == 0)
  {
    q13 = (q < 0) ? -pow(-q, (1.0 / 3.0)) : pow(q, (1.0 / 3.0));
    eigVals[0] = -term1 + 2.0 * q13;
    eigVals[1] = -(q13 + term1);
    eigVals[2] = -(q13 + term1);
    return;
  }

  p *= -1.0;
  double dum1 = p * p * p;
  dum1 = acos(q / sqrt(dum1));
  q13 = 2.0 * sqrt(p);
  eigVals[0] = -term1 + q13 * cos(dum1 / 3.0);
  eigVals[1] = -term1 + q13 * cos((dum1 + 2.0 * M_PI) / 3.0);
  eigVals[2] = -term1 + q13 * cos((dum1 + 4.0 * M_PI) / 3.0);

  printf("Eigenvalues old in Angstrom: %.12lf %.12lf %.12lf\n\n", eigVals[0], eigVals[1], eigVals[2]);

  return;
}
void eigenValues_new(double tensor[3][3], double eigVals[3])
{
  // store the tensor elements in different variables for better readability 
  /**
   *        [a b c]
   * T =    [d e f]
   *        [g h i]
   */

  double a = tensor[0][0];
  double b = tensor[0][1];
  double c = tensor[0][2];
  double d = tensor[1][0];
  double e = tensor[1][1];
  double f = tensor[1][2];
  double g = tensor[2][0];
  double h = tensor[2][1];
  double i = tensor[2][2];

  /**
   * Eigenvalue equation det(d*I-T) = 0
   * Transform into cubic equation:
   * A*x^3 + B*x^2 + C*x + D = 0
   */


  // Cubic coefficient: minus sign is pulled out
  double A = -1.0;

  // Quadratic coefficient: Sum of diagonal elements
  double B = i + e + a;

  // Linear coefficient: 
  double C = -e * i - a * i - a * e + g * c + b * d + h * f;

  // Constant Coefficient: 
  double D = a * e * i + c * d * h + b * f * g - g * c * e - d * b * i - h * f * a;

  /**
   * Solve Cubic equations: see therefore
   * https://stackoverflow.com/questions/13328676/c-solving-cubic-equations
   *
  */

  /**
   * Cubic discrimant: D = B^2*C^2-4*A*C^3-4*B^3*D+18*ABCD-27*A^2*D^2
   * D = q^3+r^2 with q = (3*C - B^2)/9, r = (-(27*D) + B * (9*C - 2*B^2))) / 54
   * 
   */

  // Transform to normal form x^3 + B*x^2 + C*x^ + D = 0

  B /= A;
  C /= A;
  D /= A;

  /** 
   * The cubic equation can be transformed into a depressed form:
   * z^3 + p*z + q = 0, with p = (3*C - B^2 / 9 and q = (-(27*D) + B*(9*C - 2*B^2)) / 54)
   * The discriminant is then D = p^3 + q^2
   */

  double p = (3.0 * C - (B * B)) / 9.0;
  double q = (-(27.0 * D) + B * (9.0 * C - 2.0 * (B * B))) / 54.0;
  double disc = p * p * p + q * q;

  double term1 = B / 3.0;

  if (disc > 0)
  {
    double s = q + sqrt(disc);
    s = (s < 0) ? -pow(-s, (1.0 / 3.0)) : pow(s, 1.0 / 3.0);
    double t = q - sqrt(disc);
    t = (t < 0) ? -pow(-t, (1.0 / 3.0)) : pow(t, 1.0 / 3.0);
    eigVals[0] = -term1 + s + t;
    term1 += (s + t) / 2.0;
    eigVals[1] = -term1;
    eigVals[2] = -term1;
    return;
  }

  double q13;

  if (disc == 0)
  {
    q13 = (q < 0) ? -pow(-q, (1.0 / 3.0)) : pow(q, (1.0 / 3.0));
    eigVals[0] = -term1 + 2.0 * q13;
    eigVals[1] = -(q13 + term1);
    eigVals[2] = -(q13 + term1);
    return;
  }

  p *= -1.0;
  double dum1 = p * p * p;
  dum1 = acos(q / sqrt(dum1));
  q13 = 2.0 * sqrt(p);
  eigVals[0] = -term1 + q13 * cos(dum1 / 3.0);
  eigVals[1] = -term1 + q13 * cos((dum1 + 2.0 * M_PI) / 3.0);
  eigVals[2] = -term1 + q13 * cos((dum1 + 4.0 * M_PI) / 3.0);

  printf("Eigenvalues new in Angstrom: %.12lf %.12lf %.12lf\n\n", eigVals[0], eigVals[1], eigVals[2]);

  return;
}

/**
 * @brief Rotational contribution
 *
 * @param eigVal eigenvalues computed by void eigenValues()
 * @param Sr parameter of computed entropy contribution
 * @param Er parameter of computed energy contribution
 * @param Cr parameter of computed heat capacity contribution
 */

void rot_contribution(double eigVal[3], double *Sr, double *Er, double *Cr)
{
  // Declaration of variables:
  double eigValSI[3];       // eigenvalues in SI-units in x-,y-,z-direction
  double rotTemp[3];        // rotational temperature in x-,y-,z-direction
  double rotConst[3];       // rotational constants in x-,y-,z-direction
  double rotTempXYZ;        // total rotational temperature
  double qr;                // partitional function contribution

  // Computation of eigenvalue in SI-units, rotational temperature and rotational constants for every direction
  for (size_t i = 0; i < 3; i++)
  {
    eigValSI[i] = eigVal[i] * u * pow(10, -20);                                     // in unit kg m^2
    rotTemp[i] = (pow(h, 2) / (8.0 * pow(M_PI, 2) * kB)) / eigValSI[i];             // in unit Kelvin
    rotConst[i] = ((pow((h / (2 * M_PI)), 2) / 2) / eigValSI[i]) / h * pow(10, -9); // in unit GHz
  }
  // total rotational temperature
  rotTempXYZ = rotTemp[0] * rotTemp[1] * rotTemp[2];

  // check
  printf("rotTempXYZ in K^3: %.12lf\n\n", rotTempXYZ);

  // partitional function contribution
  qr = (pow(M_PI, 1.0 / 2) / sigmaR) * (pow(T, 3.0 / 2) / (pow(rotTempXYZ, 1.0 / 2)));
  // entropy contribution
  *Sr = R * (log(qr) + 3.0 / 2) / cal; // in unit cal/ Mol K
  // energy contribution
  *Er = ((3.0 / 2) * R * T) / H / Na;  // in unit H/Particle
  // heat capacity contribution
  *Cr = (3.0 / 2.0) * R / cal;         // in unit cal/ Mol K

  // check
  printf("rotation temperature in K: %.12lf %.12lf %.12lf\n\n", rotTemp[0], rotTemp[1], rotTemp[2]);
  printf("rotation constants in GHz: %.12lf %.12lf %.12lf\n\n", rotConst[0], rotConst[1], rotConst[2]);
  printf("qr: %.12lf, Sr: %.12lf cal/ mol K, Er:  %.12lf H/Particle, Cr: %.12lf cal/ mol K\n\n", qr, *Sr, *Er, *Cr);

  return;
}

/**
 * @brief Vibrational contribution
 *
 * @param modes computed vibrational mode provided by an input file
 * @param n degree of freedom
 * @param Sv parameter of computed entropy contribution
 * @param Ev parameter of computed energy contribution
 * @param Cv parameter of computed heat capacity contribution
 * @param EZP
 */

void vib_contribution(double modes[], const int n, double *Sv, double *Ev, double *Cv, double *EZP)
{

  double vibTemp[n];
  double qvi[n];
  double qv = 1;
  double Svi[n];
  double Evi[n];
  double Cvi[n];
  // Check if modes are not imaganary
  for (int i = 0; i < n; i++) 
  {
    if (modes[i] <= 0)
    {
      return EXIT_FAILURE;
    }
  }

  for (int i = 0; i < n; i++)
  {
    vibTemp[i] = h * modes[i] * c * pow(10, 2) / kB; // in unit K
  }

  for (int i = 0; i < n; i++)
  {
    // qvi[i] = 1/ (1 - exp(-vibTemp[i] / T));//BOT //
    qvi[i] = exp(-vibTemp[i] / (2 * T)) / (1 - exp(-vibTemp[i] / T)); // V0 //
    Svi[i] = ((vibTemp[i] / T) / (exp(vibTemp[i] / T) - 1)) - log(1 - exp(-vibTemp[i] / T));
    Evi[i] = (vibTemp[i] * (1.0 / 2 + 1 / (exp(vibTemp[i] / T) - 1)));
    Cvi[i] = exp(-vibTemp[i] / T) * pow(((vibTemp[i] / T) / (exp(-vibTemp[i] / T) - 1)), 2);
  }

  for (int i = 0; i < n; i++)
  {
    qv = qv * qvi[i];
    *Sv = (*Sv + Svi[i]);             // in unit cal / mol K
    *Ev = (*Ev + Evi[i]);             // in unit H/Mol
    *Cv = (*Cv + Cvi[i]);             // in unit cal/mol K
    *EZP = (*EZP + (vibTemp[i] / 2)); // in cal
  }
  *Sv = R * *Sv / cal;                 // in unit cal/mol K
  *Ev = R * *Ev / H / Na;              // in unit H/Particle
  *Cv = R * *Cv / cal;                 // in unit cal/mol K
  *EZP = R * *EZP / H / Na;            // in  H/Particle

  
  /*for (int i = 0; i < n; i++)
  {
  printf("vibration temperature in K: %.5e\n", vibTemp[i]);
  }*/

  // Check
  printf("qv: %.12lf, Sv: %.12lf cal/ mol K, Ev:  %.12lf H/Particle, Cv: %.12lf cal/ mol K, EZP: %.12lf\n\n", qv, *Sv, *Ev, *Cv, *EZP);

  return;
}

/**
 * @brief Translation contribution
 *
 * @param mSum
 * @param St parameter of computed entropy contribution
 * @param Et parameter of computed energy contribution
 * @param Ct parameter of computed heat capacity contribution
 */

void transl_contribution(double mSum, double *St, double *Et, double *Ct)
{

  // Declaration partitional function contribution
  double qt;

  // Convert mass into kg unit
  mSum = mSum * u; // in unit kg

  // partitional function contribution
  qt = pow((2.0 * M_PI * mSum * kB * T) / (pow(h, 2.0)), (3.0 / 2)) * (kB * T / P);
  // energy contribution
  *Et = (3.0 / 2 * R * T) / H / Na;    // in unit H/Particle
  // entropy contribution
  *St = R * (log(qt) + 5.0 / 2) / cal; // in unit cal/mol K
  // heat capacity contribution
  *Ct = 3.0 / 2 * R / cal;             // in unit cal/ mol K

  // Check
  printf("qt: %.12lf, St: %.12lf cal/ mol K, Et:  %.12lf H/Particle, Ct: %.12lf cal/ mol K\n\n", qt, *St, *Et, *Ct);

  return;
}

/**
 * @brief Electronic contribution
 *
 * @param Se parameter of computed entropy contribution
 * @param Ee parameter of computed energy contribution
 * @param Ce parameter of computed heat capacity contribution
 */

void elec_contribution(double *Se, double *Ee, double *Ce)
{
  // Declaration and intialization of spin multiplicity
  int omega0 = 2 * s + 1; // s = spin: user input

  // Declaration of partitional function contribution
  double qe;

  // partitional function contribution
  qe = omega0;
  // energy contribution
  *Ee = 0;
  // entropy contribution
  *Se = R * log(qe) / cal; // in unit cal/mol K
  // heat capacity contribution
  *Ce = 0;

  // Check
  printf("qe: %.12lf, Se: %.12lf cal/ mol K, Ee:  %.12lf H/Particle, Ce: %.12lf cal/ mol K\n\n", qe, *Se, *Ee, *Ce);
  return;
}


int main(int argc, char *argv[]) 
{

  // Output of temperature and pressure setting
  printf("T: %f K and P: %f Pa\n\n", T, P);

  /** Delcaration and initilization of variable symcal
   * It is used as tag to provide the user the ability 
   * to choose between symmetry calculation or without
   */

  int symcalc = 0; // default value 0
  char molecule[20];
  int charge = 0;
  /*** user input spin multiplicity and symcalc tag***/

  while (--argc > 0)
  {
    if (!strncmp(argv[argc], "-symcalc=", 9))
    {
      if (!strcasecmp(argv[argc] + 9, "true") || !strcasecmp(argv[argc] + 9, "t") || !strcasecmp(argv[argc] + 9, "1"))
      {
        symcalc = 1;
      }
    }
    else if (!strncmp(argv[argc], "-spinmul=", 9))
    {
      char *endptr = NULL;
      s = strtol(argc[argv] + 9, &endptr, 10);
      if (!endptr)
      {
        fprintf(stderr, "Could not read spin multiplicity\n\n");
      }
    }
    else if (!strncmp(argv[argc], "-molecule=", 10))
    {
      // char *endptr = NULL;
      strncpy(molecule,argv[argc]+10,strlen(argv[argc]));
      // char molecule[]=argc[argv]+10;
      printf("%s\n",molecule);
      // if (!endptr)
      // {
      //   fprintf(stderr, "Could not read molecule name\n\n");
      // }
    }
        else if (!strncmp(argv[argc], "-charge=", 8))
    {
      char *endptr = NULL;
      charge = strtol(argc[argv] + 8, &endptr, 9);
      if (!endptr)
      {
        fprintf(stderr, "Could not read charge\n\n");
      }
    }
  }
  
  // Output of choosen setting
  printf("Spin multiplicity: %i\n\n", s);
  printf("symcalc=%s\n\n", symcalc ? "true" : "false");

  /*** Reading input files ***/

  FILE *frequency = fopen("frequency.txt", "r");
  if (frequency == NULL)
  {
    fprintf(stderr, "\n(x) Frequency file (frequency.txt) is not available!\n\n");
    return EXIT_FAILURE;
  }

  FILE *geometry = fopen("geo_opt.xyz", "r");
  if (geometry == NULL)
  {
    fprintf(stderr, "\n(x) Geometry file (geo_opt.xyz) is not available!\n\n");
    return EXIT_FAILURE;
  }

  char line[256];

  // vibrational modes

  double *modes = malloc(sizeof(double));
  double garbage;
  int numModes = 0;

  while (fgets(line, sizeof(line), frequency) != NULL)
  {
    modes = realloc(modes, (numModes + 1) * sizeof(double));
    if (sscanf(line, "%lf %lf", &garbage, modes + numModes) != 2)
    {
      fprintf(stderr, "\n(x) Error in frequency.txt (line %d)\n\n", numModes + 1);
      return EXIT_FAILURE;
    }
    numModes++;
  }

  // Number of Atoms 

  int numAtoms;

  if (fscanf(geometry, "%d", &numAtoms) != 1)
  {
    fprintf(stderr, "\n(x) Number of Atoms in geo_opt.xyz is not available!\n\n");
    return EXIT_FAILURE;
  }

  fclose(geometry);

  printf("Number of atoms: %d\n\n", numAtoms);

  // Geometry

  geometry = fopen("geo_opt.xyz", "r");

  int *atom_types = malloc(numAtoms * sizeof(int));
  double *mass = malloc(numAtoms * sizeof(double));
  double *xAxis = malloc(numAtoms * sizeof(double));
  double *yAxis = malloc(numAtoms * sizeof(double));
  double *zAxis = malloc(numAtoms * sizeof(double));
  char atomName[8];

  int i = -2;
  double mSum = 0;

  while (fgets(line, sizeof(line), geometry) != NULL)
  {
    if (i < 0)
    {
      i++;
      continue;
    }

    if (sscanf(line, "%s %lf %lf %lf", atomName, xAxis + i, yAxis + i, zAxis + i) != 4)
    {
      fprintf(stderr, "\n(x) Error in geo_opt.xyz (line %d)\n\n", i + 3);
      return EXIT_FAILURE;
    }

    mass[i] = init_mass(atomName);                    // mass per atom in Molar Mass
    atom_types[i] = periodic_number(atomName);        // atom type
    mSum = mSum + mass[i];                            // total mass in Molar Mass
    printf("%s (%d): Masses: %f u, ", atomName, atom_types[i], mass[i]);
    i++;
  }
  printf("\n\n");
  printf("Sum of the Mass: %f u\n\n", mSum);

  /*** Checking if molecule is linear or not linear ***/

  /**
   * degree of freedom (dof)
   * dof = 3N-6 for non-linear molecules
   * dof = 3N-4 for linear molecules
   */

  int dof_l = 3 * numAtoms - 5; // dof for linear molecules
  int linear;
  /*
   * if number of modes is equal with linear dof
   * set the variable linear at 1 (true) else 0 (false)
   */
  if (numModes == dof_l) 
  {
    linear = 1;
  }
  else
  {
    linear = 0;
  }
  /*** Caclulating rotational symmetry number ***/

  // Allocate coordinaten variable coords
  double *coords = malloc(numAtoms * 3 * sizeof(double));
  // Initialize coords with coordinats
  for (int i = 0; i < numAtoms; i++)
  {
    coords[3 * i] = xAxis[i];
    coords[3 * i + 1] = yAxis[i];
    coords[3 * i + 2] = zAxis[i];
  }
  /*
   * if symcal tag is set at 1 (true) compute rotational symmetry 
   * else set rotational symmetry at 1
   */
  if (symcalc == 1)
  {
    sigmaR = initialize_sigma(numAtoms, coords, atom_types, linear);
  }
  else
  {
    sigmaR = 1;
  }

  // Freeing allocated variables
  free(coords);
  free(atom_types);

  // Print rotational symmetry number output
  printf("Rotational symmetry number: %d\n\n", sigmaR);

  // Closing files
  fclose(geometry);
  fclose(frequency);



  // Check: Geometry output
  printf("Geometry output:\n");
  printf("x axis coordinates in Angstrom:\n");
  for (int i = 0; i < numAtoms; i++)
  {
    printf("%.6e ", xAxis[i]);
  }
  printf("\n\n");
  printf("y axis coordinates in Angstrom:\n");
  for (int i = 0; i < numAtoms; i++)
  {
    printf("%.6e ", yAxis[i]);
  }
  printf("\n\n");
  printf("z axis coordinates in Angstrom:\n");
  for (int i = 0; i < numAtoms; i++)
  {
    printf("%.6e ", zAxis[i]);
  }
  printf("\n\n");

  /*** Calculate mass center ***/

  double massCenter[3] = {0};
  massCenterPosition(mass, xAxis, yAxis, zAxis, massCenter, numAtoms);

  /*** Relocate coordinates to MassCenter ***/

  relocate(xAxis, yAxis, zAxis, massCenter, numAtoms);

  // Check: Geometry output after relocation of the coordinates to mass center
  printf("Geometry relocated to the mass center:\n\n");
  printf("x axis coordinates in Angstrom:\n");
  for (int i = 0; i < numAtoms; i++)
  {
    printf("%.6e ", xAxis[i]);
  }
  printf("\n\n");
  printf("y axis coordinates in Angstrom:\n");
  for (int i = 0; i < numAtoms; i++)
  {
    printf("%.6e ", yAxis[i]);
  }
  printf("\n\n");
  printf("z axis coordinates in Angstrom:\n");
  for (int i = 0; i < numAtoms; i++)
  {
    printf("%.6e ", zAxis[i]);
  }
  printf("\n\n");

  /*** Calculate the inertia tensor ***/

  double tensor[3][3] = {{0}};
  inertiaTensor(xAxis, yAxis, zAxis, mass, numAtoms, tensor);

  // Check inertia tensor
  printf("inertia tensor in Angstrom:\n");
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      printf("%f ", tensor[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  /*** Calculate the eigenvalues ***/

  double eigVal_old[3] = {0};
  eigenValues_old(tensor, eigVal_old);

  double eigVal_new[3] = {0};
  eigenValues_new(tensor, eigVal_new);
  /*** Calculation of the rotation contribution ***/

  double Sr_old = 0;
  double Er_old = 0;
  double Cr_old = 0;

  double Sr_new = 0;
  double Er_new = 0;
  double Cr_new = 0;
  
  printf("Rot contribution old: \n");
  rot_contribution(eigVal_old, &Sr_old, &Er_old, &Cr_old);
  printf("Rot contribution new: \n");
  rot_contribution(eigVal_new, &Sr_new, &Er_new, &Cr_new);
  /*** Calculation of the vibration contribution ***/

  double Sv = 0;
  double Ev = 0;
  double Cv = 0;
  double EZP = 0;

  /**
   * pointer arithmetic, do not use the first 6 lines
   * to get the right number of vibrational modes
   * corresponding of the degree of freedom (dof)
   */
  if (linear == 1)
  {
    vib_contribution(modes + 5, numModes - 5, &Sv, &Ev, &Cv, &EZP);
  } 
  else 
  {
    vib_contribution(modes + 6, numModes - 6, &Sv, &Ev, &Cv, &EZP);
  }
  
   

  /*** Calculation of the translation contribution ***/

  double St = 0;
  double Et = 0;
  double Ct = 0;

  transl_contribution(mSum, &St, &Et, &Ct);

  /*** Calculation of the electronic contribution ***/

  double Se = 0;
  double Ee = 0;
  double Ce = 0;

  elec_contribution(&Se, &Ee, &Ce);

  /*** Calculating the total correction parts ***/

  double Etot_old;
  double Stot_old;
  double Ctot_old;

  Etot_old = (Er_old + Ev + Et + Ee); // in unit H/Particle
  Stot_old = (Sr_old + Sv + St + Se); // in unit cal/mol K
  Ctot_old = (Cr_old + Cv + Ct + Ce); // in unit cal/mol K

  double Etot_new;
  double Stot_new;
  double Ctot_new;

  Etot_new = (Er_new+ Ev + Et + Ee); // in unit H/Particle
  Stot_new = (Sr_new + Sv + St + Se); // in unit cal/mol K
  Ctot_new = (Cr_new + Cv + Ct + Ce); // in unit cal/mol K


  printf("Old:\n");
  printf("Etot: %.12lf H/Particle, Stot: %.12lf cal/ mol K,Ctot: %.12lf cal/ mol K\n\n", Etot_old, Stot_old, Ctot_old);
  printf("New:\n");
  printf("Etot: %.12lf H/Particle, Stot: %.12lf cal/ mol K,Ctot: %.12lf cal/ mol K\n\n", Etot_new, Stot_new, Ctot_new);

  /*** Calculating the entropy and gibbs free energy ***/

  double Hcorr_old;
  double Gcorr_old;


  double Hcorr_new;
  double Gcorr_new;

  Hcorr_old = Etot_old + (R * T) / H / Na; // in unit H/Particle

  Gcorr_old = Hcorr_old - (Stot_old * cal / Na / H) * T; // in unit H/Particle

  Hcorr_new = Etot_new + (R * T) / H / Na; // in unit H/Particle

  Gcorr_new = Hcorr_new - (Stot_new * cal / Na / H) * T; // in unit H/Particle


  // Check
  printf("Old:\n");
  printf("Hcorr: %.12lf in H/Particle, Gcorr: %.12lf in H/Particle\n\n", Hcorr_old, Gcorr_old);
  // Check
  printf("New:\n");
  printf("Hcorr: %.12lf in H/Particle, Gcorr: %.12lf in H/Particle\n\n", Hcorr_new, Gcorr_new);
   /*** Read external calculated total electronic energy from a file ***/

  FILE *electronic_energy = fopen("electronic_energy.txt", "r"); // open file
  if (electronic_energy == NULL)
  {
    // check if file exists
    fprintf(stderr, "\n(x) Electronic Energy file (electronic_energy.txt) is not available!\n\n");
    return EXIT_FAILURE;
  }

  double E0;
  fscanf(electronic_energy, "%lf", &E0); // read value

  fclose(electronic_energy); // close file

  // Check
  printf("Old:\n");
  printf("Electronic Energy E0: %lf H/Particle\n\n", E0);

  /*** Summarized thermochemical properties ***/

  // Declaration of variables
  double E0EZP;
  double E0Etot_old;
  double E0Hcorr_old;
  double E0Gcorr_old;

  double E0Etot_new;
  double E0Hcorr_new;
  double E0Gcorr_new;


  E0EZP = (E0 + EZP);           // Sum of electronic and zero-point energies
  E0Etot_old = (E0 + Etot_old);         // Sum of electronic and thermal energies
  E0Hcorr_old = (E0 + Hcorr_old);       // Sum of electronic and thermal enthalpies
  E0Gcorr_old = (E0 + Gcorr_old);       // Sum of electronic and thermal Free Energies


  E0Etot_new = (E0 + Etot_new);         // Sum of electronic and thermal energies
  E0Hcorr_new = (E0 + Hcorr_new);       // Sum of electronic and thermal enthalpies
  E0Gcorr_new = (E0 + Gcorr_new);       // Sum of electronic and thermal Free Energies
  /*** Summarized ouput ***/

  printf("Zero-Point Energy: %.12lf in H/Particle\n\n", EZP);
  printf("Sum of electronic and zero-point Energy: %.12lf in H/Particle\n\n", E0EZP);
  printf("Old\n");
  printf("Total Energy: %.12lf in H/Particle\n\n", E0Etot_old);
  printf("Total Enthalpy: %.12lf in H/Particle\n\n", E0Hcorr_old);
  printf("Total Free Energy: %.12lf in H/Particle\n\n", E0Gcorr_old);

  printf("New\n");
  printf("Total Energy: %.12lf in H/Particle\n\n", E0Etot_new);
  printf("Total Enthalpy: %.12lf in H/Particle\n\n", E0Hcorr_new);
  printf("Total Free Energy: %.12lf in H/Particle\n\n", E0Gcorr_new);
  // Freeing the allocated arrays.
  FILE *sF = fopen("/home/stk/masterthesis/ubuntu_stuff/masterarbeit/calculations/thermochemistry/DMF_COSMO/all_compounds/control_results_old_new.csv","a");
  if(sF==NULL){
    fprintf(stderr, "\n(x) File 'control_results_old_new.csv' for save is not available!\n\n");
    return EXIT_FAILURE;
  }
  fprintf(sF,"%s",molecule);
  fprintf(sF," ");
  fprintf(sF,"%d",charge);
  fprintf(sF," ");
  fprintf(sF,"%f",E0Gcorr_old);
  fprintf(sF," ");
  fprintf(sF,"%f",E0Gcorr_new);
  fprintf(sF," ");
  fprintf(sF,"%f",E0Gcorr_new-E0Gcorr_old);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_old[0]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_new[0]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_new[0]-eigVal_old[0]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_old[1]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_new[1]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_new[1]-eigVal_old[1]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_old[2]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_new[2]);
  fprintf(sF," ");
  fprintf(sF,"%f",eigVal_new[2]-eigVal_old[2]);
  fprintf(sF,"\n");
  fclose(sF);
  free(xAxis);
  free(yAxis);
  free(zAxis);
  free(mass);

  return EXIT_SUCCESS;
}
