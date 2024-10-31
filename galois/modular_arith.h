/**********************************************************************
   Project: C++ Library for General Galois Field Arithmetic

   Language: C++ 2007	   
   Author: Saied H. Khayat
   Date:   Feb 2013
   URL: https://github.com/saiedhk
   
   Copyright Notice: Free use of this library is permitted under the
   guidelines and in accordance with the MIT License (MIT).
   http://opensource.org/licenses/MIT

**********************************************************************/

/*
   Modular arithmetic operations
*/

#include <iostream>
using namespace std;

#include "typedefs.h"

#ifndef MODULAR_ARITH_H
#define MODULAR_ARITH_H

namespace shk_galoiscpp
{

/**
   Implements the Extended Euclid Algorithm.

   It computes the greatest common divisor d of two positive integers a and b
   and determine two integers x and y such that d = ax + by.

   (See Chapter 2, Algorithm 2.107 in Handbook of Applied Cryptography,
   by Menezes, van Oorschot, Vanstone, CRC Press, 1996.)

   @param a input integer (a > 0)
   @param b input integer (b > 0)
   @param x output integer
   @param y output integer
   @param d output integer (greatest common divisor of a and b)
*/
// void extendedEuclid(Fint a, Fint b, Fint& x, Fint& y, Fint& d);


//------------------------------------------------
// modular arithmetic functions
//------------------------------------------------

/**
   Determines the multiplicative inverse of an integer modulo mod
   @param a input integer (a>0)
   @param mod input integer (modulus, mod>0)
   @returns b the inverse of a such that a * b = 1 modulo mod
   @throws ErrorNoInverse if mod is not prime relative to a
*/
inline Fint inverseModular(Fint a, Fint mod)
{
   Fint result;
   int flag = mpz_invert(result.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());

   if (flag == 0)
   {
      throw ErrorNoInverse;
   }

   return result;
}

} // namespace shk_galoiscpp

#endif