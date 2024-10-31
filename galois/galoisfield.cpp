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

#include <iostream>
using namespace std;

#include <cassert>
#include "galoisfield.h"

namespace shk_galoiscpp
{

//------------------------------------------------------
GaloisField::GaloisField(Fint mod, Int dim, Fint omega)
{
   assert (mod>1);
   assert (dim>0);
   assert (omega>0);

   modulus = mod;
   dimension = dim;
   rou = omega;
}


//------------------------------------------------------
ostream& operator<<(ostream& output, const GaloisField& gf)
{
   output << "\nPrime Modulus: " << gf.modulus;
   output << "\nDimension: " << gf.dimension;

   output << "\nRoot of Unity: " << gf.rou << endl;

   return output;
}


} // namespace shk_galoiscpp
