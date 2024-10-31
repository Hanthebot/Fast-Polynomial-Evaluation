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

#include "gfelement.h"
#include "modular_arith.h"

namespace shk_galoiscpp
{


//-----------------------------------------------------------------------------
// Constructors and destructors
//-----------------------------------------------------------------------------

GFelement::GFelement(): field{ nullptr }, value{0} {}

GFelement::GFelement(GaloisField* gf): field{gf}, value{0} {}

//------------------------------------------------------
GFelement::GFelement(GaloisField* gf, const Fint val): field{gf}, value{val} {}


//------------------------------------------------------
GFelement::GFelement(const GFelement& gfe): field{gfe.field}, value{gfe.value} {}


//------------------------------------------------------
// GFelement::~GFelement()
// {
//    delete [] polynomial;
// }


//-----------------------------------------------------------------------------
// mutators
//-----------------------------------------------------------------------------

void GFelement::setZero(void)
{
   value = 0;
}



//------------------------------------------------------
void GFelement::setOne(void)
{
   value = 1;
}

//------------------------------------------------------


//-----------------------------------------------------------------------------
// operators
//-----------------------------------------------------------------------------

// Fint& GFelement::operator[](Int index)
// {
//    return polynomial[index];
// }



//------------------------------------------------------
// Fint GFelement::operator[](Int index) const
// {
//    return polynomial[index];
// }



//------------------------------------------------------
GFelement& GFelement::operator=(const GFelement& right)
{
   if (this != &right)
   {
      field = right.field;
      value = right.value;
   }
   return *this;
}



//------------------------------------------------------
GFelement GFelement::operator-()
{
   GFelement result(field);
   
   result.value = (value>0) ? (field->getModulus() - value) : (Fint) 0;

   return result;
}



//------------------------------------------------------
bool operator==(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   return left.value == right.value;
}



//------------------------------------------------------
bool operator!=(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   return left.value != right.value;
}



//------------------------------------------------------
bool operator<(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   return left.value < right.value;
}


//------------------------------------------------------
bool operator>(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   return left.value > right.value;
}



//------------------------------------------------------
bool operator<=(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   return ((left==right) | (left<right)) ? true : false;
}



//------------------------------------------------------
bool operator>=(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   return ((left==right) | (left>right)) ? true : false;
}


//------------------------------------------------------
bool operator<(const GFelement& left, const Fint& right)
{
   return left.getX() < right;
}


//------------------------------------------------------
bool operator>(const GFelement& left, const Fint& right)
{
   return left.getX() > right;
}



//------------------------------------------------------
bool operator<=(const GFelement& left, const Fint& right)
{
   return ((left.getX()==right) | (left<right)) ? true : false;
}



//------------------------------------------------------
bool operator>=(const GFelement& left, const Fint& right)
{
   return ((left.getX()==right) | (left>right)) ? true : false;
}


//------------------------------------------------------
bool operator<(const GFelement& left, const unsigned int& right)
{
   return left.getX() < right;
}


//------------------------------------------------------
bool operator>(const GFelement& left, const unsigned int& right)
{
   return left.getX() > right;
}



//------------------------------------------------------
bool operator<=(const GFelement& left, const unsigned int& right)
{
   return ((left.getX()==right) | (left<right)) ? true : false;
}



//------------------------------------------------------
bool operator>=(const GFelement& left, const unsigned int& right)
{
   return ((left.getX()==right) | (left>right)) ? true : false;
}



//------------------------------------------------------
GFelement operator+(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   GFelement result(left.field);
   result.value = left.value + right.value;
   mpz_mod(result.value.get_mpz_t(), result.value.get_mpz_t(), left.field->getModulus().get_mpz_t());

   return result;
}



//------------------------------------------------------
GFelement operator-(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   GFelement result(left.field);
   result.value = left.value - right.value;
   mpz_mod(result.value.get_mpz_t(), result.value.get_mpz_t(), left.field->getModulus().get_mpz_t());

   return result;
}



//------------------------------------------------------
GFelement operator*(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;

   GFelement result(left.field);
   result.value = left.value * right.value;
   result.value %= left.field->getModulus();

   return result;
}



//------------------------------------------------------
GFelement operator/(const GFelement& left, const GFelement& right)
{
   if (left.field != right.field) throw ErrorIncompatibleFields;
   
   return left * GFmultInverse( right );
}



//------------------------------------------------------
GFelement operator*(Fint left, const GFelement& right)
{
   GFelement result(right.field);
   result.value = left * right.value;
   result.value %= right.field->getModulus();

   return result;
}



//------------------------------------------------------
GFelement operator*(const GFelement& left, Fint right)
{

   GFelement result(left.field);
   result.value = left.value * right;
   result.value %= left.field->getModulus();

   return result;
}



//------------------------------------------------------
GFelement operator/(const GFelement& left, Fint right)
{
   GFelement result(left.field);

   return left * inverseModular(right, left.field->getModulus());
}



//------------------------------------------------------
ostream& operator<<(ostream& output, const GFelement& right)
{
   output << right.value; // << endl;

   return output;
}



//------------------------------------------------------
GFelement GFmultInverse(const GFelement& gfe)
{
   GFelement result(gfe.field);
   result.value = inverseModular(gfe.value, gfe.field->getModulus());

   return result;
}



//------------------------------------------------------
GFelement GFexp(const GFelement& G, Fint m)
{
   GFelement S( G.getField(), 1 ); // create S to be in the same field as G

   if (m == 0) return S;

   GFelement Gtemp( G.getField() ); // create G in the same field as g
   Gtemp = G; // set Gtemp(x) = G(x)

   Fint mask = 1;

   if ((m & mask) != 0) // edit: GMP doesn't support bool conversion naturally
   {
      S = G;
   }

   for (Int i=0; i < static_cast<signed>(8*sizeof(Fint)); i++)
   {
      Gtemp = Gtemp * Gtemp;
      mask <<= 1;
      if ((m & mask) != 0) // edit: GMP doesn't support bool conversion naturally
      {
         S = Gtemp * S;
      }
   }

   return S;
}

//------------------------------------------------------


} // namespace
