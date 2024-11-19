#ifndef GFELEMENT_H
#define GFELEMENT_H

#include <iostream>
using namespace std;

#include "galoisfield.h"
// #include "polynomial_arith.h"


namespace shk_galoiscpp
{

/**
   Defines a Galois Field element. Also defines arithmetic operations on Galois Field elements.

   In mathematics, an element of the Galois Field GF(p^k) is a polynomial of degree less than k
   with coefficients chosen from the set {0,1,2,..., p-1}.

   In this program, a GFelement contains the following private data:
   @param field a pointer to the Galois Field this GFelement belongs to
   // @param modulus the prime modulus of the Galois Field
   // @param dimension the dimension of the Galois Field
   @param value the array representing the polynomial describing this GFelement
*/
class GFelement
{

   public:

      //---------------------------
      // constructors
      //---------------------------

      GFelement();
      /**
         Constructs GFelement without initializing its coefficients
         @param gf pointer to a Galois Field
      */
      GFelement(GaloisField* gf);


      /**
         Constructs a GFelement and initializes its coefficients
         @param gf pointer to a Galois Field
         @param poly an external array containing coefficients (length of array = dimension)
      */
      GFelement(GaloisField* gf, const Fint value);


      /**
         Constructs and clones a GFelement.
         @param gfe reference to a GFelement
      */
      GFelement(const GFelement& gfe);


      /**
         Destructs a GFelement and frees its allocated memory
      */
      // ~GFelement();


      //---------------------------
      // accessors
      //---------------------------

      /**
         Returns the pointer to the GF this GFelement is associated with.
      */
      GaloisField* getField() const;


      /**
         Returns the i-th coefficient of this GFelement's polynomial
      */
      Fint getX() const;


      //---------------------------
      // mutators
      //---------------------------

      /**
         Sets GFelement to the ZERO polynomial.
         It sets all coefficients of GFelement polynomial to zero.
      */
      void setZero(void);


      /**
         Sets GFelement to the ONE polynomial.
         It sets all coefficients of GFelement polynomial to zero except the zero-th
         degree coefficient which is set to 1.
      */
      void setOne(void);


      /**
         Sets the i-th coefficient of GFelement's polynomial to a specified constant
         @param X the value of the value
      */
      void setX(Fint X);
      void setField(GaloisField* field);

      GFelement getRootOfUnity() const;
      GFelement getZero() const;
      GFelement getOne() const;

      // /**
      //    Sets all coefficients of GFelement's polynomial to the same specified constant
      //    @param value input, the value for all coefficients
      // */
      // void setCoefficients(Fint value);


      //---------------------------
      // operators
      //---------------------------

      // /**
      //    Allows accessing the i-th coefficient of GFelement polynomial.

      //    This can be used as an l-value. It sets the i-th coefficient to a specified value.
      //    @param i (input), index
      //    @returns Reference of the i-th coefficient of GFelement polynomial
      // */
      // Fint& operator[](Int i);


      // /**
      //    Returns the i-th coefficient of GFelement polynomial
      //    @param i (input), index
      //    @returns the i-th coefficient of GFelement polynomial
      // */
      // Fint operator[](Int i) const;


      /**
         Assigns/copies the values of another GFelement to this GFelement
         @param right GFelement on the right of = sign
         @returns GFelement on the left of = sign (this GFelement)
      */
      GFelement& operator=(const GFelement& right);

      GFelement& operator=(const Fint& right);


      /**
         Returns the negated copy of this GFelement
      */
      GFelement operator-();

      /**
         Increments
      */
      GFelement& operator++();


      /**
         Compares two GFelements, checks for equality
         @param right GFelement on the right of == sign
         @param left GFelement on the left of == sign
         @returns true if equal, false otherwise
      */
      friend bool operator==(const GFelement& left, const GFelement& right);


      /**
         Compares two GFelements, checks for inequality
         @param right GFelement on the right of != sign
         @param left GFelement on the left of != sign
         @returns true if equal, false otherwise
      */
      friend bool operator!=(const GFelement& left, const GFelement& right);


      /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of < sign
         @param left GFelement on the left of < sign
         @returns true if left < right, false otherwise
      */
      friend bool operator<(const GFelement& left, const GFelement& right);


      /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of > sign
         @param left GFelement on the left of > sign
         @returns true if left > right, false otherwise
      */
      friend bool operator>(const GFelement& left, const GFelement& right);


      /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of <= sign
         @param left GFelement on the left of <= sign
         @returns true if left <= right, false otherwise
      */
      friend bool operator<=(const GFelement& left, const GFelement& right);


      /**
         Compares two GFelements based on their lexicographical ordering.
         @param right GFelement on the right of >= sign
         @param left GFelement on the left of >= sign
         @returns true if left >= right, false otherwise
      */
      friend bool operator>=(const GFelement& left, const GFelement& right);

      friend bool operator<(const GFelement& left, const Fint& right);
      friend bool operator>(const GFelement& left, const Fint& right);
      friend bool operator<=(const GFelement& left, const Fint& right);
      friend bool operator>=(const GFelement& left, const Fint& right);

      friend bool operator<(const GFelement& left, const unsigned int& right);
      friend bool operator>(const GFelement& left, const unsigned int& right);
      friend bool operator<=(const GFelement& left, const unsigned int& right);
      friend bool operator>=(const GFelement& left, const unsigned int& right);


      /**
         Adds two GFelements
         @param right GFelement on the right of + sign
         @param left GFelement on the left of + sign
         @returns sum of two operands
      */
      friend GFelement operator+(const GFelement& left, const GFelement& right);


      /**
         Subtracts two GFelements
         @param right GFelement on the right of - sign
         @param left GFelement on the left of - sign
         @returns difference of two operands
      */
      friend GFelement operator-(const GFelement& left, const GFelement& right);


      /**
         Multiplies two GFelements
         @param right GFelement on the right of * sign
         @param left GFelement on the left of * sign
         @returns product of two operands
      */
      friend GFelement operator*(const GFelement& left, const GFelement& right);


      /**
         Divides two GFelements
         @param right GFelement on the right of / sign
         @param left GFelement on the left of / sign
         @returns quotient of two operands
      */
      friend GFelement operator/(const GFelement& left, const GFelement& right);


      /**
         Multiplies a GFelement by a scalar (i.e. a Fint)
         @param right GFelement on the right of * sign
         @param left Fint on the left of * sign
         @returns product of two operands
      */
      friend GFelement operator*(Fint left, const GFelement& right);


      /**
         Multiplies a GFelement by a scalar (i.e. a Fint)
         @param right Fint on the right of * sign
         @param left GFelement on the left of * sign
         @returns product of two operands
      */
      friend GFelement operator*(const GFelement& left, Fint right);


      /**
         Divides a GFelement by a scalar (i.e. a Fint)
         @param right Fint on the right of / sign
         @param left GFelement on the left of / sign
         @returns quotient of two operands
         @throws ErrorDivideByZero
      */
      friend GFelement operator/(const GFelement& left, Fint right);


      /**
         Outputs a GFelement to the standard output
         @param right GFelement on the right of << sign
         @param output an output stream on the left of << sign
         @returns output stream
      */
      friend ostream& operator<<(ostream& output, const GFelement& right);


      /**
         Determines the multiplicative inverse of a GFelement.
         @param gfe input GFelement
         @returns GFelement (the multiplicative inverse of gfe)
         @throws ErrorNoInverse, ErrorInvalidArraySize
      */
      friend GFelement GFmultInverse(const GFelement& gfe);


      /**
         Raises a GFelement G to a positive integer power m

         It uses the repeated square-and-multiply algorithm for exponentiation in GF(p^k).

         (See Chapter 2, Algorithm 2.227 in Handbook of Applied Cryptography,
         by Menezes, van Oorschot, Vanstone, CRC Press, 1996.)

         @param G an input GFelement
         @param m an input positive integer
         @returns G raised to the power m
         @throws ErrorInvalidArraySize, ErrorDivideByZero, ErrorIncompatibleFields
      */
      friend GFelement GFexp(const GFelement& G, Fint m);


   private:

      GaloisField* field;  // Galois Field associated with this polynomial
      Fint value;    // value for this GFelement

};



//------------------------------------------------
// Inline class functions
//------------------------------------------------

inline GaloisField* GFelement::getField() const
{
   return field;
}


inline Fint GFelement::getX() const
{
   return value;
}


inline void GFelement::setX(Fint value)
{
   if (field)
      this->value = value % field->getModulus();
   else this->value = value;
}

inline void GFelement::setField(GaloisField* gf) {
   field = gf;
   this->value = value % field->getModulus();
}

inline GFelement GFelement::getRootOfUnity() const
{
   GFelement temp{field, field->getRootOfUnity()};
   return temp;
}


inline GFelement GFelement::getZero() const
{
   Fint num = 0;
   GFelement temp{field, num};
   return temp;
}


inline GFelement GFelement::getOne() const
{
   Fint num = 1;
   GFelement temp{field, num};
   return temp;
}

inline GFelement& GFelement::operator++()
{
   ++value;
   if (value >= field->getModulus())
      value %= field->getModulus();
   return *this;
}

} // namespace shk_galoiscpp

#endif
