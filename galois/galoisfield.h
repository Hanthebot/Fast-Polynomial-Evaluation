#ifndef GALOISFIELD_H
#define GALOISFIELD_H

#include "typedefs.h"
#include <cassert>


namespace shk_galoiscpp
{

/**
   This class defines a Galois Field of characteristic p and dimension k (for GF(p^k)).
   It stores the following (private) data about the field:
   @param modulus The modulus (aka the characteristics p), which MUST be a prime number.
   @param dimension The dimension (k) of Galois Field (a nonzero positive integer).
   @param rou The root of unity 

   The attributes of a GaloisField are set when it is instantiated.
   Once instantiated, the attributes of cannot be changed.
*/
class GaloisField
{
   public:

      //---------------------------
      // constructors
      //---------------------------
      /**
         Constructs a Galois Field (GF)
         @param mod a positive prime integer (the modulus of GF)
         @param dim a positive integer (the dimension of GF)
         @param rou
      */
      GaloisField(Fint mod, Int dim, Fint rou);
      GaloisField(Fint mod, Int dim, Fint rou, GFelement* dlog, GFelement* w);

      /**
         Destructs GaloisField, freeing allocated space for reduction polynomial in it.
      */
      // ~GaloisField();

      //---------------------------
      // accessors
      //---------------------------

      /**
         Returns modulus of GF.
      */
      Fint getModulus() const;

      /**
         Returns dimension of GF
      */
      Int getDimension() const;

      /**
         Returns roof of unity
      */
      Fint getRootOfUnity() const;

      /**
         Returns dlog (exponent cache)
      */
      GFelement* const getDlog() const;

      /**
         Returns w (exponent cache)
      */
      GFelement* const getW() const;

      /**
          
      */
     void setDlogW(GFelement* dlog, GFelement* w);

      /**
         Outputs information about this Galois Field to the standard output.
         @param gf GaloisField on the right of << sign
         @param output an output stream on the left of << sign
         @returns output stream
      */
      friend ostream& operator<<(ostream& output, const GaloisField& gf);

   private:
      Fint modulus;      // prime modulus (characteristic) of Galois Field
      Int dimension;     // dimension of Galois Field
      Fint rou;  // reduction polynomial of Galois Field
      GFelement* dlog; // precomputed values for inverse multiplication
      GFelement* w;    // precomputed values for inverse multiplication
};




//--------------------------------------
// Inline class functions
//--------------------------------------

inline Fint GaloisField::getModulus() const
{
   return modulus;
}


inline Int GaloisField::getDimension() const
{
   return dimension;
}


inline Fint GaloisField::getRootOfUnity() const
{
   return rou;
}

inline GFelement* const GaloisField::getDlog() const
{
   return dlog;
}

inline GFelement* const GaloisField::getW() const
{
   return w;
}

inline void GaloisField::setDlogW(GFelement* dlog, GFelement* w)
{
   assert(dlog != nullptr && w != nullptr);
   this->dlog = dlog;
   this->w = w;
}

} // namespace shk_galoiscpp

#endif