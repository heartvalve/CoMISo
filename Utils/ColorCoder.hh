/*===========================================================================*\
 *                                                                           *
 *                              CoMISo                                       *
 *      Copyright (C) 2008-2009 by Computer Graphics Group, RWTH Aachen      *
 *                      www.graphics.rwth-aachen.de                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is a part of CoMISo.                                          *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 1    $                                                       *
 *   $Author: zimmer $                                                       *
 *   $Date: 2009-08-05 16:35:37 +0200 (Wed, 05 Aug 2009) $                   *
 *                                                                           *
\*===========================================================================*/


//=============================================================================
//
//  CLASS ColorCoder
//
//=============================================================================


#ifndef ACG_COLORCODER_HH
#define ACG_COLORCODER_HH

#ifdef WIN32
#undef min
#undef max
#endif


//== INCLUDES =================================================================


#include <ACG/Math/VectorT.hh>


//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================


class ColorCoder
{
public:
   
  /// Default constructor.
  ColorCoder(float _min=0.0, float _max=1.0, bool _signed=false) 
  { set_range(_min, _max, _signed); }
 


  /// set the color coding range for unsigned coding
  void set_range(float _min, float _max, bool _signed) 
  {
    if (_min == _max)
    {
      val0_ = val1_ = val2_ = val3_ = val4_ = _min;
    }
    else
    {
      if (_min > _max) std::swap(_min, _max);
      val0_ = _min;
      val4_ = _max;
      val2_ = 0.5f * (val0_ + val4_);
      val1_ = 0.5f * (val0_ + val2_);
      val3_ = 0.5f * (val2_ + val4_);
    }
    signed_mode_ = _signed;
  }


  /// color coding
  ACG::Vec3uc color(float _v) const {
    return signed_mode_ ? color_signed(_v) : color_unsigned(_v);
  }

  /// color coding
  ACG::Vec3f color_float(float _v) const {
    ACG::Vec3uc c;
    if(signed_mode_) c=color_signed(_v); else c=color_unsigned(_v);
    return (ACG::Vec3f(c[0],c[1],c[2])/255.f);
  }


  /// min scalar value
  float min() const { return val0_; }
  /// max scalar value
  float max() const { return val4_; }


private:


  ACG::Vec3uc color_unsigned(float _v) const 
  {
    if (val4_ <= val0_) return ACG::Vec3uc(0, 0, 255);

    unsigned char u;

    if (_v < val0_) return ACG::Vec3uc(0, 0, 255);
    if (_v > val4_) return ACG::Vec3uc(255, 0, 0);


    if (_v <= val2_) 
    {
      // [v0, v1]
      if (_v <= val1_) 
      {
	u = (unsigned char) (255.0 * (_v - val0_) / (val1_ - val0_));
	return ACG::Vec3uc(0, u, 255);
      }
      // ]v1, v2]
      else 
      {
	u = (unsigned char) (255.0 * (_v - val1_) / (val2_ - val1_));
	return ACG::Vec3uc(0, 255, 255-u);
      }
    }
    else 
    {
      // ]v2, v3]
      if (_v <= val3_) 
      {
	u = (unsigned char) (255.0 * (_v - val2_) / (val3_ - val2_));
	return ACG::Vec3uc(u, 255, 0);
      }
      // ]v3, v4]
      else 
      {
	u = (unsigned char) (255.0 * (_v - val3_) / (val4_ - val3_));
	return ACG::Vec3uc(255, 255-u, 0);
      }
    }
  }


  ACG::Vec3uc color_signed(float _v) const 
  {
    if (val4_ <= val0_) return ACG::Vec3uc(0,255,0);

    unsigned char r,g,b;

    if      (_v < val0_) _v = val0_;
    else if (_v > val4_) _v = val4_;

    if (_v < 0.0) 
    {
      r = val0_ ? (unsigned char)(255.0 * _v / val0_) : 0;
      b = 0;
    }
    else 
    {
      r = 0;
      b = val4_ ? (unsigned char)(255.0 * _v / val4_) : 0;
    }
    g = 255 - r - b;

    return ACG::Vec3uc(r, g, b);
  }


  
  
  float  val0_, val1_, val2_, val3_, val4_;
  bool   signed_mode_;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_COLORCODER_HH defined
//=============================================================================

