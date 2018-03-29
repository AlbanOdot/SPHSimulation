#ifndef FIELD_3D_H
#define FIELD_3D_H

#include <cstdlib>
#include <iostream>
#include "assert.h"
#include "particle.h"

using namespace std;

template < class T = particle>
class FIELD_3D {
  
  typedef vector<T> TVector;
  
public:
  
  FIELD_3D():_xRes(0), _yRes(0), _zRes(0), _data(NULL){}

  FIELD_3D(int xRes, int yRes, int zRes) : _xRes(xRes), _yRes(yRes), _zRes(zRes), _cellCount(xRes*yRes*zRes)
  {
    _data = new TVector[_xRes * _yRes * _zRes];
  }

  virtual ~FIELD_3D()
  {
    if (_data) delete[] _data;
  }
  
  inline TVector& operator()(int x, int y, int z) {
    
    /*
     
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes); // i*length*width + j*width + k
    
     */
    return _data[x + y*_xRes + z*_xRes*_yRes];
  }
  
  // accessors
  int xRes() const { return _xRes; }
  int yRes() const { return _yRes; }
  int zRes() const { return _zRes; }
  int cellCount() const { return _cellCount; }
  TVector* data() const { return _data; }
  
private:
  
  int _xRes;
  int _yRes;
  int _zRes;
  int _cellCount;
  
  TVector* _data;
  
};

#endif
