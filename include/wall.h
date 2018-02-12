#ifndef WALL_H
#define WALL_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "particle.h"
#include "VEC3F.h"
#include "field_3D.h"

using namespace std;

class wall {
public:
  wall();
  wall(const VEC3F& normal, const VEC3F& point);

  // draw to OGL
  void draw();
  void createwall( double , double, vector<wall>& );

  // accessors
  VEC3F& getNormal() { return _normal; }
  VEC3F& getPoint()  { return _point; }

private:

  VEC3F _normal;
  VEC3F _point;
};

#endif
