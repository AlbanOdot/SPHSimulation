#include "../include/particle.h"

VEC3F red(1,0,0);
VEC3F blue(0,0,1); 
VEC3F black(0,0,0);
VEC3F green(0,1,0);
VEC3F lightBlueColor(0.01,0.25,1.0);
VEC3F purpleColor(0.88,0.08,0.88);

int count = 0;

bool particle::isSurfaceVisible = false;
bool particle::showArrows = false;
bool particle::showSplash = false;
bool particle::display = true;
unsigned int particle::count = 0;

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////

particle::particle()
{
}

particle::particle(const VEC3F& position) :
  _position(position),_velocity(VEC3F()),_acceleration(VEC3F()),_mass(0.0457)
{
  myQuadric = NULL;
  _id = count++;
}

particle::particle(const VEC3F& position, const VEC3F& velocity) :
_position(position), _velocity(velocity), _acceleration(VEC3F()),_mass(0.0457)
{
  myQuadric = NULL;
  _id = count++;
}

///////////////////////////////////////////////////////////////////////////////
// OGL drawing
///////////////////////////////////////////////////////////////////////////////
void particle::draw()
{
  if(!display)
      return;

  if (_flag && isSurfaceVisible)
    glMaterialfv(GL_FRONT, GL_DIFFUSE, purpleColor);
  else
    glMaterialfv(GL_FRONT, GL_DIFFUSE,blue);//1.5f * VEC3F(_position.x,_position.y,_position.z)

  //Since splash are surface red == surface
  if(_splash && showSplash)
         glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
  glPushMatrix();
  glTranslated(_position[0], _position[1], _position[2]);
  glutSolidSphere(PARTICLE_DRAW_RADIUS, 10, 10);
  glPopMatrix();
}

void particle::clearParameters() {
  _position = VEC3F();
  _velocity = VEC3F();
  _acceleration = VEC3F();
  _density = 0.0;
  _pressure = 0.0;
  
}

