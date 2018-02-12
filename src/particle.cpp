#include "../include/particle.h"

VEC3F red(1,0,0);
VEC3F blue(0,0,1); 
VEC3F black(0,0,0);
VEC3F green(0,1,0);
VEC3F lightBlueColor(0.01,0.25,1.0);
VEC3F purpleColor(0.88,0.08,0.88);

int count = 0;

#define PARTICLE_DRAW_RADIUS 0.015 //0.015//0.01 //
#define h 0.0457
bool particle::isSurfaceVisible = false;
bool particle::showArrows = false;
bool particle::showSplash = false;
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
  /*
  if (_flag) 
    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
  else 
    glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
  */
  

  if (_flag && isSurfaceVisible)
    glMaterialfv(GL_FRONT, GL_DIFFUSE, purpleColor);
  else
    glMaterialfv(GL_FRONT, GL_DIFFUSE,blue);//1.5f * VEC3F(_position.x,_position.y,_position.z)

  //Since splash are surface red == surface
  if(_splash && showSplash)
         glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
  glPushMatrix();
    glTranslated(_position[0], _position[1], _position[2]);
  
    
    if (showArrows) {
      // scale
      
      //glColor3f(0.2f,0.3f,0.6f);
      if (!myQuadric) {
        myQuadric = gluNewQuadric();
        gluQuadricDrawStyle(myQuadric, GLU_FILL); 
        gluQuadricNormals(myQuadric, GLU_SMOOTH);
      }
      
      
      
      float angle1 = asin(_velocity[0]) * 180.0 / M_PI;
      float angle2 = asin(_velocity[1]) * 180.0 / M_PI;
      //float angle3 = asin(_velocity[2]) * 180.0 / M_PI;
      
      
      glRotatef(-angle1, 0, 1, 0);
      glRotatef(-angle2, 1, 0, 0);
      //glRotatef(-angle3, 0, 0, 1);
      
      gluCylinder(myQuadric, 0.001, 0.001, 0.01, 10, 10);
      glTranslated(0.00, 0.01, 0.00);
      glutSolidCone(0.003, 0.01, 10, 10);
      
      glFlush();        
      
      //;
    }
    else {
        if( _flag && isSurfaceVisible)
            glutWireSphere(PARTICLE_DRAW_RADIUS, 5, 5);
        else
            glutSolidSphere(PARTICLE_DRAW_RADIUS, 10, 10);

    }
  
  glPopMatrix();
}

void particle::clearParameters() {
  _position = VEC3F();
  _velocity = VEC3F();
  _acceleration = VEC3F();
  _density = 0.0;
  _pressure = 0.0;
  
}

