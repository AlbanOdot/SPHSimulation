#ifndef PARTICLE_H
#define PARTICLE_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "vec3f.h"
#include <vector>
#include <iostream>

#define PARTICLE_DRAW_RADIUS 0.015 //0.015//0.01 //
#define h 0.0457

using namespace std;

class particle {
      
public:
  
  static bool isSurfaceVisible;
  static bool showArrows;
  static bool showSplash;
  static bool display;
  
  //static unsigned int count;
  particle();
  particle(const VEC3F& position);
  particle(const VEC3F& position, const VEC3F& velocity);
  //~PARTICLE();
  
  // draw to OGL
  void draw();

  // clear all previous accumulated forces
  void clearForce() { _force *= 0; }
  // accumulate forces
  void addForce(VEC3F newForce) { _force += newForce; }
  
  void calculateAcceleration();

  // accessors
  inline VEC3F& position() { return _position; }
  inline VEC3F& velocity() { return _velocity; }
  inline VEC3F& acceleration() { return _acceleration; }
  inline VEC3F& force() { return _force; }
  inline VEC3F& surfaceTension() { return _surfaceTension;}
  inline VEC3F& normal() { return _normal;}
  inline float density()  const { return _density; }
  inline float pressure() const { return _pressure; }
  inline bool& flag() { return _flag; }
  inline bool& splash() { return _splash;}
  inline int id() const { return _id; }
  inline float mass() { return _mass;}

  
  //setters
  inline void setPosition(const VEC3F& pos){_position = pos;}
  inline void setVelocity(const VEC3F& vel){_velocity = vel;}
  inline void setAcceleration(const VEC3F& acc){ _acceleration = acc;}
  inline void setForce(const VEC3F& force){_force = force;}
  inline void setDensity(const float& density){_density = density;}
  inline void setPressure(const float& pressure){_pressure = pressure;}

  void clearParameters();
  

  //operators
  bool operator !=(const particle& p){
      return _id != p.id();
  }
  bool operator ==(const particle& p){
      return _id == p.id();
  }
  friend std::ostream& operator <<(std::ostream& stream,particle& p){
     return stream << "position : "<<p.position()<<"  velocity : "<<p.velocity()<<"  acceleration :"<<p.acceleration();
  }


  static unsigned int count;
  
private:  
  VEC3F _position;
  VEC3F _velocity;
  VEC3F _force;
  VEC3F _acceleration;
  VEC3F _normal;
  VEC3F _surfaceTension;
  float _density;
  float _pressure;
  float _mass;
  bool _flag;
  bool _splash;
  int _id;
  GLUquadricObj* myQuadric;
    
};
#endif
