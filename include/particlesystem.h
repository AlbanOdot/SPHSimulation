#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

//#include <tr1/tuple>
//#include <map>
#include "particle.h"
#include "wall.h"
#include <vector>
#include "field_3D.h"
#include "simulation.h"

#define h 0.0457 //0.0457 0.02 //0.045

#define GAS_STIFFNESS 3.0 //20.0 // 461.5  // Nm/kg is gas constant of water vapor
#define REST_DENSITY 998.29 // kg/m^3 is rest density of water particle
#define CRITICAL_DENSITY 1300
#define SURFACE_TENSION 0.0728 // N/m
#define SURFACE_THRESHOLD 19.065
#define KERNEL_PARTICLES 20.0

#define GRAVITY_ACCELERATION -9.80665

#define WALL_K 10000.0 // wall spring constant
#define WALL_DAMPING -0.9 // wall damping constant

#define BOX_SIZE 0.4
#define MAX_PARTICLES 10000

#define INITIAL_SCENARIO SCENARIO_CUBE

using namespace std;

class particlesystem {

public:
    particlesystem();
    ~particlesystem();

    void updateGrid();

    // draw to OGL
    void draw();

    void addParticle(const VEC3F& position);

    void addParticle(const VEC3F& position, const VEC3F& velocity);

    void stepVerlet();

    void collisionForce(particle& particle, VEC3F& f_collision);

    float Wpoly6(float radiusSquared);

    void Wpoly6Gradient(VEC3F& diffPosition, float radiusSquared, VEC3F& gradient);

    float Wpoly6Laplacian(float radiusSquared);

    void WspikyGradient(VEC3F& diffPosition, float radiusSquared, VEC3F& gradient);

    float WviscosityLaplacian(float radiusSquared);

    void toggleGridVisble();

    void toggleSurfaceVisible();

    void toggleGravity();

    void toggleArrows();

    void toggleTumble();

    void toggleShowSpash();

    void generateFaucetParticleSet();

    void generateCubeParticleSet();

    void generateDamParticleSet();
    void fatCube();

    void makeItRain();

    void setGravityVectorWithViewVector(VEC3F viewVector);

    void densityAndPressureComputation();

    void accelerationComputation();

    void smoothTension();
    float C( float);
    //setters
    inline void scenario(const int scenario){ _scenario = scenario;}

    //getters
    inline int scenario() const { return _scenario;}
    void loadScenario(int scenario);

    FIELD_3D* grid;
    FIELD_3D* nextGrid;
    float surfaceThreshold;
    VEC3F gravityVector;


private:
    // list of particles, walls, and springs being simulated
    vector<particle> *_oldParticles;//Store current data
    vector<particle> *_newParticles;//Store next step data
    vector<wall>     _walls;
    wall     boundary;

    float viscosity;
    float particleMass;
    float dt;

    //unsigned int _particleCount;
    bool _isGridVisible;
    bool _tumble;

    VEC3F boxSize;

    int _scenario = INITIAL_SCENARIO;
};

#endif
