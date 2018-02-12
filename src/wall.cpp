#include "../include/wall.h"

const double thickness= 0.02;
///////////////////////////////////////////////////////////////////////////////
// Empty Constructor
///////////////////////////////////////////////////////////////////////////////
wall::wall()
{}
///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
wall::wall(const VEC3F& normal, const VEC3F& point) :
  _normal(normal), _point(point)
{
  // just in case, normalize the normal vector
  _normal.normalize();
}
///////////////////////////////////////////////////////////////////////////////
// Create the boundary condition wall
///////////////////////////////////////////////////////////////////////////////
void wall::createwall( double boundarysize, double k , vector<wall>& _walls)
{
    // prepare the 3d grid dimension
    VEC3F boxSize;
    FIELD_3D* grid;

    boxSize.x = boundarysize*2.0;
    boxSize.y = boundarysize;
    boxSize.z = boundarysize/2.0;

    int gridXRes = (int)ceil(boxSize.x/k);
    int gridYRes = (int)ceil(boxSize.y/k);
    int gridZRes = (int)ceil(boxSize.z/k);

    grid = new FIELD_3D(gridXRes, gridYRes, gridZRes);

   _walls.push_back(wall(VEC3F(0,0,1), VEC3F(0,0,-boxSize.z/2.0)));  // back
   _walls.push_back(wall(VEC3F(0,0,-1), VEC3F(0,0,boxSize.z/2.0)));  // front
   _walls.push_back(wall(VEC3F(1,0,0), VEC3F(-boxSize.x/2.0,0,0)));  // left
   _walls.push_back(wall(VEC3F(-1,0,0), VEC3F(boxSize.x/2.0,0,0)));  // right
   _walls.push_back(wall(VEC3F(0,1,0), VEC3F(0,-boxSize.y/2.0,0)));  // bottom

    cout << "Create the boundary condations" << endl;
    cout << "Grid size is " << (*grid).xRes() << "x" << (*grid).yRes() << "x" << (*grid).zRes() << endl;

}
///////////////////////////////////////////////////////////////////////////////
// OGL drawing
///////////////////////////////////////////////////////////////////////////////
void wall::draw()
{
  glPushMatrix();
  // translate to the point
  glTranslated(_point[0], _point[1], _point[2]);
    
    // apply a rotation
  double angle1 = asin(_normal[0]) / (2 * M_PI) * 360.0;
  double angle2 = asin(_normal[1]) / (2 * M_PI) * 360.0;
  //double angle3 = asin(_normal[2]) / (2 * M_PI) * 360.0;

  
    glRotatef(-angle1, 0, 1, 0);
  //cout << "1: " << angle1 << " 2: " << angle2 << " 3: " << angle3 << endl;
  glRotatef(-angle2, 1, 0, 0);

    // make it a plane at 0,0
    glTranslated(0, 0, thickness/2.0);
    glScalef(20,20,1);
    glutSolidCube(thickness);
  
    
  glPopMatrix();
}

