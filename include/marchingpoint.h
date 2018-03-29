#ifndef MARCHINGPOINT_H
#define MARCHINGPOINT_H

#include "vec3f.h"


class MarchingPoint {
public:
    MarchingPoint(){}
    MarchingPoint(float x, float y, float z) : _position(VEC3F(x,y,z)), _color(VEC3F(0,0,0)), _colorFieldValue(0.0){}

    VEC3F& getColor() { return _color;}
    VEC3F& getColorFieldValue() { return _colorFieldValue;}
    VEC3F& getPosition() { return _position;}
    void setColor(float r) { _color = VEC3F(r,0,0);}
    void setColorFieldValue( float cfv) { _colorFieldValue = cfv;}
    void updateColorInfo( float cfv ) { _color = VEC3F(cfv,0,0); _colorFieldValue = cfv;}

private:
    VEC3F _position;
    VEC3F _color;
    VEC3F _colorFieldValue;
};
#endif // MARCHINGPOINT_H
