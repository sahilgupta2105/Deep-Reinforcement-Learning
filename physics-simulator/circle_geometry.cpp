#include "circle_geometry.h"
#include <cmath>

// relative to object's local frame of reference
// NOTE: center of circle at (0,0)
bool Circle::is_inside(const VEC2F& pt) const{
	return (pt.norm() <= _radius);
}

double Circle::get_signed_dist(const VEC2F& pt) const{
	return (pt.norm() - _radius);
}

void Circle::project_out(VEC2F& pt){
	// check if point is inside the body
	if(is_inside(pt)){
		VEC2F v = pt - VEC2F(0,0);
		// calculate the angle of vector from x-axis
		double angle =0;
		if(pt[1]>0)
			angle = acos(v[1]/v.norm());
		else if (pt[1]<0)
			angle = 2*M_PI - acos(v[1]/v.norm());
		else{
			// this means 0 or pi
			if(pt[0]>=0)
				angle = 0;
			else
				angle = M_PI;
		}
		// the projection point is (rcos,rsin)
		pt[0] = _radius*cos(angle);
		pt[1] = _radius*sin(angle);
	}
}

double Circle::get_mass(double d){
	return d*M_PI*_radius*_radius;
}

void Circle::get_inertia(double d, double& i){
	const double mass = get_mass(d);
	i = (mass/2.0)*_radius*_radius;
}

void Circle::get_vertices(std::vector<VEC2F>& v){
	// for now we use an arbitrary number: 20
	for(int i=0;i<20;i++){
		double ang = 0 + i*M_PI/10;
		v.push_back(VEC2F(_radius*cos(ang),_radius*sin(ang)));
	}
}