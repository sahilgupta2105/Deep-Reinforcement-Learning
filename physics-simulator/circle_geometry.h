#ifndef CIRCLE_GEOMETRY_H
#define CIRCLE_GEOMETRY_H

#include "eigen_def.h"
#include <vector>

// defines the geometry for circle: will be used for representing
// the rigid body in simulation

class Circle{
private:
	double _radius; // origin: at center of circle
public:
	Circle(double r): _radius(r) {}
	bool is_inside(const VEC2F&) const; // checks if given point inside body
	double get_signed_dist(const VEC2F&) const; // represents the SDF for circle
	void project_out(VEC2F&); // in case of collisions, projects the point to the boundary of rigid body
	double get_mass(double);
	void get_inertia(double, double&);
	void get_vertices(std::vector<VEC2F>&);
};

#endif