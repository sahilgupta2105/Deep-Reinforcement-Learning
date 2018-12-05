#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "eigen_def.h"
#include "circle_geometry.h"
#include "coord_frame.h"

// defines class for handling rigid bodies in the simulation

class RigidBody{
private:
	double _density;
	double _mass;
	Circle _geometry;
	CoordFrame _localFrame; // coordinate frame to track the rigid body
	VEC2F _linearVelocity;
	double _angularMomentum;
	double _angularVelocity;
	double _inertia;
	double _inertia_inv;
	//////////////
	void recomputeAngularVelocity(void);
public:
	RigidBody(const double, Circle&);
	// computes the inverse of inertia matrix at point in world frame
	void getWorldInvInertia (double&) const;
	
	CoordFrame getCoordFrame(void) const;
	void setCoordFrame(CoordFrame);

	const double getMass(void) const;
	const double getDensity(void) const;

	const double getInertiaModulus(void) const;
	const double getInvInertiaModulus(void) const;

	void getCOM(VEC2F&) const;
	void setCOM(const VEC2F&);

	double getAngle(void);
	void setAngle(double);

	void getLinearVelocity(VEC2F&) const;
	void setLinearVelocity(const VEC2F&);

	void getAngularMomentum(double& ) const;
	void setAngularMomentum(const double);

	VEC2F getPointVelocity(const VEC2F&);

	// point is relative to world frame
	const bool isInside(const VEC2F&) const;
	const bool testCollisionAndProject(const VEC2F&, VEC2F&);
	const double getSignedDist(const VEC2F&) const;
	void getVertices(std::vector<VEC2F>&);

	void advance(const double);

	void updatePosition(const double);
	void updateOrientation(const double);

	void applyForce(const VEC2F&, const double);
	void applyTorque(double, const double);

	void getAngularVelocity(double&);
};

#endif