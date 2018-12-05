#include "rigid_body.h"

RigidBody::RigidBody(const double d, Circle& c): _density(d), _geometry(c),
	_linearVelocity(0,0), _angularMomentum(0), _angularVelocity(0){

	_localFrame.set_position(VEC2F(0,0));
	_localFrame.set_orientation(0);

	_inertia = 0;
	_geometry.get_inertia(_density, _inertia);
	_mass = _geometry.get_mass(_density);
	_inertia_inv = 1.0/_inertia;
}

void RigidBody::getWorldInvInertia(double& worldInv) const{
	worldInv = _inertia_inv;
}

CoordFrame RigidBody::getCoordFrame(void) const{
	return _localFrame;
}

void RigidBody::setCoordFrame(CoordFrame f){
	_localFrame = f;
}

const double RigidBody::getMass(void) const{
	return _mass;
}
const double RigidBody::getDensity(void) const{
	return _density;
}

const double RigidBody::getInertiaModulus(void) const{
	return _inertia;
}

const double RigidBody::getInvInertiaModulus(void) const{
	return _inertia_inv;
}

void RigidBody::getCOM(VEC2F& v) const{
	v = _localFrame.get_position();
}

void RigidBody::setCOM(const VEC2F& v){
	_localFrame.set_position(v);
}

double RigidBody::getAngle(void){
	return _localFrame.get_orientation();
}

void RigidBody::setAngle(double a){
	_localFrame.set_orientation(a);
}

void RigidBody::getLinearVelocity(VEC2F& v) const{
	v = _linearVelocity;
}

void RigidBody::setLinearVelocity(const VEC2F& v){
	_linearVelocity = v;
}

void RigidBody::getAngularMomentum(double& d) const{
	d = _angularMomentum;
}

void RigidBody::setAngularMomentum(const double d){
	_angularMomentum = d;
}

VEC2F RigidBody::getPointVelocity(const VEC2F& pt){
	VEC2F linV, cen;
	double omega;
	getLinearVelocity(linV); getAngularVelocity(omega); getCOM(cen);

	VEC3F r(pt[0]-cen[0],pt[1]-cen[1],0);
	VEC3F pt_vel = VEC3F(linV[0],linV[1],0) + (VEC3F(0,0,omega)).cross(r);
	return VEC2F(pt_vel[0],pt_vel[1]);
}

const bool RigidBody::isInside(const VEC2F& pt) const{
	return _geometry.is_inside(_localFrame.global_to_local(pt));
}

const bool RigidBody::testCollisionAndProject(const VEC2F& pt, VEC2F& new_pt){
	VEC2F v = _localFrame.global_to_local(pt);
	assert(!isnan((new_pt)(0)));
	_geometry.project_out(v);
	assert(!isnan((new_pt)(0)));
	new_pt = _localFrame.local_to_global(v);
	assert(!isnan((new_pt)(0)));
	return true;
}

const double RigidBody::getSignedDist(const VEC2F& pt) const{
	VEC2F v = _localFrame.global_to_local(pt);
	return _geometry.get_signed_dist(v);
}

void RigidBody::getVertices(std::vector<VEC2F>& v){
	_geometry.get_vertices(v);
	for(auto it = v.begin();it!=v.end();++it){
		*it = _localFrame.local_to_global(*it);
	}
}

void RigidBody::advance(const double dt){
	updatePosition(dt);
	updateOrientation(dt);

	applyForce(VEC2F(0,-9.81),dt);
	applyTorque(0.0,dt);
}

void RigidBody::updatePosition(const double dt){
	_localFrame.set_position(_localFrame.get_position() + _linearVelocity*dt);
}
void RigidBody::updateOrientation(const double dt){
	_localFrame.rotate(_angularVelocity*dt);
	recomputeAngularVelocity();
}

void RigidBody::applyForce(const VEC2F& f, const double dt){
	_linearVelocity += f*dt;
}

void RigidBody::applyTorque(double T, const double dt){
	_angularMomentum += T*dt;
	recomputeAngularVelocity();
}

void RigidBody::getAngularVelocity(double& w){
	recomputeAngularVelocity();
	w = _angularVelocity;
}

void RigidBody::recomputeAngularVelocity(void){
	_angularVelocity = _inertia_inv * _angularMomentum;
}


