#ifndef MAC_H
#define MAC_H

#include <vector>
#include <algorithm>
#include <cmath>
#include "scalar_field.h"
#include "rigid_body.h"

const bool VERBOSE = 0;

class MAC_Grid{
private:
	int _nX, _nY;
	double  _h, _net_rigid_u, _net_rigid_v;
	ScalarF _fluid_u;
	ScalarF _fluid_u_valid;
	ScalarF _fluid_v;
	ScalarF _fluid_v_valid;
	ScalarF _fluid_p;
	ScalarF _fluid_phi;
	ScalarF _rigid_phi;
	ScalarF _boundary_phi;
	ScalarF _boundary_u, _boundary_v;
	ScalarF _rigid_u_weight;
	ScalarF _rigid_v_weight;
	ScalarF _fluid_u_weight;
	ScalarF _fluid_v_weight;
	//////////////
	double interpolate_value(const VEC2F&, const ScalarF&);
	double interpolate_gradient(VEC2F& ,const VEC2F&, const ScalarF&);
	void get_barycentric(double, int&, double&, int, int);
	double bilerp(double, double, double, double,double, double);
	double lerp(double, double, double);
	double fraction_in(double, double);
	double clamp(double, double, double);
	void set_domain_boundary(void);
	double sdf_domain(const VEC2F&);
	double U(const VEC2F&);
	double V(const VEC2F&);
	VEC2F rk2_tracer(const VEC2F&, double);
	void compute_pressure_weights(void);
	void solve_pressure(double, RigidBody*);
	void velocity_extrapolator(ScalarF&, ScalarF&);
	VEC2F velocity_boundary(const VEC2F&);
	MATRIX3 construct_star_matrix(const VEC3F&);
public:
	MAC_Grid(int , int , double);
	~MAC_Grid();
	// interpolate the value of velocity at any arbitrary position
	VEC2F velocity(const VEC2F&);
	// lookup for velocity field
	const ScalarF& getU(void);
	const ScalarF& getV(void);
	// accessor for velocity field
	ScalarF& modifyU(void);
	ScalarF& modifyV(void);
	// interpolate the value of 'phi' at any arbitrary position
	double boundary_phi(const VEC2F&, VEC2F&);

	void update_rigid_body_properties(const RigidBody*);
	void update_fluid_properties(const std::vector<VEC2F>&, const double);
	void advect_fluid_on_MAC(double);
	void force_on_fluid_MAC(double);
	void apply_projection(double, RigidBody*);
	void extrapolate_velocities(void);
	void recompute_boundary_velocity(RigidBody*);
	void constrain_velocity(void);
	void process_collisions_rigid_boundary(RigidBody*);
	void enforce_controller_bc(const std::vector<int>&, const double, const std::vector<double>&, double);
};

#endif