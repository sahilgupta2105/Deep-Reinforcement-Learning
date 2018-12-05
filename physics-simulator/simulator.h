#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "MAC.h"
#include "rigid_body.h"
#include "poisson_disk_sampling.h"

class Simulator{
private:
	int _nX, _nY;
	double _h;
	int _frame;

	// MAC grid for the simulator
	MAC_Grid _grid;

	// collection of fluid particles
	std::vector<VEC2F> _fluid_particles;
	double _fluid_radius;

	// controller info
	std::vector<double> _controller;

	RigidBody* _rbd;
	Circle* _geometry;

	// current time in simulator
	double sim_time;

	/////////////////
	double cfl(void);
	void advect_particles(double);
	void generate_fluid_particles(float , float , float , float , float, std::vector<int>&);
	void update_controller(double, std::vector<double>&);
	void add_particle(VEC2F);
	void controller(double, std::vector<double>&);
	void remove_particles(void);
	double rand_number(double, double);
	VEC2F rotate_matrix(const VEC2F&, double);
public:
	Simulator(int, int, double, double);
	~Simulator();
	std::vector<double> advance(double, std::vector<double>);
	std::vector<double> reset_simulator(void);
	std::vector<double > sample_velocity(void);
	void write_to_file(void);
};

#endif