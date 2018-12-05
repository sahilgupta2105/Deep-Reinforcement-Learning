#include "simulator.h"
#include "scalar_field.h"
#include <Partio.h>
#include <iostream>
#include <string>
#include <ctime>
#include <set>
#include <fstream>
#include <random>

Simulator::Simulator(int nX, int nY, double h, double r): _grid(nX,nY,h), _controller(6,0.0){
	_nX = nX;
	_nY = nY;
	_h = h;
	_fluid_particles.clear();
	_fluid_radius = r;

	_frame = 0;
	// rigid body initializing goes here
	_geometry = new Circle(0.1);
	_rbd = new RigidBody(1.5, *_geometry); // 1st argument is the density of solid
	_rbd->setCOM(VEC2F(0.8f,1.0f));
	_rbd->setAngularMomentum(0);
	_rbd->setLinearVelocity(VEC2F(0.0f,0.0f));
	_controller[0]=0.5;
}

Simulator::~Simulator() {};

double Simulator::cfl(void){
	ScalarF field = _grid.getU();
	double max_vel = 0;
	for (int i=0; i < field.size(); i++)
		max_vel = std::max(max_vel,std::abs(field(i)));
	field = _grid.getV();
	for (int i=0; i < field.size(); i++)
		max_vel = std::max(max_vel,std::abs(field(i)));
	return _h / max_vel;
}

std::vector<double> Simulator::reset_simulator(void){
	_frame = 0;
	_fluid_particles.clear();
	_grid = MAC_Grid(_nX,_nY,_h);
	_controller.clear();
	_rbd->setCOM(VEC2F(0.8f,1.0f));
	_rbd->setAngularMomentum(0);
	_rbd->setLinearVelocity(VEC2F(0.0f,0.0f));
	// generate a random state for ball
	std::vector<double> obs_state;
	double x_vel = rand_number(-0.01,0.01);
	double y_vel = rand_number(-0.01,0.01);
	double x_pos = rand_number(0.8-0.01,0.8+0.01);
	double y_pos = rand_number(1.0-0.01,1.0+0.01);
	double b_omg =rand_number(-0.01,0.01) ;
	_rbd->setCOM(VEC2F(x_pos,y_pos));
	_rbd->setLinearVelocity(VEC2F(x_vel,y_vel));
	obs_state.push_back(x_vel); // ball x-vel
	obs_state.push_back(y_vel);
	obs_state.push_back(b_omg);
	obs_state.push_back(x_pos); // ball pos
	obs_state.push_back(y_pos);
	_controller[0] = rand_number(0.4,1.2);
	obs_state.push_back(_controller[0]); // controller x_pos
	_controller[1] = rand_number(-0.1,0.1);
	obs_state.push_back(_controller[1]);
	_controller[3] = rand_number(-60*(M_PI/180),60*(M_PI/180));
	obs_state.push_back(_controller[3]);
	_controller[4] = rand_number(-1,1);
	obs_state.push_back(_controller[4]);
	// before starting the simulation we want to write the current state of system
	write_to_file();
	_frame +=1;
	return obs_state;
}

double Simulator::rand_number(double l, double u){
	double lower_bound = l;
	double upper_bound = u;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(lower_bound,upper_bound);
	double a_random_double =  dis(gen) ;
	return a_random_double;
}

std::vector<double> Simulator::advance(double dt, std::vector<double> control_action){
	double t = 0;
	double frame_t = (1.0/24.0); // fps: 24
	// till the time step is less than dt, we keep advancing
	clock_t begin, end;
	if(VERBOSE)
		std::cout<<"-------I am inside advance, starting to step in time.------"<<std::endl;
	controller(dt, control_action);
	// to visualize the shooting let's write data here as well
	write_to_file();
	_frame+=1;
	while(t < dt){
		if(VERBOSE)
			std::cout<<"-----Time History------"<<std::endl;
		begin = clock();
		double substep = cfl();
		end = clock();
		if(VERBOSE)
			std::cout<<"Cfl computation: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		if(t + substep > dt)
			substep = dt - t;
		if(t + 2*substep > frame_t)
			substep = (frame_t - t)/2;
		else if (t + substep > frame_t)
			substep = (frame_t - t);
		begin = clock();
		advect_particles(substep);
		end = clock();
		if(VERBOSE)
			std::cout<<"Fluid Particle advection: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		
		// integrate the forces for rigid body
		begin = clock();
		_rbd->advance(substep);
		end = clock();
		if(VERBOSE)
			std::cout<<"Ball advancing: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		
		// process collision between domain boundary and rigid body
		begin = clock();
		_grid.process_collisions_rigid_boundary(_rbd);
		end = clock();
		if(VERBOSE)
			std::cout<<"Ball collisions: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		// recompute signed distances and face areas for rigid body
		begin = clock();
		_grid.update_rigid_body_properties(_rbd);
		end = clock();
		if(VERBOSE)
			std::cout<<"Ball phi computation: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		
		// compute liquid signed distances
		begin = clock();
		_grid.update_fluid_properties(_fluid_particles,_fluid_radius);
		end = clock();
		if(VERBOSE)
			std::cout<<"Fluid phi computation: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		
		// advance the fluid velocity and add any forces on MAC
		begin = clock();
		_grid.advect_fluid_on_MAC(substep);
		end = clock();
		if(VERBOSE)
			std::cout<<"Advect fluid on MAC: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		begin = clock();
		_grid.force_on_fluid_MAC(substep);
		end = clock();
		if(VERBOSE)
			std::cout<<"Force on fluid: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
			
		// project the solution on a subspace such that it 
		// satisfies the incompressibility condition and 
		// the boundary conditions
		begin = clock();
		_grid.apply_projection(substep,_rbd);
		end = clock();
		if(VERBOSE)
			std::cout<<"Pressure projection: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		
		// extrapolate velocities from fluid domain to zero face area faces
		begin = clock();
		_grid.extrapolate_velocities();
		end = clock();
		if(VERBOSE)
			std::cout<<"Extrapolate velocities: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		
		// MAC grid approx. velocities for rigid body to be used in 
		// constrainted extrapolation
		begin = clock();
		_grid.recompute_boundary_velocity(_rbd);
		end = clock();
		if(VERBOSE)
			std::cout<<"Recompute boundary velocity: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		begin = clock();
		_grid.constrain_velocity();
		end = clock();
		if(VERBOSE)
			std::cout<<"Constrain velocity: "<<double(end-begin)/CLOCKS_PER_SEC<<std::endl;
		t += substep;
		remove_particles();
	}
	sim_time += t;
	if(_frame*(1/24) - sim_time <= 1e-3){
		write_to_file();
		_frame+=1;
	}
	std::vector<double> obs_state;
	VEC2F body_vel;
	double ang_vel;
	VEC2F body_com;
	_rbd->getLinearVelocity(body_vel);  // 2
	_rbd->getAngularVelocity(ang_vel);  // 1
	_rbd->getCOM(body_com); // 2
	// controller 4: pos, vel, angle, omega
	obs_state.push_back(body_vel(0));
	obs_state.push_back(body_vel(1));
	obs_state.push_back(ang_vel);
	obs_state.push_back(body_com(0));
	obs_state.push_back(body_com(1));
	obs_state.push_back(_controller[0]);
	obs_state.push_back(_controller[1]);
	obs_state.push_back(_controller[3]);
	obs_state.push_back(_controller[4]);
	std::vector<double> vel_state;
	vel_state = sample_velocity();
	obs_state.reserve(obs_state.size()+vel_state.size());
	obs_state.insert(obs_state.end(),vel_state.begin(),vel_state.end());
	return obs_state;
}

void Simulator::advect_particles(double dt){
	// we want to advect each particle
	for(auto it = _fluid_particles.begin(); it!= _fluid_particles.end();++it){
		// 2nd order Runge Kutta time stepping
		VEC2F pos_before = *it;
		VEC2F start_vel = _grid.velocity(pos_before);
		VEC2F pos_mid = pos_before + 0.5*dt*start_vel;
		VEC2F mid_vel = _grid.velocity(pos_mid);
		(*it) += dt*mid_vel;
		// particles can leave the domain
		// correction for that follows
		VEC2F normal;
		double phi = _grid.boundary_phi(*it,normal);
		// phi is 'negative' outside the boundary and 'positive' inside
		if(phi < 0){
			normal.normalize();
			(*it) -= phi*normal;
		}
		// check for if any fluid particle doesn't end up inside the rigid body
		_rbd->testCollisionAndProject((*it),(*it));
	}
	// BOTTLENECK: O(n2)
	for(auto it = _fluid_particles.begin();it!=_fluid_particles.end();++it){
		for(auto it1 =_fluid_particles.begin();it1!=_fluid_particles.end();++it1){
			if(it==it1)
				continue;
			double target_d = 0.5*_h;
			if((*it - *it1).norm() < target_d){
				// this means particles are too close
				VEC2F sep_d = *it - *it1;
				double current_d = sep_d.norm();
				sep_d.normalize();
				(*it) -= 0.5*(current_d - target_d)*sep_d;
				(*it1) += 0.5*(current_d - target_d)*sep_d;
			}
		}
	}
}

void Simulator::add_particle(VEC2F p){
	_fluid_particles.push_back(p);
}

void Simulator::remove_particles(void){

	if(VERBOSE)
		std::cout<<"Particle removal in action!"<<std::endl;
	for(auto it = _fluid_particles.begin(); it!= _fluid_particles.end();){
		if((*it)(1) < 0.3)
			it = _fluid_particles.erase(it);
		else
			++it;
	}
}

void Simulator::controller(double dt, std::vector<double>& action){

	// move the controller based on action
	update_controller(dt, action);

	// shoot only if required
	
	if(action[2]!=0){

		double x_loc = _controller[0];
		double y_loc = 0.3; 
		// hose dimensions
		double length = 0.1;
		double height = 0.1;
		std::vector<int> cell_ids; // 1-D id's: id_x + nX*id_y
		generate_fluid_particles(x_loc-length/2,x_loc+length/2,y_loc,y_loc+(height)*action[2],2*_fluid_radius, cell_ids);
		// apply force on fluid particles
		// its a two-step process:
		// a) find unique cells where particles are present
		// b) transfer velocities from cell to particles

		//duplicate removal from cell id vector
		std::set<int> s(cell_ids.begin(),cell_ids.end());
		cell_ids.clear();
		cell_ids.assign(s.begin(),s.end());
		
		// enforce boundary condition at every cell which contains newly created particles
		std::vector<double> control_state(2,0);
		control_state.push_back(_controller[1]);
		control_state.push_back(_controller[3]);
		_grid.enforce_controller_bc(cell_ids, dt, control_state, 1);

	}

}

void Simulator::update_controller(double dt, std::vector<double>& action){
	
	bool b = false;
	if(b){
		if(_frame==0){
			std::cout<<"----------------- NEW EPISODE--------------"<<std::endl;
			VEC2F dum_vel;
			_rbd->getLinearVelocity(dum_vel);
			std::cout<<"BALL_START: "<<dum_vel(0)<<", "<<dum_vel(1)<<std::endl;
		}
		std::cout<<"-------"<<std::endl;
		std::cout<<"DB: "<<_controller[0]<<", "<<_controller[1]<<", "<<_controller[2]<<", "<<_controller[3]<<", "<<_controller[4]<<", "<<_controller[5]<<std::endl;
		std::cout<<"DB(action): "<<action[0]<<", "<<action[1]<<", "<<action[2]<<", "<<action[3]<<std::endl;
	}
	_controller[2] = 1000*action[0];
	_controller[5] = 1325*action[1];

	double controller_x_max = 1.2;
	double controller_x_min = 0.4;

	_controller[1] += _controller[2]*dt;
	_controller[0] += _controller[1]*dt;
	if(_controller[0] >= controller_x_max){
		//  this means the controller collided with wall so velocity should also be forced to zero
		_controller[1] = 0;
		_controller[0] = controller_x_max;
	}
	if(_controller[0] <= controller_x_min){
		_controller[1] = 0;
		_controller[0] = controller_x_min;
	}

	_controller[4] += _controller[5]*dt;
	_controller[3] += _controller[4]*dt;
	if(_controller[3] >= (M_PI/180)*60){
		_controller[4] = 0;
		_controller[3] = (M_PI/180)*60;
	}
	if(_controller[3] <= -(M_PI/180)*60){
		_controller[4] = 0;
		_controller[3] = -(M_PI/180)*60;
	}

	if(b){
		std::cout<<"DB: "<<_controller[0]<<", "<<_controller[1]<<", "<<_controller[2]<<", "<<_controller[3]<<", "<<_controller[4]<<", "<<_controller[5]<<std::endl;
		std::cout<<"-------"<<std::endl;
	}
}

void Simulator::generate_fluid_particles(float x_min, float x_max, float y_min, float y_max, float r, std::vector<int>& id){
	float radius = r;
	if(y_max-y_min<r)
		return;
  	const auto bound_min = std::array<float, 2>{{ x_min, y_min }};
  	const auto bound_max = std::array<float, 2>{{ x_max, y_max }};
  	const auto samples = thinks::poisson_disk_sampling::PoissonDiskSampling(radius, bound_min, bound_max);
  	int i=1;
  	for(auto it=samples.begin();it!=samples.end();++it){
  		VEC2F v((*it)[0],(*it)[1]);
  		add_particle(v);
  		id.push_back(int(std::floor(v(0)/_h) + std::floor(v(1)/_h)*_nX));
  		i+=1;
  	}

  	if(VERBOSE)
  		std::cout<<"??????? Created "<<i<< " particles!"<<std::endl;
}

std::vector<double > Simulator::sample_velocity(void){
	std::vector<double> velX(196,0.0);
	std::vector<double> velY(196,0.0);
	std::vector<double> velMag(196, 0.0);
	
	// moving window sampling for fluid velocity
	// size: 14x14
	// get index of cell containing ball
	VEC2F ball_pos;
	_rbd->getCOM(ball_pos);

	int ball_x = std::floor(ball_pos(0)/_h);
	int ball_y = std::floor(ball_pos(1)/_h);

	int x_min = (ball_x - 3>=0)?ball_x - 3:0;
	int x_max = (ball_x + 3<=_nX-1)?ball_x + 3:_nX-1;
	int y_min = (ball_y - 3>=0)?ball_y - 3:0;
	int y_max = (ball_y + 3>=_nY-1)?ball_y + 3:_nY-1;

	// adjust for boundaries
	if(x_max-x_min!=6){
		if(x_max==_nX-1)
			x_min -= 3 - (x_max-ball_x);
		else
			x_max += 3 - (ball_x-x_min);
	}

	if(y_max-y_min!=6){
		if(y_max==_nY-1)
			y_min -= 3 - (y_max-ball_y);
		else
			y_max += 3 - (ball_y-y_min);
	}
	for(auto it = _fluid_particles.begin(); it!= _fluid_particles.end(); ++it){
		// find cell id to which the particle belongs
		int idxX = int(std::floor((*it)(0)/_h));
		int idxY = int(std::floor( (*it)(1) / _h ));
		int idxT = idxX + _nX*idxY;
		if((idxX>=x_min && idxX<=x_max) && (idxY>=y_min && idxY<=y_max)){
			int xidxX = idxX - x_min;
			int yidxY = idxY - y_min;
			int tidxT = xidxX + 14*yidxY;
			VEC2F vel = _grid.velocity(*it);
			velX.at(tidxT) = vel(0);
			velY.at(tidxT) = vel(1);
			velMag.at(tidxT) = std::sqrt(vel(0)*vel(0) + vel(1)*vel(1) );
		}
	}

	// normalize the values
	for(auto i = 0; i<196; i++){
		if(velMag[i]!=0){
			velX[i] = velX[i]/velMag[i];
			velY[i] = velY[i]/velMag[i];
		}
	}

	velX.reserve(2*velX.size());
	velX.insert(velX.end(),velY.begin(),velY.end());
	return velX;
}

VEC2F Simulator::rotate_matrix(const VEC2F& pt, double a) {
		// some coord transform w.r.t rotation origin: (_controller[0], 0.3)
		double pt_x = pt(0) - _controller[0];
		double pt_y = pt(1) - 0.3;
		double c = std::cos(a);
		double s = std::sin(a);
		VEC2F r_pt(c*pt_x - s*pt_y, s*pt_x + c*pt_y);
		return VEC2F(r_pt(0)+_controller[0],r_pt(1)+0.3);
}


void Simulator::write_to_file(void){

	Partio::ParticlesDataMutable* parts = Partio::create();
	Partio::ParticleAttribute posH, vH, mH;
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	for (auto it = _fluid_particles.begin(); it!=_fluid_particles.end(); ++it){
		int idx = parts->addParticle();
		float* p = parts->dataWrite<float>(posH, idx);
		assert(!isnan((*it)(0)));
		assert(!isnan((*it)(1)));
		for (int k = 0; k < 2; k++)
			p[k] = (float)(*it)[k];
		p[2] = (float)(0); //otherwise writes random numbers for z-direction
	}
	std::string filename = "data/fluid" + std::to_string(_frame) + ".bgeo";
	Partio::write(filename.c_str() , *parts);
	//////// RIGID BODY
	parts = Partio::create();
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	std::vector<VEC2F> rbd_vert;
	VEC2F rbd_com;
	_rbd->getCOM(rbd_com);
	rbd_vert.push_back(rbd_com);
	for (auto it = rbd_vert.begin(); it!=rbd_vert.end(); ++it) {
		int idx = parts->addParticle();
		float* p = parts->dataWrite<float>(posH, idx);
		for (int k = 0; k < 2; k++)
			p[k] = (float)(*it)[k];
		p[2] = (float)(0);
	}
	filename = "data/rigid" + std::to_string(_frame) + ".bgeo";
	Partio::write(filename.c_str(), *parts);
	/////////// CONTROLLER
	parts = Partio::create();
	// we are only pushing the position and velocity for now
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	std::vector<VEC2F> control_vert;
	// first compute corners
	VEC2F c1, c2, c3, c4;
	c1(0) = _controller[0] - 0.1/2; c1(1) = 0.3;
	c2(0) = _controller[0] + 0.1/2; c2(1) = 0.3;
	c3(0) = _controller[0] - 0.1/2; c3(1) = 0.4;
	c4(0) = _controller[0] + 0.1/2; c4(1) = 0.4;
	// now rotate the points
	double ang = -_controller[3];
	//c1 = rotate_matrix(c1,ang);
	//c2 = rotate_matrix(c2,ang);
	//c3 = rotate_matrix(c3, ang);
	//c4 = rotate_matrix(c4,ang);
	//control_vert.push_back(c1);
	//control_vert.push_back(c2);
	//control_vert.push_back(c3);
	//control_vert.push_back(c4);
	//control_vert.push_back();
	// also add the center of controller
	control_vert.push_back(VEC2F(_controller[0],0.3));
	for (auto it = control_vert.begin(); it!=control_vert.end(); ++it) {
		int idx = parts->addParticle();
		float* p = parts->dataWrite<float>(posH, idx);
		for (int k = 0; k < 2; k++)
			p[k] = (float)(*it)[k];
		p[2]= (float)(0.0);
	}
	filename = "data/control" + std::to_string(_frame) + ".bgeo";
	Partio::write(filename.c_str(), *parts);

	parts->release();
	parts = NULL;
}