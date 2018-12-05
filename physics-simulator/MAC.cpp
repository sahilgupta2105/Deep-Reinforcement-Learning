#include "MAC.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <set>

MAC_Grid::MAC_Grid(int nX, int nY, double h): _fluid_u(nX+1,nY,h),_fluid_v(nX,nY+1,h),
_fluid_p(nX,nY,h), _fluid_phi(nX,nY,h), _rigid_phi(nX+1,nY+1,h), _rigid_u_weight(nX+1, nY, h),
_rigid_v_weight(nX,nY+1,h), _fluid_u_weight(nX+1,nY,h), _fluid_v_weight(nX,nY+1,h),
_fluid_u_valid(nX+1,nY,h),_fluid_v_valid(nX,nY+1,h), _boundary_u(nX+1,nY,h),
_boundary_v(nX,nY+1,h), _boundary_phi(nX+1,nY+1,h){
	//////////////////////////////
	_nX = nX;
	_nY = nY;
	_h = h;
	_net_rigid_u = 0;
	_net_rigid_v = 0;

	// initialize the boundary here
	// uses a rectangle sdf function based on nX, nY and h
	set_domain_boundary();
}

MAC_Grid::~MAC_Grid() {};

void MAC_Grid::set_domain_boundary(void){
	for(int j=0; j < _boundary_phi.nY(); ++j){
		for(int i=0; i < _boundary_phi.nX(); ++i){
			_boundary_phi.update_field_val_at(i,j,sdf_domain(VEC2F(i*_h,j*_h)));
		}
	}
}

double MAC_Grid::sdf_domain(const VEC2F& pt){
	// rectangular sdf
	std::vector<double> negatives;
	double length = 1.0;
	double height = 1.4;
	VEC2F center_loc = VEC2F(0.8,0.8);
	VEC2F delta = pt - center_loc;
	double x_l = -length/2 - delta(0);
	double x_r = delta(0) - length/2;
	double y_b = -height/2 - delta(1);
	double y_t = delta(1) - height/2;
	if(x_l>=0)
		negatives.push_back(x_l);
	if(x_r>=0)
		negatives.push_back(x_r);
	if(y_b>=0)
		negatives.push_back(y_b);
	if(y_t>=0)
		negatives.push_back(y_t);

	if(negatives.size()==0)
		return -std::max(std::max(x_l,x_r),std::max(y_b,y_t));
	else if(negatives.size()==1)
		return -negatives[0];
	else if (negatives.size()==2)
		return -std::sqrt(negatives[0]*negatives[0] + negatives[1]*negatives[1]);
	else{
		std::cout<<"Crap: something wierd happened in the boundary SDF function!"<<std::endl;
		return 0;
	}
	
	/*// circle sdf
	VEC2F c(0.5,0.5);
	double r = 0.45;
	return -((pt - c).norm() - r);*/
	
}

double MAC_Grid::U(const VEC2F& pos) {return interpolate_value(pos / _h - VEC2F(0,0.5),_fluid_u);}


double MAC_Grid::V(const VEC2F& pos) {return interpolate_value(pos / _h - VEC2F(0.5,0),_fluid_v);}

VEC2F MAC_Grid::velocity(const VEC2F& pos){
	return VEC2F(U(pos),V(pos));
}

void MAC_Grid::recompute_boundary_velocity(RigidBody* r){
	for(int i =0; i < _nX+1; ++i){
		for(int j=0; j<_nY-1; ++j){
			VEC2F pos(i*_h,(j+0.5)*_h);
			if(0.5*(_boundary_phi(i,j)+_boundary_phi(i,j+1)) < 0.5*(_rigid_phi(i,j)+_rigid_phi(i,j+1)))
				_boundary_u.update_field_val_at(i,j,0);
			else{
				_boundary_u.update_field_val_at(i,j,r->getPointVelocity(pos)(0));
			}
		}
	}

	for(int i =0; i < _nX-1; ++i){
		for(int j=0; j<_nY+1; ++j){
			VEC2F pos((i+0.5)*_h,j*_h);
			if(0.5*(_boundary_phi(i,j)+_boundary_phi(i+1,j)) < 0.5*(_rigid_phi(i,j)+_rigid_phi(i+1,j)))
				_boundary_v.update_field_val_at(i,j,0);
			else{
				_boundary_v.update_field_val_at(i,j,r->getPointVelocity(pos)(1));
			}
		}
	}
}

void MAC_Grid::constrain_velocity(void){

	// for extrapolated points from the previous step in Simulator::advance,
	// we are replacing the normal component of velocity with the Ball's
	// velocity (here it's zero)

	for(int j=0; j < _fluid_u.nY(); ++j){
		for(int i=0; i < _fluid_u.nX(); ++i){
			if(_fluid_u_weight(i,j)==0){
				// apply constraint
				VEC2F pos(i*_h, (j+0.5)*_h);
				VEC2F vel = velocity(pos);
				VEC2F vel_boundary = velocity_boundary(pos);

				VEC2F boundary_normal(0,0);
				double boundary_phi = interpolate_gradient(boundary_normal, pos/_h, _boundary_phi);
				boundary_normal.normalize();

				VEC2F rigid_normal(0,0);
				double rigid_phi = interpolate_gradient(rigid_normal, pos/_h, _rigid_phi);
				rigid_normal.normalize();

				if(boundary_phi < rigid_phi){
					vel -= vel.dot(boundary_normal) * boundary_normal;
					vel += vel_boundary.dot(boundary_normal) * boundary_normal;
					_fluid_u.update_field_val_at(i,j,vel(0));
				}else{
					vel -= vel.dot(rigid_normal) * rigid_normal;
					vel += vel_boundary.dot(rigid_normal) * rigid_normal;
					_fluid_u.update_field_val_at(i,j,vel(0));
				}
			}
		}
	}

	for(int j=0; j < _fluid_v.nY(); ++j){
		for(int i=0; i < _fluid_v.nX(); ++i){
			if(_fluid_v_weight(i,j)==0){
				// apply constraint
				VEC2F pos((i+0.5)*_h, j*_h);
				VEC2F vel = velocity(pos);
				VEC2F vel_boundary = velocity_boundary(pos);

				VEC2F boundary_normal(0,0);
				double boundary_phi = interpolate_gradient(boundary_normal, pos/_h, _boundary_phi);
				boundary_normal.normalize();

				VEC2F rigid_normal(0,0);
				double rigid_phi = interpolate_gradient(rigid_normal, pos/_h, _rigid_phi);
				rigid_normal.normalize();

				if(boundary_phi < rigid_phi){
					vel -= vel.dot(boundary_normal) * boundary_normal;
					vel += vel_boundary.dot(boundary_normal) * boundary_normal;
					_fluid_v.update_field_val_at(i,j,vel(1));
				}else{
					vel -= vel.dot(rigid_normal) * rigid_normal;
					vel += vel_boundary.dot(rigid_normal) * rigid_normal;
					_fluid_v.update_field_val_at(i,j,vel(1));
				}
			}
		}
	}
}

VEC2F MAC_Grid::velocity_boundary(const VEC2F& pos){
	double val_u = interpolate_value(pos/_h - VEC2F(0,0.5), _boundary_u);
	double val_v = interpolate_value(pos/_h - VEC2F(0.5,0), _boundary_v);
	return VEC2F(val_u,val_v);
}

void MAC_Grid::process_collisions_rigid_boundary(RigidBody* r){
	// as coolisions are not a focus of this simulator,
	// a hacky way is used to handle impulse responses.
	// body is projected out of the walls
	std::vector<VEC2F> vertices;
	r->getVertices(vertices);
	double min_phi = _h;
	int min_ind = -1;
	VEC2F min_normal;
	std::vector<unsigned int> penetrating_pts;
	for(unsigned int i =0; i < vertices.size(); ++i){
		double phi;
		VEC2F normal;
		VEC2F vertex(vertices[i](0),vertices[i](1));
		phi = interpolate_gradient(normal, (vertex-VEC2F(0,0))/_h,_boundary_phi);
		normal.normalize();
		if(phi < 0){
			VEC2F pt_vel = r->getPointVelocity(vertex);
			if(pt_vel.dot(normal) < 1e-2){
				penetrating_pts.push_back(i);
			}
			if(phi < min_phi){
				min_phi = phi;
				min_ind = i;
				min_normal = normal;
			}
		}
	}

	int itr = 0;
	while(penetrating_pts.size()>0){
		itr++;

		int p = penetrating_pts[0];
		penetrating_pts.pop_back();

		VEC2F pt = vertices[p];
		VEC2F com;
		r->getCOM(com);
		VEC2F rad = pt - com;

		VEC2F vel = r->getPointVelocity(pt);
		VEC2F normal;
		interpolate_gradient(normal, (pt - VEC2F(0,0))/_h, _boundary_phi);
		normal.normalize();

		// projecting the normal velocity: inelastic collision
		double normal_vel = vel.dot(normal);
		// Compute impulse
		VEC3F rad3d(rad(0),rad(1),0);
		VEC3F normal3d(normal(0),normal(1),0);
		MATRIX3 K = MATRIX3::Zero();
		K(0,0) = 1/r->getMass();
		K(1,1) = 1/r->getMass();
		K(2,2) = 1/r->getMass();
		MATRIX3 R = construct_star_matrix(rad3d);
		MATRIX3 Iinv;
		double imodinv;
		imodinv = r->getInvInertiaModulus();
		Iinv(0,0) = 0;
		Iinv(1,1) = 0;
		Iinv(2,2) = imodinv;
		K = K + R.transpose()*Iinv*R;
		double j = -normal_vel / (normal3d.dot(K*normal3d));

		//update linear velocity
		VEC2F lin_vel;
		r->getLinearVelocity(lin_vel);
		lin_vel += j/r->getMass() * normal;
		r->setLinearVelocity(lin_vel);
		//update angular velocity
		VEC3F momen;
		r->getAngularMomentum(momen(2));
		momen += j * rad3d.cross(normal3d);
		r->setAngularMomentum(momen(2));

		//recompute the list of separating points
		std::vector<int> new_non_separating;
		for(unsigned int i=0; i<penetrating_pts.size();++i){
			double phi;
			VEC2F normal;
			int index = penetrating_pts[i];
			VEC2F vertex(vertices[index](0),vertices[index](1));
			phi = interpolate_gradient(normal, (vertex - VEC2F(0.5*_h,0.5*_h))/_h, _boundary_phi);
			normal.normalize();

			VEC2F pt_vel = r->getPointVelocity(vertex);
			if(pt_vel.dot(normal) < 1e-2){
				new_non_separating.push_back(index);
			}
		}

		penetrating_pts.clear();
		for(unsigned int i=0; i< new_non_separating.size(); ++i){
			penetrating_pts.push_back(new_non_separating[i]);
		}
	}
	// project the whole body to be outside the wall
	bool found_penetrating = false;
	do{
		found_penetrating = false;
		vertices.clear();
		r->getVertices(vertices);
		min_phi = 10*_h;
		min_ind = -1;
		min_normal = VEC2F(0,0);

		for(unsigned int i=0; i<vertices.size(); ++i){
			double phi;
			VEC2F normal;
			VEC2F vertex(vertices[i](0),vertices[i](1));
			phi = interpolate_gradient(normal,(vertex-VEC2F(0,0))/_h,_boundary_phi);
			normal.normalize();
			if(phi < min_phi){
				min_phi = phi;
				min_normal = normal;
			}
		}

		if(min_phi < 1e-3){
			found_penetrating = true;
			VEC2F com0;
			r->getCOM(com0);
			r->setCOM(com0 - (min_phi - 1e-2)*min_normal);
		}
	}while(found_penetrating);
}

void MAC_Grid::update_rigid_body_properties(const RigidBody* r){
	// used to recompute rigid body distribution on MAC grid after
	// integrating forces on rigid body
	for(int j=0; j < _rigid_phi.nY(); ++j){
		for(int i=0; i < _rigid_phi.nX(); ++i){
			_rigid_phi.update_field_val_at(i,j, r->getSignedDist(VEC2F(i*_h,j*_h)));
		}
	}
	// Use the updated distance fields to re-compute area fractions
	for(int i=0; i < _rigid_u_weight.nX(); ++i){
		for(int j=0; j < _rigid_u_weight.nY(); ++j){
			_rigid_u_weight.update_field_val_at(i,j, fraction_in(_rigid_phi(i,j), _rigid_phi(i,j+1)));
		}
	}
	for(int j=0; j < _rigid_v_weight.nY(); ++j){
		for(int i=0; i < _rigid_v_weight.nX(); ++i){
			_rigid_v_weight.update_field_val_at(i,j, fraction_in(_rigid_phi(i,j), _rigid_phi(i+1,j)));
		}
	}
	_net_rigid_u = r->getDensity()*_rigid_u_weight.field_total();
	_net_rigid_v = r->getDensity()*_rigid_v_weight.field_total();
	if(VERBOSE)
		std::cout<<"<<<< Rigid Body Net mass U,V: "<<_net_rigid_u<<", "<<_net_rigid_v<<std::endl;
}

void MAC_Grid::update_fluid_properties(const std::vector<VEC2F>& p, const double fluid_r){
	_fluid_phi.assign_field(3*_h);
	// value is estimated from nearby particles
	for(auto it=p.begin(); it!=p.end(); ++it){
		VEC2F this_pt = *it;
		// compute which cell this point belongs to
		int i,j;
		double fX, fY;
		get_barycentric(this_pt(0)/_h - 0.5, i, fX, 0, _nX);
		get_barycentric(this_pt(1)/_h - 0.5, j, fY, 0, _nY);
		// compute the distance from neighbouring points
		// and keep it if it is the minimum
		for(int j1 = j-2; j1 <= j+2; ++j1){
			for(int i1 = i-2; i1 <= i+2; ++i1){
				if(i1<0 || i1>= _nX || j1<0 || j1>=_nY)
					continue;
				VEC2F nbd_pt((i1+0.5)*_h, (j1+0.5)*_h);
				// '1.02' is a buffer value
				double phi_nbd = (nbd_pt - this_pt).norm() - 1.02*fluid_r;
				_fluid_phi.update_field_val_at(i1,j1,std::min(_fluid_phi(i1,j1),phi_nbd));
			}
		}
	}
	// extrapolating phi value to nearby boundary
	for(int j=0; j<_nY; ++j){
		for(int i=0; i<_nX; ++i){
			if(_fluid_phi(i,j) < 0.5*_h){
				double boundary_phi = 0.25*(_boundary_phi(i,j)+_boundary_phi(i+1,j)+_boundary_phi(i,j+1)+_boundary_phi(i+1,j+1));
				if(boundary_phi < 0)
					_fluid_phi.update_field_val_at(i,j, -0.5*_h);
			}
		}
	}
}

void MAC_Grid::advect_fluid_on_MAC(double dt){
	// 1st order semi-Lagrangian advection of fluid velocities
	// based on Stam's paper
	// u-component
	for(int j=0; j<_fluid_u.nY(); ++j){
		for(int i=0; i<_fluid_u.nX(); ++i){
			VEC2F this_pt(i*_h, (j+0.5)*_h);
			this_pt = rk2_tracer(this_pt, -dt);
			_fluid_u.update_field_val_at(i,j, velocity(this_pt)(0));
		}
	}
	//v-component
	for(int j=0; j<_fluid_v.nY(); ++j){
		for(int i=0; i<_fluid_v.nX(); ++i){
			VEC2F this_pt((i+0.5)*_h, j*_h);
			this_pt = rk2_tracer(this_pt, -dt);
			_fluid_v.update_field_val_at(i,j, velocity(this_pt)(1));
		}
	}
}

VEC2F MAC_Grid::rk2_tracer(const VEC2F& pos, double dt){
	// back traces the particle using RK2 time stepping
	VEC2F vel_pos = velocity(pos);
	vel_pos = velocity(pos + 0.5*dt*vel_pos);
	return (pos + dt*vel_pos);
}

void MAC_Grid::force_on_fluid_MAC(double dt){
	// for now we just add a body force: gravity
	for(int j=0; j<_fluid_v.nY(); ++j){
		for(int i=0; i<_fluid_v.nX(); ++i){
			_fluid_v.update_field_val_at(i,j, _fluid_v(i,j) - 9.81*dt);
		}
	}
}

void MAC_Grid::compute_pressure_weights(void){
	// compute FVM style face-weights for fluid using signed distances
	// account for both boundary and rigid body
	for(int j=0; j<_fluid_u_weight.nY(); ++j){
		for(int i=0; i<_fluid_u_weight.nX(); ++i){
			double val = 1 - fraction_in(_boundary_phi(i,j+1),_boundary_phi(i,j)) - _rigid_u_weight(i,j);
			// make sure value is between 0 and 1 as it is a fraction
			val = clamp(val,0,1);
			_fluid_u_weight.update_field_val_at(i,j, val);
		}
	}
	for(int j=0; j<_fluid_v_weight.nY(); ++j){
		for(int i=0; i<_fluid_v_weight.nX(); ++i){
			double val = 1 - fraction_in(_boundary_phi(i+1,j),_boundary_phi(i,j)) - _rigid_v_weight(i,j);
			// make sure value is between 0 and 1 as it is a fraction
			val = clamp(val,0,1);
			_fluid_v_weight.update_field_val_at(i,j, val);
		}
	}
}

void MAC_Grid::apply_projection(double dt, RigidBody* rbd){
	compute_pressure_weights();
	solve_pressure(dt,rbd);
}

void MAC_Grid::solve_pressure(double dt, RigidBody* rbd){
	// assemble data for the J vectors
	// translation coupling in x and y
	// rotation coupling in z
	Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> base_trans_x, base_trans_y, base_rot_z;
	base_trans_x = Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(_nX,_nY);
	base_trans_y = Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(_nX,_nY);
	base_rot_z = Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(_nX,_nY);
	VEC2F center_of_mass;
	rbd->getCOM(center_of_mass);
	for(int j=0; j<_nY; ++j){
		for(int i=0; i<_nX; ++i){
			double u_term = (_rigid_u_weight(i+1,j) - _rigid_u_weight(i,j)) / _h;
			double v_term = (_rigid_v_weight(i,j+1) - _rigid_v_weight(i,j)) / _h;
			base_trans_x(i,j) = u_term;
			base_trans_y(i,j) = v_term;
			VEC2F del = VEC2F((i+0.5)*_h,(j+0.5)*_h) - center_of_mass;
			base_rot_z(i,j) = del[0]*v_term - del[1]*u_term;
		}
	}

	int nX = _fluid_v.nX();
	int nY = _fluid_u.nY();
	int sys_size = nX*nY + 3; // pressure samples + 3 entries for rigid body (two linear vel and one rotation)

	std::vector<TRIPLET> triplets;
	VECTOR RHS(sys_size);

	bool any_liquid_surf = false;

	// building the linear system for Pressure
	// fluid part
	for(int j=1; j < nY-1; ++j){
		for(int i=1; i < nX-1; ++i){
			int idx = i + nX*j;
			RHS(idx) = 0;
			double center_phi = _fluid_phi(i,j);
			if(center_phi < 0){
				// means fluid is present, so process

				// right neighbour
				double term = _fluid_u_weight(i+1,j)*dt / (_h*_h);
				if(term > 0){
					double right_phi = _fluid_phi(i+1,j);
					if(right_phi < 0){
						triplets.push_back(TRIPLET(idx,idx,term));
						triplets.push_back(TRIPLET(idx,idx+1,-term));
					}else{
						double theta = fraction_in(center_phi, right_phi);
						if (theta < 0.01) theta = 0.01;
						triplets.push_back(TRIPLET(idx,idx,term/theta));
						any_liquid_surf = true;
					}
					RHS(idx) -= _fluid_u_weight(i+1,j)*_fluid_u(i+1,j) / _h;
				}

				// left neighbour
				term = _fluid_u_weight(i,j)*dt / (_h*_h);
				if(term > 0){
					double left_phi = _fluid_phi(i-1,j);
					if(left_phi < 0){
						triplets.push_back(TRIPLET(idx,idx,term));
						triplets.push_back(TRIPLET(idx,idx-1,-term));
					}else{
						double theta = fraction_in(center_phi, left_phi);
						if (theta < 0.01) theta = 0.01;
						triplets.push_back(TRIPLET(idx,idx,term/theta));
						any_liquid_surf = true;
					}
					RHS(idx) += _fluid_u_weight(i,j)*_fluid_u(i,j) / _h;
				}

				// top neighbour
				term = _fluid_v_weight(i,j+1)*dt / (_h*_h);
				if(term > 0){
					double top_phi = _fluid_phi(i,j+1);
					if(top_phi < 0){
						triplets.push_back(TRIPLET(idx,idx,term));
						triplets.push_back(TRIPLET(idx,idx+nX,-term));
					}else{
						double theta = fraction_in(center_phi, top_phi);
						if (theta < 0.01) theta = 0.01;
						triplets.push_back(TRIPLET(idx,idx,term/theta));
						any_liquid_surf = true;
					}
					RHS(idx) -= _fluid_v_weight(i,j+1)*_fluid_v(i,j+1) / _h;
				}

				// bottom neighbour
				term = _fluid_v_weight(i,j)*dt / (_h*_h);
				if(term > 0){
					double bottom_phi = _fluid_phi(i,j-1);
					if(bottom_phi < 0){
						triplets.push_back(TRIPLET(idx,idx,term));
						triplets.push_back(TRIPLET(idx,idx-nX,-term));
					}else{
						double theta = fraction_in(center_phi, bottom_phi);
						if (theta < 0.01) theta = 0.01;
						triplets.push_back(TRIPLET(idx,idx,term/theta));
						any_liquid_surf = true;
					}
					RHS(idx) += _fluid_v_weight(i,j)*_fluid_v(i,j) / _h;
				}

			}
		}
	}
	// fluid build completes here!

	VEC2F rigidLinVel;
	double rigidAngVel;
	rbd->getLinearVelocity(rigidLinVel);
	rbd->getAngularVelocity(rigidAngVel);

	const double Jinv = rbd->getInvInertiaModulus();

	int rigid_start_idx = nX*nY;
	// add rigid body components to system
	triplets.push_back(TRIPLET(rigid_start_idx,rigid_start_idx,-_net_rigid_u/dt));
	triplets.push_back(TRIPLET(rigid_start_idx+1,rigid_start_idx+1, -_net_rigid_v/dt));
	triplets.push_back(TRIPLET(rigid_start_idx+2,rigid_start_idx+2, -rbd->getInertiaModulus()/dt));
	RHS(rigid_start_idx) = -_net_rigid_u*rigidLinVel(0)/dt;
	RHS(rigid_start_idx+1) = -_net_rigid_v*rigidLinVel(1)/dt;
	RHS(rigid_start_idx+2) = -rbd->getInertiaModulus()*rigidAngVel/dt;

	for(int j=0; j<nY; ++j){
		for(int i=0; i<nX; ++i){
			int idx = i + nX*j;
			double center_phi = _fluid_phi(i,j);
			if(center_phi<0){
				// adding J'V to top-right block and
				// Jp to bottom left block
				double val_x = base_trans_x(i,j);
				double val_y = base_trans_y(i,j);
				double val_z = base_rot_z(i,j);
				if(std::fabs(val_x) > 1e-10){
					triplets.push_back(TRIPLET(idx,rigid_start_idx,val_x));
					triplets.push_back(TRIPLET(idx,rigid_start_idx+1,val_y));
					triplets.push_back(TRIPLET(idx,rigid_start_idx+2,val_z));
					///////////
					triplets.push_back(TRIPLET(rigid_start_idx,idx,val_x));
					triplets.push_back(TRIPLET(rigid_start_idx+1,idx,val_y));
					triplets.push_back(TRIPLET(rigid_start_idx+2,idx,val_z));
				}
			}
		}
	}

	// system building complete!

	// solve the system using Eigen Sparse solver
	SpMat LHS(sys_size,sys_size);
	LHS.setFromTriplets(triplets.begin(),triplets.end());

	// replace empty rows/cols to make +ve semi-definite
	for(int row=0; row < LHS.outerSize(); ++row){
		if(LHS.innerVector(row).nonZeros()==0){
			LHS.coeffRef(row,row) =1;
			RHS[row] =0 ;
		}
	}
	
	if(!any_liquid_surf){
		int del_index = -1;
		for(int j=0; j < nY && del_index<0;++j){
			for(int i=0;i<nX &&del_index<0;++i){
				int index = i + nX*j;
				double center_phi = _fluid_phi(i,j);
				if(center_phi<0 && (_fluid_u_weight(i+1,j)>0 || _fluid_u_weight(i,j)>0 || _fluid_v_weight(i,j+1)>0 || _fluid_v_weight(i,j)>0)){
					del_index = index;
					break;
				}
			}
		}
		if(del_index>=0){
			RHS[del_index] = 0;

			for(SpMat::InnerIterator it(LHS,del_index); it; ++it){
				if(it.col() == it.row()){
					LHS.coeffRef(it.row(), it.col()) = 1;
				}else{
					LHS.coeffRef(it.row(), it.col()) = 0;
					LHS.coeffRef(it.col(), it.row()) = 0;
				}
			}
		}
	}

	LHS.makeCompressed();
	PSolver solver;
	solver.compute(LHS);
	if (solver.info() != Eigen::Success) {
      std::cout << "Crap: Eigen factorization failed.\n";
      exit(0);
   	}
   VECTOR pressure = solver.solve(RHS);
   if (solver.info() != Eigen::Success) {
      std::cout << "Crap: Eigen solve failed.\n";
      exit(0);
   }
   // pressure solve successfully completed!

   // update velocity fields
   _fluid_u_valid.assign_field(0);
   _fluid_v_valid.assign_field(0);
   for(int j=0; j < _fluid_u.nY();++j){
   		for(int i=1; i < _fluid_u.nX()-1;++i){
   			int idx = i + nX*j;
   			if(_fluid_u_weight(i,j)>0 && ( _fluid_phi(i,j)<0 || _fluid_phi(i-1,j)<0 ) ){
   				double theta = 1;
   				if(_fluid_phi(i,j)>=0 || _fluid_phi(i-1,j)>=0)
   					theta = fraction_in(_fluid_phi(i-1,j),_fluid_phi(i,j));
   				if(theta < 0.01)
   					theta = 0.01;
   				_fluid_u.update_field_val_at(i,j,_fluid_u(i,j) - dt*( (pressure(idx) - pressure(idx-1))/_h/theta ) );
   				_fluid_u_valid.update_field_val_at(i,j,1);
   			}else
   				_fluid_u.update_field_val_at(i,j,0);
   		}
   }
   for(int j=1; j < _fluid_v.nY()-1;++j){
   		for(int i=0; i < _fluid_v.nX();++i){
   			int idx = i + nX*j;
   			if(_fluid_v_weight(i,j)>0 && ( _fluid_phi(i,j)<0 || _fluid_phi(i,j-1)<0 ) ){
   				double theta = 1;
   				if(_fluid_phi(i,j)>=0 || _fluid_phi(i,j-1)>=0)
   					theta = fraction_in(_fluid_phi(i,j-1),_fluid_phi(i,j));
   				if(theta < 0.01)
   					theta = 0.01;
   				_fluid_v.update_field_val_at(i,j,_fluid_v(i,j) - dt*( (pressure(idx) - pressure(idx-nX))/_h/theta ) );
   				_fluid_v_valid.update_field_val_at(i,j,1);
   			}else
   				_fluid_v.update_field_val_at(i,j,0);
   		}
   }

   // finally update rigid body using pressure
   VEC2F updated_rigidLinVel;
   double updated_rigidAngMomen;
   rbd->getLinearVelocity(updated_rigidLinVel);
   rbd->getAngularMomentum(updated_rigidAngMomen);
   for(int j=0; j<nY; ++j){
   		for(int i=0; i<nX; ++i){
   			int idx = i + nX*j;
   			double center_phi = _fluid_phi(i,j);
   			if(center_phi < 0){
   				updated_rigidLinVel(0) += (dt*base_trans_x(i,j)*pressure(idx)/_net_rigid_u);
   				updated_rigidLinVel(1) += (dt*base_trans_y(i,j)*pressure(idx)/_net_rigid_v);
   				updated_rigidAngMomen += (dt*base_rot_z(i,j)*pressure(idx));
   			}
   		}
   }
   rbd->setLinearVelocity(updated_rigidLinVel);
   rbd->setAngularMomentum(updated_rigidAngMomen);
}

void MAC_Grid::extrapolate_velocities(void){
	velocity_extrapolator(_fluid_u,_fluid_u_valid);
	velocity_extrapolator(_fluid_v,_fluid_v_valid);
}

void MAC_Grid::velocity_extrapolator(ScalarF& f, ScalarF& f_valid){
	// we apply several iterations of a "Jacobi"- style propogation
	// of valid field data in all directions
	ScalarF f_old_valid(f_valid.nX(),f_valid.nY(),_h);
	for(int times =0; times < 10; ++times){
		f_old_valid = f_valid;
		ScalarF f_temp = f;
		for(int j=1; j < f.nY()-1; ++j){
			for(int i=1; i <f.nX()-1; ++i){
				double sum =0;
				int count =0;
				if(!f_old_valid(i,j)){

					if(f_old_valid(i+1,j)){
						sum += f(i+1,j);
						++count;
					}

					if(f_old_valid(i-1,j)){
						sum += f(i-1,j);
						++count;
					}

					if(f_old_valid(i,j+1)){
						sum += f(i,j+1);
						++count;
					}

					if(f_old_valid(i,j-1)){
						sum += f(i,j-1);
						++count;
					}

					// for valid nbd cells, we assign this cell their avg.
					// and tag it as valid
					if(count > 0){
						f_temp.update_field_val_at(i,j,sum / (double) count);
						f_valid.update_field_val_at(i,j,1);
					}
				}
			}
		}
		f = f_temp;
	}
}

double MAC_Grid::boundary_phi(const VEC2F& pos, VEC2F& n){
	// interpolate value and also compute the normal at that point
	double val = interpolate_gradient(n, pos/ _h, _boundary_phi);
	return interpolate_value(pos/ _h, _boundary_phi);
}

void MAC_Grid::enforce_controller_bc(const std::vector<int>& ids, const double dt, const std::vector<double>& state, double a){

	VEC2F dir(std::sin(state[1]),std::cos(state[1]));
	double scaling_const = 200;
	VEC2F vel(scaling_const*std::abs(1)*dir(0)+state[0],scaling_const*std::abs(1)*dir(1)); // relative velocity
	std::set<int> ids_set(ids.begin(), ids.end());
	for(auto it = ids_set.begin(); it!= ids_set.end(); ++it){
		int i = (*it)%_nX;
		int j = (*it)/_nX;
		_fluid_v.update_field_val_at(i,j,_fluid_v(i,j)+vel(1));
		_fluid_v.update_field_val_at(i,j+1,_fluid_v(i,j+1)+vel(1));
		_fluid_u.update_field_val_at(i,j,_fluid_u(i,j)+vel(0));
		_fluid_u.update_field_val_at(i+1,j,_fluid_u(i,j+1)+vel(0));
	}
}

const ScalarF& MAC_Grid::getU(void) {return _fluid_u;}

const ScalarF& MAC_Grid::getV(void) {return _fluid_v;}

ScalarF& MAC_Grid::modifyU(void) {return _fluid_u;}

ScalarF& MAC_Grid::modifyV(void) {return _fluid_v;}

double MAC_Grid::interpolate_gradient(VEC2F& g, const VEC2F& pos, const ScalarF& grid){
	int i,j;
	double fX, fY;
	get_barycentric(pos(0), i, fX, 0, grid.nX());
	get_barycentric(pos(1), j, fY, 0, grid.nY());

	double v00 = grid(i,j);
	double v01 = grid(i,j+1);
	double v10 = grid(i+1,j);
	double v11 = grid(i+1,j+1);

	double ddy0 = (v01 - v00);
	double ddy1 = (v11 - v10);

	double ddx0 = (v10 - v00);
	double ddx1 = (v11 - v01);

	g(0) = lerp(ddx0,ddx1,fY);
	g(1) = lerp(ddy0,ddy1,fX);
	
	return bilerp(v00,v10,v01,v11,fX,fY);
}

double MAC_Grid::interpolate_value(const VEC2F& pos, const ScalarF& grid){
	int i,j;
	double fX, fY;

	get_barycentric(pos(0), i, fX, 0, grid.nX());
	get_barycentric(pos(1), j, fY, 0, grid.nY());

	return bilerp(
		grid(i,j), grid(i+1,j),
		grid(i,j+1), grid(i+1,j+1),
		fX, fY);
}

void MAC_Grid::get_barycentric(double x, int& i, double& f, int i_low, int i_high){
	double s = std::floor(x);
	i = (int)s;
	if(i < i_low){
		i = i_low;
		f = 0;
	}else if(i>i_high-2){
		i = i_high - 2;
		f = 1;
	}else
		f = x - s;
}

double MAC_Grid::bilerp(double v00, double v10, double v01, double v11, double fX, double fY) {return lerp(lerp(v00, v10, fX), lerp(v01, v11, fX), fY);}


double MAC_Grid::lerp(double v0, double v1, double f) {return ((1-f)*v0 + f*v1);}

double MAC_Grid::fraction_in(double phi_l, double phi_r){
	// compute fraction inside the surface along a 1D line joining phi_l
	// and phi_r using the signed distance field
	return (phi_l >=0 && phi_r >=0)? // means everything is empty, return '0'
			0 :
			( (phi_l < 0 && phi_r < 0)? // means all in, return '1'
				1:
				(
					(phi_l >= 0)?
					(1 - phi_l / (phi_l - phi_r)): //means right is in
					(phi_l / (phi_l - phi_r)) // means left is in
				)
			);
}

double MAC_Grid::clamp(double v, double l, double u){
	if(v<l) return l;
	else if (v>u) return u;
	else return v;
}

MATRIX3 MAC_Grid::construct_star_matrix(const VEC3F& v){
	MATRIX3 m = MATRIX3::Zero();
	m(0,1) = -v(2); m(0,2) = v(1);
	m(1,0) = v(2); m(1,2) = -v(0);
	m(2,0) = -v(1); m(2,1) = v(0);
	return m;
}
