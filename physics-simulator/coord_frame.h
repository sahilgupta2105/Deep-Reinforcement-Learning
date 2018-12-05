#ifndef COORD_FRAME_H
#define COORD_FRAME_H

#include "eigen_def.h"
#include <cmath>

class CoordFrame{
private:
	VEC2F _pos;
	double _orient;
	VEC2F rotate_matrix(const VEC2F& pt, double a) const{
		double c = std::cos(a);
		double s = std::sin(a);
		return VEC2F(c*pt(0) - s*pt(1), s*pt(0) + c*pt(1));
	}
public:
	CoordFrame(): _pos(0,0), _orient(0){}
	CoordFrame(VEC2F p, double a): _pos(p), _orient(a){}

	void set_position(const VEC2F& p){
		_pos = p;
	}

	void set_orientation(double a){
		_orient = a;
	}

	VEC2F get_position(void) const{
		return _pos;
	}

	double get_orientation(void){
		return _orient;
	}

	void translate(const VEC2F& v){
		_pos += v;
	}

	void rotate(float a){
		_orient += a;
	}

	VEC2F global_to_local(const VEC2F& pt) const{
		return rotate_matrix(pt - _pos,-_orient);
	}

	VEC2F local_to_global(const VEC2F& pt) const{
		return (rotate_matrix(pt ,_orient)+_pos);
	}

	void operator=(CoordFrame f){
		_pos = f.get_position();
		_orient = f.get_orientation();
	}
};

#endif