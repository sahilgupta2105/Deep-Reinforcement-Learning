#ifndef SCALAR_FIELD_H
#define SCALAR_FIELD_H

#include <vector>
#include "eigen_def.h"

class ScalarF{
private:
	int _nX, _nY, _totalCells;
	double _h;
	Eigen::Array<double, Eigen::Dynamic, 1> _data;
public:
	ScalarF(int nX, int nY, double h);
	~ScalarF();
	int nX(void) const;
	int nY(void) const;
	double operator()(int, int ) const;
	double operator()(int) const;
	void update_field_val_at(int, int, double);
	double field_total(void);
	void assign_field(double);
	Eigen::Array<double, Eigen::Dynamic, 1> field_data(void) const;
	int size(void) const;
};

#endif