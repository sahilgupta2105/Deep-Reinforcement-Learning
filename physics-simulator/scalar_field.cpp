#include "scalar_field.h"

ScalarF::ScalarF(int nX, int nY, double h){
	_nX = nX;
	_nY = nY;
	_h = h;
	_totalCells = nX*nY;
	_data = Eigen::Array<double, Eigen::Dynamic,1>::Zero(_totalCells);
}

ScalarF::~ScalarF() {};

int ScalarF::nX(void) const {return _nX;}


int ScalarF::nY(void) const {return _nY;}

int ScalarF::size(void) const {return _totalCells;}

double ScalarF::operator()(int idxX, int idxY) const{
	assert(idxX>=0);
	assert(idxX<_nX);
	assert(idxY>=0);
	assert(idxY<_nY);
	return _data[idxX + idxY*_nX];
}

double ScalarF::operator()(int idx) const {return _data[idx];}

void ScalarF::update_field_val_at(int idxX, int idxY, double v){
	assert(idxX>=0);
	assert(idxX<_nX);
	assert(idxY>=0);
	assert(idxY<_nY);
	_data[idxX + idxY*_nX] = v;
}

double ScalarF::field_total(void){
	return _data.sum();
}

void ScalarF::assign_field(double v){
	for(int i=0; i<_totalCells;++i)
		_data(i) = v;
}

Eigen::Array<double, Eigen::Dynamic,1> ScalarF::field_data(void) const{
	return _data;
}