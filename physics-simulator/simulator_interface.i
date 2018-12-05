%module physics_sim

%include "std_vector.i"

namespace std {
    %template(VectorDouble) vector<double>;
};

%{
#include "scalar_field.h"
#include "eigen_def.h"
#include "circle_geometry.h"
#include "coord_frame.h"
#include "rigid_body.h"
#include "MAC.h"
#include "poisson_disk_sampling.h"
#include "simulator.h"
%}

%include "scalar_field.h"
%include "eigen_def.h"
%include "circle_geometry.h"
%include "rigid_body.h"
%include "MAC.h"
%include "poisson_disk_sampling.h"
%include "simulator.h"