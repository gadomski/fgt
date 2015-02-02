#pragma once

#include "parameters.hpp"


namespace ifgt
{


Parameters choose_parameters(int d, double h, double epsilon);
Parameters choose_parameters(int d, double h, double epsilon, int k_limit);


}
