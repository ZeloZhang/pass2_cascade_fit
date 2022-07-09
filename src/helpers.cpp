#include "../include/helpers.h"

NuFit::helpers::par_options::par_options() 
{

	// initialize to some default
	name = std::string("name");
	seed = 0.0;
	stepsize = 0.01;
	limit_low = -1.e5;
	limit_high = 1.e5;
	return;

}

NuFit::helpers::par_options::par_options(std::string name_, double seed_, double stepsize_, double limit_low_, double limit_high_)
{
	// initialize to arguments
        name = name_;
        seed = seed_;
        stepsize = stepsize_;
        limit_low = limit_low_;
        limit_high = limit_high_; 
	return;

}

void NuFit::helpers::par_options::set_options(std::string name_, double seed_, double stepsize_, double limit_low_, double limit_high_)
{
	// set parameters anytime
	name = name_;
	seed = seed_;
	stepsize = stepsize_;
	limit_low = limit_low_;
	limit_high = limit_high_;
	return;
}

NuFit::helpers::scan_options::scan_options()
{

        // initialize to some default
        name = std::string("name"); 
        nsteps = 10;
        range_low = -1.e5;
        range_high = 1.e5;
        return;

}

NuFit::helpers::scan_options::scan_options(std::string name_, unsigned int nsteps_, double range_low_, double range_high_)
{
        // initialize to arguments
        name = name_; 
        nsteps = nsteps_;
        range_low = range_low_;
        range_high = range_high_;
        return;

}

void NuFit::helpers::scan_options::set_options(std::string name_, unsigned int nsteps_, double range_low_, double range_high_)
{
        // set parameters anytime
        name = name_; 
        nsteps = nsteps_;
        range_low = range_low_;
        range_high = range_high_;
        return;
}



	 
	


