#include "../../include/models/astro_model_plaw_singlep.h"

NuFit::astro_model_plaw_singlep::astro_model_plaw_singlep()
{
	// model name
	model_name = std::string("powerlawsinglep");

	// specify names and ordering (indices) of model parameters
	std::string par0("astro_index");
	par_names.push_back(par0);

	// store the ordering in a map
	for (unsigned int i=0; i<par_names.size(); ++i)
		parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));
	
	npars = parameters.size();

	return;
} 

double NuFit::astro_model_plaw_singlep::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above
	
	//std::cout << " par0=" << pars[0] << " flux=" << 1.7 * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[0])  << std::endl;
	return 1.7 * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[0]);
}
