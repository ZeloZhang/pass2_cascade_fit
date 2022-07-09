#include "../../include/models/astro_model_double_plaw.h"

NuFit::astro_model_double_plaw::astro_model_double_plaw()
{
	// model name
	model_name = std::string("double_powerlaw");

	// specify names and ordering (indices) of model parameters
	std::string par0("astro_norm");
	std::string par1("astro_index_soft");
	std::string par2("astro_index_hard");
	std::string par3("astro_fraction_hard");
	par_names.push_back(par0);
	par_names.push_back(par1);
	par_names.push_back(par2);
	par_names.push_back(par3);

	// store the ordering in a map
	for (unsigned int i=0; i<par_names.size(); ++i)
		parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));
	
	npars = parameters.size();

	return;
} 

double NuFit::astro_model_double_plaw::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above

	return pars[0] * 1.e-18 * ( (1.-pars[3]) * TMath::Power(energy / 1.e5, -1.0 * pars[1]) + pars[3] * TMath::Power(energy / 1.e5, -1.0 * pars[2]) );
}

/* sometimes it can be helpful to reparameterize index_hard = index_soft - delta_index,
 * i.e. replace index_hard by delta_index as parameter. */
