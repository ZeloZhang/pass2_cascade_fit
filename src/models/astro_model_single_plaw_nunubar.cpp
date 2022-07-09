#include "../../include/models/astro_model_single_plaw_nunubar.h"

NuFit::astro_model_single_plaw_nunubar::astro_model_single_plaw_nunubar()
{
	// model name
	model_name = std::string("single_powerlaw");

	// specify names and ordering (indices) of model parameters
	std::string par0("astro_norm");
	std::string par1("astro_index");
	par_names.push_back(par0);
	par_names.push_back(par1);


	// store the ordering in a map
	for (unsigned int i=0; i<par_names.size(); ++i)
		parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));
	
	npars = parameters.size();

	return;
} 

double NuFit::astro_model_single_plaw_nunubar::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above
	
	double tf = pars[0] * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[1]);
	if (ptype == NuFit::helpers::NuE) 
		return 2.0 * (1. - 0.22) * tf;
	else if (ptype == NuFit::helpers::NuEBar)
		return 2.0 * 0.22 * tf;
	else if (ptype == NuFit::helpers::NuMu)
		return 2.0 * 0.61 * tf;
	else if (ptype == NuFit::helpers::NuMuBar)
		return 2.0 * 0.39* tf;
        else if (ptype == NuFit::helpers::NuTau)
                return 2.0 * 0.61 * tf;
        else if (ptype == NuFit::helpers::NuTauBar)
                return 2.0 * 0.39 * tf;
	else
		return 0.;

		
}
