#include "../../include/models/astro_model_broken_plaw.h"

NuFit::astro_model_broken_plaw::astro_model_broken_plaw()
{
	// model name
	model_name = std::string("broken_powerlaw");

	// specify names and ordering (indices) of model parameters
	std::string par0("astro_norm");
	std::string par1("astro_index1");
	std::string par2("astro_index2");
	std::string par3("astro_logEbreak");
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

double NuFit::astro_model_broken_plaw::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above
	
	//double ebreak = 60.e4;
	double ebreak = TMath::Power(10.0, pars[3]);
	double norm = 0.0;
	double flux = 0.0;

	if (1.e5 < ebreak) norm = pars[0] * 1.e-18 * TMath::Power(1.e5 / ebreak, pars[1]);
	else norm = pars[0] * 1.e-18 * TMath::Power(1.e5 / ebreak, pars[2]);


	if (energy < ebreak) flux = norm * TMath::Power(energy / ebreak, -1.0 * pars[1]);
	else flux = norm * TMath::Power(energy / ebreak, -1.0 * pars[2]);

	return flux;
}
