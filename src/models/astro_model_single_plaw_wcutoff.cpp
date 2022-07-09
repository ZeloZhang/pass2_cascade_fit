#include "../../include/models/astro_model_single_plaw_wcutoff.h"

NuFit::astro_model_single_plaw_wcutoff::astro_model_single_plaw_wcutoff()
{
	// model name
	model_name = std::string("single_powerlaw_wcutoff");

	// specify names and ordering (indices) of model parameters
	std::string par0("astro_norm");
	std::string par1("astro_index");
	std::string par2("astro_cutoff"); // this is Ecutoff in PeV
	par_names.push_back(par0);
	par_names.push_back(par1);
	par_names.push_back(par2);

	// store the ordering in a map
	for (unsigned int i=0; i<par_names.size(); ++i)
		parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));
	
	npars = parameters.size();

	return;
} 

double NuFit::astro_model_single_plaw_wcutoff::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above
	
	//std::cout << " par0=" << pars[0] << " pars1=" << pars[1] << " flux=" << pars[0] * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[1])  << std::endl;
	//return pars[0] * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[1]) * TMath::Exp(-1.0/1.e6 * energy/pars[2]) * TMath::Exp(0.1/pars[2]);
        double tf = pars[0] * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[1]) * TMath::Exp(-energy / TMath::Power(10, pars[2]));
        //if (ptype == NuFit::helpers::NuE)
        //        return 2.0 * (1. - 0.22) * tf;
        //else if (ptype == NuFit::helpers::NuEBar)
        //        return 2.0 * 0.5 * tf;
        //else
        return tf;
}
