#include "../include/astro_model_base.h"

NuFit::astro_model_base::astro_model_base() { } 

void NuFit::astro_model_base::get_par_names(std::map<std::string, unsigned int> &map)
{
	// add model parameters to map
	for (std::map<std::string, unsigned int>::iterator it=parameters.begin(); it!=parameters.end(); ++it)
		map.insert(std::pair<std::string, unsigned int>(it->first, it->second));

	return;
}

void NuFit::astro_model_base::get_par_names(std::vector<std::string> &names)
{
        // add parameter names to vector
        for (unsigned int i=0; i<par_names.size(); ++i)
		names.push_back(par_names[i]);
           
        return;
}

unsigned int NuFit::astro_model_base::get_npars() const
{
	return npars;
}

std::string NuFit::astro_model_base::get_model_name() const
{
	return model_name;
}
