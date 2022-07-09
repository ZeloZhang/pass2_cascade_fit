#include "../../include/models/astro_model_segmented.h"

NuFit::astro_model_segmented::astro_model_segmented()
{
	// model name
	model_name = std::string("segmented");

	// segmented model works in bins
	for (int i=0; i<20; ++i) 
		flux_bins.push_back(2.0 + 0.3333333334 * i);

	nbins = flux_bins.size()-1;
	flux_hist = new TH1D("flux_hist", "flux_hist", nbins, &(flux_bins[0]));
		

	// specify names and ordering (indices) of model parameters
	for (int i=0; i<nbins; ++i) {
		std::string par=std::string("astro_norm_")+std::to_string(i);
		par_names.push_back(par);
	}
	


	// store the ordering in a map
	for (unsigned int i=0; i<par_names.size(); ++i)
		parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));
	
	npars = parameters.size();


	return;
} 

double NuFit::astro_model_segmented::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above
		
	int index = flux_hist->FindBin(TMath::Log10(energy))-1;
	if (index<0 || index>=nbins)
		return 0.0;

	//std::cout << index << " " << pars[index] * 1.e-8 * TMath::Power(energy, -2.0) << std::endl;
	return pars[index] * 1.e-8 * TMath::Power(energy, -2.0);
}

/* as faster implementation would simply scale the E^-2 histograms corresponding to each bin.
 * we chose this for reasons of compatibility */
