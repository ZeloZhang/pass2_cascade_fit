#include "../include/model_base.h"

NuFit::model_base::model_base(std::vector<hists> fitdata, NuFit::astro_model_base *astro) 
{ 
	
        input = fitdata;
        astro_model = astro;

        std::string par0("muon_norm");	
        std::string par1("conv_norm");
        std::string par2("prompt_norm");
	//std::string par3("muon_norm_mlb");
	//std::string par4("muon_norm_as");
        par_names.push_back(par0);
        par_names.push_back(par1);
        par_names.push_back(par2);
	//par_names.push_back(par3);
	//par_names.push_back(par4);

        store_parameters();

        std::cout << "... jointly analyze " << fitdata.size() << " event selections." << std::endl;
        std::cout << std::endl;
        return;
}

NuFit::model_base::~model_base() {
	delete astro_model;
	std::cout << std::endl;
	std::cout << "... end of analysis. all models cleaned from memory." << std::endl;
}

void NuFit::model_base::store_parameters()
{
	// store the ordering in a map
	parameters.clear();
        for (unsigned int i=0; i<par_names.size(); ++i)
                parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));

        npars_base = parameters.size();
        npars_astro = astro_model -> get_npars();
        npars = npars_base + npars_astro;

        // add names for astro parameters
        astro_model -> get_par_names(parameters);
        astro_model -> get_par_names(par_names);

        // have to increase parameter indices for astro parameters by npars_base
        // to ensure that base parameters are always in the front of the containers
        for (unsigned int i=npars_base; i<npars; ++i)
                parameters[par_names[i]]+=npars_base;

        // proceed with book keeping
        ndatasets = input.size();

        // keep pointers to all available histograms in this analysis
	    hist_ptrs.clear();
        for (unsigned int i=0; i<ndatasets; ++i) {

                hists &dataset = input[i];
                hist_ptrs.push_back(&dataset.nue.astro);
                hist_ptrs.push_back(&dataset.numu.astro);
                hist_ptrs.push_back(&dataset.nutau.astro);
                hist_ptrs.push_back(&dataset.nue.conv);
                hist_ptrs.push_back(&dataset.numu.conv);
                hist_ptrs.push_back(&dataset.nutau.conv);
                hist_ptrs.push_back(&dataset.nue.prompt);
                hist_ptrs.push_back(&dataset.numu.prompt);
                hist_ptrs.push_back(&dataset.nutau.prompt);
                hist_ptrs.push_back(&dataset.muon.hist);
                hist_ptrs.push_back(&dataset.astro);
                hist_ptrs.push_back(&dataset.atm_conv);
                hist_ptrs.push_back(&dataset.atm_prompt);
                hist_ptrs.push_back(&dataset.data.hist);
                hist_ptrs.push_back(&dataset.mcsum);
                hist_ptrs.push_back(&dataset.sigmasq);
                hist_ptrs.push_back(&dataset.nue.sigmasq);
                hist_ptrs.push_back(&dataset.numu.sigmasq);
                hist_ptrs.push_back(&dataset.nutau.sigmasq);
                hist_ptrs.push_back(&dataset.muon.sigmasq);

		input_indices.insert(std::pair<std::string, unsigned int>(dataset.name, i));
        }

	// and force caching of log factorials for subsequent computations
	cache_logfactorials();

	for (unsigned int i=0; i<par_names.size(); ++i)
		std::cout << i << " " << par_names[i] << " " << parameters[par_names[i]] << std::endl;

}

void NuFit::model_base::cache_logfactorials()
{
	double current_val = 0.0;
	log_factorials.insert(std::pair<int, double>(0, 0.0));
	log_factorials.insert(std::pair<int, double>(1, 0.0));
	for(int n=2; n<=100; ++n)
	{
		current_val += TMath::Log(n);
		log_factorials.insert(std::pair<int, double>(n, current_val));
	}

	stirling_constant = 0.5 * TMath::Log(2. * TMath::Pi());

	return;
}

double NuFit::model_base::log_factorial(int value)
{
	if(value<=100)
		return log_factorials[value];
	else
		return (value + 0.5) * TMath::Log(value) - value + stirling_constant;
}


void NuFit::model_base::change_astro_model(NuFit::astro_model_base *astro)
{

	std::cout << std::endl;	
	std::cout << "... remove current astro model: " << astro_model->get_model_name() << " (" << astro_model->get_npars() << " parameters)" << std::endl;
	std::cout << "... replace with new astro model: " << astro->get_model_name() << " (" << astro->get_npars() << " parameters)" << std::endl;
 		
        // change in astro model means we need to remove the prevous astro parameters first
	// from par_names and parameters containers
	for(unsigned int i=npars_base; i<npars; ++i) 
		parameters.erase(par_names[i]);

	par_names.erase(par_names.begin()+npars_base, par_names.end());

	delete astro_model;

	astro_model = astro;

	npars_astro = astro_model -> get_npars();
	npars = npars_base + npars_astro;

	// add names for astro parameters
	astro_model -> get_par_names(parameters);
	astro_model -> get_par_names(par_names);

	// have to increase parameter indices for astro parameters by npars_base
	// to ensure that base parameters are always in the front of the containers
	for (unsigned int i=npars_base; i<npars; ++i)
		parameters[par_names[i]]+=npars_base;

	std::cout << "... done" << std::endl;

        return;
}


double NuFit::model_base::likelihood(const double *pars)
{

	// calculates a joint poisson-likelihood over all bins
	// ignores terms that are independent of the model parameters

        double neglogl=0.0;

	//test if any of the parameters is NaN
	//for(unsigned int i=0; i<npars; ++i) {
	//	if(std::isnan(pars[i])) {
	//		std::cout << "!!! FATAL: parameter " << i << " is NaN " << std::endl;
	//		std::cout << "!!! parameter " << i << ":" << pars[i] << std::endl;
	//		exit( 1 );
	//	}
	//}


	// set all histogram to current parameters
	update_hists(pars);

	// need to loop over all eventselections and all bins
	for (unsigned int i=0; i<ndatasets; ++i) {
		hists &dataset = input[i];

		double observed = 0.0;
		double expected = 0.0;
		for (unsigned int k=0; k<dataset.get_nbinsx(); ++k) {
			for (unsigned int l=0; l<dataset.get_nbinsy(); ++l) {
				for (unsigned int m=0; m<dataset.get_nbinsz(); ++m) {
					
					// ROOT convention: bin indices run from 1 to N
					expected = dataset.mcsum.GetBinContent(k+1, l+1, m+1);
					observed = dataset.data.hist.GetBinContent(k+1, l+1, m+1);
					if (expected<1.e-20) // assume bin contents to be non-zero
						expected = 1.e-20 ;

					neglogl += (expected - observed * TMath::Log(expected));

				}	
			}
		}

	}


	if(std::isnan(neglogl)) {
		std::cout << " !!!!! FATAL ERROR !!!! - Likelihood evaluated to NaN. " << std::endl;
		for(unsigned int i=0; i<npars; ++i) std::cout << "parameter " << i << ":" << pars[i] << std::endl;
		exit( 1 );
	}

	return neglogl;
}

double NuFit::model_base::likelihood_gof(const double *pars) 
{

	// calculates full likelihood (including data-only dependent terms)
	double neglogl = likelihood(pars);

	neglogl = likelihood_gof(neglogl); // adds the parameter independent terms
						
	return neglogl;
}

double NuFit::model_base::likelihood_gof(double neglogl) 
{
       
        // need to loop over all eventselections and all bins
	// adds terms to likelihood that are parameter independent
	// idential to the saturated gof (Baker, Cousins 1984)
	// note that delta-llh values cancel these terms again

	//std::cout << neglogl << std::endl;

        for (unsigned int i=0; i<ndatasets; ++i) {
                hists &dataset = input[i];

                double observed = 0.0;
                for (unsigned int k=0; k<dataset.get_nbinsx(); ++k) {
                        for (unsigned int l=0; l<dataset.get_nbinsy(); ++l) {
                                for (unsigned int m=0; m<dataset.get_nbinsz(); ++m) {

                                        observed = dataset.data.hist.GetBinContent(k+1, l+1, m+1);
                                        if(observed) neglogl += (observed * TMath::Log(observed) - observed);

                                }
                        }
                }
        }

        return neglogl;
}

double NuFit::model_base::likelihood_abs(const double *pars)
{

        // calculates full likelihood (including data-only dependent terms)
        double neglogl = likelihood(pars);

        neglogl = likelihood_abs(neglogl); // adds the parameter independent terms

        return neglogl;
}

double NuFit::model_base::likelihood_abs(double neglogl)
{

        // need to loop over all eventselections and all bins
        // adds terms to likelihood that are parameter independent 
	
	//std::cout << neglogl << std::endl;

        for (unsigned int i=0; i<ndatasets; ++i) {
                hists &dataset = input[i];

                double observed = 0.0;
                for (unsigned int k=0; k<dataset.get_nbinsx(); ++k) {
                        for (unsigned int l=0; l<dataset.get_nbinsy(); ++l) {
                                for (unsigned int m=0; m<dataset.get_nbinsz(); ++m) {

                                        observed = dataset.data.hist.GetBinContent(k+1, l+1, m+1);	
                                        neglogl += log_factorial((int)observed); 

                                }
                        }
                }
        }

        return neglogl;
}


void NuFit::model_base::get_par_names(std::map<std::string, unsigned int> &names)
{
        // add model parameters to map
        for (std::map<std::string, unsigned int>::iterator it=parameters.begin(); it!=parameters.end(); ++it)
                names.insert(std::pair<std::string, unsigned int>(it->first, it->second));

        return;
}

void NuFit::model_base::get_par_names(std::vector<std::string> &names)
{
	for (unsigned int i=0; i<par_names.size(); ++i)
		names.push_back(par_names[i]);
	
	return;
}

void NuFit::model_base::update_hists(const double *pars)
{
	double astro_pars[npars_astro];
	for (unsigned int i=0; i<npars_astro; ++i) {
		//std::cout << pars[i] << std::endl;	
		astro_pars[i] = pars[npars_base+i];
	}

	update_astro(astro_pars);
	update_atmospherics(pars);	
	update_sum();
	update_sigmasq(astro_pars,pars);

}


void NuFit::model_base::update_astro(const double *astro_pars)
{
	// adjust histograms according to parameter values
	// this needs to be performed on each dataset that enters the likelihood
	for (unsigned int i=0; i<ndatasets; ++i) {
		
		hists &dataset = input[i];

		// first reset all histograms 
		dataset.nue.astro.Reset();
		dataset.numu.astro.Reset();
		dataset.nutau.astro.Reset();	

		// adjust astro histograms

		// nue 
		for (unsigned int j=0; j<dataset.nue.get_size(); ++j) 
			dataset.nue.astro.Fill(dataset.nue.logenergy_rec[j], dataset.nue.coszenith_rec[j], dataset.nue.ra_rec[j], dataset.nue.astro_weight[j]*astro_model->get_flux(astro_pars, dataset.nue.energy_prim[j], dataset.nue.coszenith_prim[j], dataset.nue.ra_prim[j], dataset.nue.ptype[j]));
			
		
		// numu
		for (unsigned int j=0; j<dataset.numu.get_size(); ++j)
                        dataset.numu.astro.Fill(dataset.numu.logenergy_rec[j], dataset.numu.coszenith_rec[j], dataset.numu.ra_rec[j], dataset.numu.astro_weight[j]*astro_model->get_flux(astro_pars, dataset.numu.energy_prim[j], dataset.numu.coszenith_prim[j], dataset.numu.ra_prim[j], dataset.numu.ptype[j]));

		// nutau
		for (unsigned int j=0; j<dataset.nutau.get_size(); ++j)
                        dataset.nutau.astro.Fill(dataset.nutau.logenergy_rec[j], dataset.nutau.coszenith_rec[j], dataset.nutau.ra_rec[j], dataset.nutau.astro_weight[j]*astro_model->get_flux(astro_pars, dataset.nutau.energy_prim[j], dataset.nutau.coszenith_prim[j], dataset.nutau.ra_prim[j], dataset.nutau.ptype[j]));

	}
	
	return;
}

void NuFit::model_base::update_sum() 
{
	for (unsigned int i=0; i<ndatasets; ++i) 
	{
		hists &dataset = input[i];		

		dataset.astro.Reset();
		dataset.astro.Add(&dataset.nue.astro);
        dataset.astro.Add(&dataset.numu.astro);
        dataset.astro.Add(&dataset.nutau.astro);

		dataset.mcsum.Reset();
		dataset.mcsum.Add(&dataset.astro);
		dataset.mcsum.Add(&dataset.atm_conv);
		dataset.mcsum.Add(&dataset.atm_prompt);
		dataset.mcsum.Add(&dataset.muon.hist);
	}
	return;
}

void NuFit::model_base::update_sigmasq(const double *astro_pars, const double *pars)
{
	// adjust histograms according to parameter values
	// this needs to be performed on each dataset that enters the likelihood
	for (unsigned int i=0; i<ndatasets; ++i) {
		
        double w;
        double w_astro;
        double w_conv;
        double w_prompt;
		hists &dataset = input[i];

		// first reset all histograms 
		dataset.nue.sigmasq.Reset();
		dataset.numu.sigmasq.Reset();
		dataset.nutau.sigmasq.Reset();	
        dataset.muon.sigmasq.Reset();

		// adjust sigmasq histograms

        // muon 

		for (unsigned int j=0; j<dataset.muon.get_size(); ++j){ 
            w = dataset.muon.muon_weight[j]*pars[0];
            dataset.muon.sigmasq.Fill(dataset.muon.logenergy_rec[j], dataset.muon.coszenith_rec[j], dataset.muon.ra_rec[j], w*w);
        }

		// nue 
		for (unsigned int j=0; j<dataset.nue.get_size(); ++j){ 
            w_astro = dataset.nue.astro_weight[j]*astro_model->get_flux(astro_pars, dataset.nue.energy_prim[j], dataset.nue.coszenith_prim[j], dataset.nue.ra_prim[j], dataset.nue.ptype[j]);
            w_conv = dataset.nue.conv_weight[j]*pars[1];
            w_prompt = dataset.nue.prompt_weight[j]*pars[2];
            w = w_astro + w_conv + w_prompt;
			dataset.nue.sigmasq.Fill(dataset.nue.logenergy_rec[j], dataset.nue.coszenith_rec[j], dataset.nue.ra_rec[j], w*w);
        }
		
		// numu
		for (unsigned int j=0; j<dataset.numu.get_size(); ++j){
            w_astro = dataset.numu.astro_weight[j]*astro_model->get_flux(astro_pars, dataset.numu.energy_prim[j], dataset.numu.coszenith_prim[j], dataset.numu.ra_prim[j], dataset.numu.ptype[j]);
            w_conv = dataset.numu.conv_weight[j]*pars[1];
            w_prompt = dataset.numu.prompt_weight[j]*pars[2];
            w = w_astro + w_conv + w_prompt;
			dataset.numu.sigmasq.Fill(dataset.numu.logenergy_rec[j], dataset.numu.coszenith_rec[j], dataset.numu.ra_rec[j], w*w);
        }

		// nutau
		for (unsigned int j=0; j<dataset.nutau.get_size(); ++j){
            w_astro = dataset.nutau.astro_weight[j]*astro_model->get_flux(astro_pars, dataset.nutau.energy_prim[j], dataset.nutau.coszenith_prim[j], dataset.nutau.ra_prim[j], dataset.nutau.ptype[j]);
            w_conv = dataset.nutau.conv_weight[j]*pars[1];
            w_prompt = dataset.nutau.prompt_weight[j]*pars[2];
            w = w_astro + w_conv + w_prompt;
			dataset.nutau.sigmasq.Fill(dataset.nutau.logenergy_rec[j], dataset.nutau.coszenith_rec[j], dataset.nutau.ra_rec[j], w*w);
        }
		dataset.sigmasq.Reset();
		dataset.sigmasq.Add(&dataset.nue.sigmasq);
		dataset.sigmasq.Add(&dataset.numu.sigmasq);
		dataset.sigmasq.Add(&dataset.nutau.sigmasq);
		dataset.sigmasq.Add(&dataset.muon.sigmasq);

	}
	
	return;
}

void NuFit::model_base::update_atmospherics(const double *pars) 
{
	for (unsigned int i=0; i<ndatasets; ++i)
	{
		hists &dataset = input[i];

		if (dataset.name.compare(std::string("cascade_mlb"))==0) {
			dataset.muon.hist.Reset();
			dataset.muon.hist = *((TH3D*) dataset.muon.hist_orig.Clone());
			dataset.muon.hist.Scale(pars[3]); // muon norm (no shape change)	
		}

		else if (dataset.name.compare(std::string("cascade_as"))==0) {
		        dataset.muon.hist.Reset();
                        dataset.muon.hist = *((TH3D*) dataset.muon.hist_orig.Clone());
                        dataset.muon.hist.Scale(pars[4]); // muon norm (no shape change)	
		}
		else {
                        dataset.muon.hist.Reset();
                        dataset.muon.hist = *((TH3D*) dataset.muon.hist_orig.Clone());
                        dataset.muon.hist.Scale(pars[0]); // muon norm (no shape change)
		}

		dataset.atm_conv.Reset();
		dataset.atm_conv = *((TH3D*) dataset.atm_conv_orig.Clone());
		dataset.atm_conv.Scale(pars[1]); // conv norm (no shape change)

		dataset.atm_prompt.Reset();
		dataset.atm_prompt = *((TH3D*) dataset.atm_prompt_orig.Clone());
		dataset.atm_prompt.Scale(pars[2]); // prompt norm (no shape change)	
	}
	return;
}



unsigned int NuFit::model_base::get_npars_base() const {
	return npars_base;
}

unsigned int NuFit::model_base::get_npars() const {
	return npars;
}


void NuFit::model_base::get_histograms(std::string outfile, std::map<std::string, double> &pars_map) 
{

	// order parameters as required by model. check for parameter names
        double pars[npars];
	fill_parameters(pars_map, pars);

	// ** this is to time speed of histogram filling
	//std::clock_t start;
	//double duration;
	
	TFile fout(outfile.c_str(), "RECREATE");

	// ** initilize histograms to requested hypothesis
	// ** and time evaluation
	
	//start = std::clock();

	update_hists(pars);

	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//std::cout<<"one cycle through data: "<< duration <<std::endl;;

	std::cout << std::endl;
	std::cout << "... writing analysis histograms" << std::endl;

	TH3D* clone;
	for (unsigned int i=0; i<hist_ptrs.size(); ++i) {
		clone = (TH3D*) hist_ptrs[i] -> Clone();
		clone -> Write();
		delete clone;
	}

	fout.Close();

	std::cout << "... done" << std::endl;
	std::cout << std::endl;

}	

std::vector<TH3D*> NuFit::model_base::get_hist_mcsum(std::map<std::string, double> &pars_map)
{
	double pars[npars];
	fill_parameters(pars_map, pars);
	update_hists(pars);

	std::vector<TH3D*> hists_mcsum;
	// now copy histograms for each event selection (dataset)
	for (unsigned int i=0; i<input.size(); ++i)
	{
		hists &dataset = input[i];
		TH3D* hist = (TH3D*) dataset.mcsum.Clone();
		hists_mcsum.push_back(hist);
	}

	return hists_mcsum;
}

void NuFit::model_base::fill_parameters(std::map<std::string, double> &pars_map, double *pars)
{
	// check if we have received the expected number of parameter values
        if(pars_map.size()!=npars) {
                std::cout << "!!! FATAL: received unexpected number of parameters: " << pars_map.size() << " instead of: " << npars << std::endl;
                std::cout << "... exiting" << std::endl;
                exit(1);
        }

        for (std::map<std::string, unsigned int>::iterator it = parameters.begin(); it != parameters.end(); ++it )
        {
                // does parameter exist?
                if(pars_map.count(it->first)==0) {
                        std::cout << "!!! FATAL: parameter "<< it->first << " not found. Did you specify it correctly?" << std::endl;
                        std::cout << "!!! only found the following parameters:" << std::endl;
                        for ( std::map<std::string, double>::iterator it2 = pars_map.begin(); it2 != pars_map.end(); ++it2 )
                                std::cout << it2->first << std::endl;
                        std::cout << "... exiting" << std::endl;
                        exit(1);
                }
                // parameter found
                pars[it->second]=pars_map.at(it->first);
        }

	return;
}

std::vector<std::string> NuFit::model_base::get_analysis_names()
{

	unsigned int size = input.size();
	std::vector<std::string> names;
	names.reserve(size);
	for (unsigned int i=0; i<size; ++i)
	{
		names.push_back(input[i].name);
	}

	return names;
}

void NuFit::model_base::cache_data_hists()
{
	data_hist_cache.clear();
	for (unsigned int i=0; i<input.size(); ++i)
		data_hist_cache.push_back(input[i].data.hist);
	
	return;
}

void NuFit::model_base::restore_data_hists()
{
	for (unsigned int i=0; i<input.size(); ++i)
		input[i].data.hist = data_hist_cache[i];


	data_hist_cache.clear();
	return;
}


void NuFit::model_base::set_hist(const std::string &analysis, TH3D *hist)
{
	unsigned int index = input_indices[analysis];
	input[index].data.hist = *hist;
	return;
}

double NuFit::model_base::gaussian_prior_penalty(const double &x, const double &mean, const double &sigma)
{
        return 1./(2.0 * sigma * sigma) * (x - mean) * (x - mean);
}

double NuFit::model_base::bivariate_prior_penalty(const double &abs, const double &scat, const double &mean_abs, const double &mean_scat)
{

        // need inverse covariance matrix. inverse is symmetric.
        // corresponding to sigma = 0.07, rho = -0.1

        double cov_inf11 = 206.143;
        double cov_inf12 = 20.614;

        //double delta_1 = abs-1.0;
        //double delta_2 = scat-1.0;

        double delta_1 = abs - mean_abs;
        double delta_2 = scat - mean_scat;

        return 0.5*(delta_1 * cov_inf11 * delta_1 + delta_2 * cov_inf11 * delta_2 + 2.0 * (delta_1 * cov_inf12 * delta_2));
}

void NuFit::model_base::update_auxillary_data(std::map<std::string, double> &point)
{
        // use only during toy fits
        // make sure means are consistent with what parameter values have been simulated as toy
        return;

}

void NuFit::model_base::reset_auxillary_data()
{
        return;
}
