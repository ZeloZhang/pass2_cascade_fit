/* author: Hans Niederhausen
 * date: 01/08/2017
 *
 * this class manages the likelihood function.
 * needs:  
 * 1) an astrophysical model (instance of astro class)
 * 2) mc input (vector of instances of hists class)
 *
 * the standard base model has 3 parameters: 
 * 	conventional norm
 * 	prompt norm
 * 	muon norm
 *
 * this class is derived from model base and adds the systematics treatment
 */

#ifndef _MODEL_BASE_SYS_H
#define _MODEL_BASE_SYS_H

#include <TROOT.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <TMath.h>

#include "../hists.h"
#include "../astro_model_base.h"
#include "../model_base.h"
#include "../helpers.h"
#include "interpolated_par.h"

#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <stdio.h>

#include <ctime>
#include <random>

namespace NuFit
{

	class model_base_sys: public model_base
	{
		public:
			model_base_sys(std::vector<hists> fitdata, astro_model_base *astro, std::map<std::string, NuFit::interpolated_par *> systematics_, unsigned int seed);
			virtual ~model_base_sys();	

			double likelihood(const double *pars); // calls likelihood from base + adds priors
			double logprior(const double *pars); // add prior penalty

			//void update_auxillary_data();
			void update_auxillary_data(std::map<std::string, double> &point);
			void reset_auxillary_data();
					
		protected:	
			void update_hists(const double *pars);
			void update_bincorrections(const double *pars);	
			void update_sum();

			double generate_normal_rv(const double &mean, const double &sigma);
			double generate_truncated_normal_rv(const double &mean, const double &sigma);
			void generate_truncated_bivnorm_rv(double means[2], double rvars[2]);
			std::mt19937 rand;

			int npars_sys;
			int npars_sys_interp; // systematics parameters that use per-bin interpolations
			int npars_sys_weights;
			std::vector<std::string> par_names_sys_weights;
			std::vector<std::string> par_names_sys_bineff;
			std::unordered_map<std::string, NuFit::interpolated_par *> systematics_interp;
			std::vector<std::string> flavors;
			std::vector<std::string> components;
			std::string skip;

                        double pivot_energy_delta_cr_conv;
                        double pivot_energy_delta_cr_prompt;
                        double pivot_energy_delta_cr_muon;

			// priors
			bool add_prior_penalty;

			double prior_domeff_mean_current;
			double prior_domeff_mean_default;
			double prior_domeff_sigma;

			double prior_deltacr_mean_current;
			double prior_deltacr_mean_default;
			double prior_deltacr_sigma;

                        double prior_muon_norm_mlb_mean_default;
                        double prior_muon_norm_mlb_mean_current;
                        double prior_muon_norm_mlb_sigma;

			double prior_scattering_current;
			double prior_scattering_default;
			double prior_absorption_current;
			double prior_absorption_default;

	};

}

#endif
