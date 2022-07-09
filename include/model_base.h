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
 * additional parameters (e.g. systematics) can be implemented in a separate class (derived from this base class)
 */

#ifndef _MODEL_BASE_H
#define _MODEL_BASE_H

#include <TROOT.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <TMath.h>

#include "hists.h"
#include "astro_model_base.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>

#include <ctime>

namespace NuFit
{

	class model_base
	{
		public:
			model_base(std::vector<hists> fitdata, NuFit::astro_model_base *astro);
			virtual ~model_base();
			
			virtual double likelihood(const double *pars); // returns negative logllh
			
			double likelihood_gof(const double *pars); // returns negative logllh according to Baker+Cousins 1984
			double likelihood_gof(double neglogl);
			double likelihood_abs(const double *pars);
			double likelihood_abs(double neglogl);
			double log_factorial(int value);
			void change_astro_model(NuFit::astro_model_base *astro);	
		
			// stuff below should not be necessary	
			void get_par_names(std::map<std::string, unsigned int> &names);
			void get_par_names(std::vector<std::string> &names);
			void get_histograms(std::string outfile, std::map<std::string, double> &pars);
			std::vector<TH3D*> get_hist_mcsum(std::map<std::string, double> &pars); // returns pointer to mcsum histogram - relevant for toyMC class
	        unsigned int get_npars_base() const;
	        unsigned int get_npars() const;
			std::vector<std::string> get_analysis_names();	

			void set_hist(const std::string &analysis, TH3D *hist);
			void cache_data_hists(); // before data hists are overwritten with toy hists. allow to cache them
            void restore_data_hists(); // after toy data fits, replace toy data with cached real exp data.

            virtual void update_auxillary_data(std::map<std::string, double> &point); // only meaningful is priors are specified
			virtual void reset_auxillary_data(); // more relevant for systematics
	
		protected:
			model_base();
			std::vector<NuFit::hists> input;
			std::map<std::string, unsigned int> input_indices; // keep track of ordering
			std::vector<std::string> par_names; // keep track of parameter names	
			std::map<std::string, unsigned int> parameters; // ... and ordering	

			std::vector<TH3D*> hist_ptrs; // convenience
	        unsigned int npars_base;
			unsigned int npars_astro;
	        unsigned int npars; // npars_base + npars_astro (depends on astro model)
	
			unsigned int ndatasets; // keep track of number of analyses added

			//std::map<std::string, TH3D*> data_hist_cache;		
			std::vector<TH3D> data_hist_cache;
	
			void store_parameters();
			void fill_parameters(std::map<std::string, double> &pars_map, double *pars);
			virtual void update_hists(const double *pars);	
			virtual void update_sum();
			NuFit::astro_model_base *astro_model;

			// cache log factorials
			std::unordered_map<int, double> log_factorials;
			void cache_logfactorials();
			double stirling_constant;
            double gaussian_prior_penalty(const double &x, const double &mean, const double &sigma);
			double bivariate_prior_penalty(const double &abs, const double &scat, const double &mean_abs, const double &mean_scat);

			// default values for priors
			//

		private:
			void update_astro(const double *astro_pars);
			void update_atmospherics(const double *pars);
            void update_sigmasq(const double *astro_pars, const double *pars);
			
	};

}

#endif
