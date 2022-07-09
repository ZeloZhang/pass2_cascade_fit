/* author: Hans Niederhausen
 * date: 03/14/2017
 *
 * wrapper class to specify analysis details
 *
 */

#ifndef _ANALYSIS_H
#define _ANALYSIS_H

#include "hists.h"
#include "model_base.h"
#include "./systematics/model_base_sys.h" // if systematics are needed
#include "./systematics/interpolated_par_domeff.h"
#include "./systematics/interpolated_par_scattering.h"
#include "./systematics/interpolated_par_absorption.h"
#include "./systematics/interpolated_par_p0.h"
#include "./systematics/interpolated_par_p1.h"
#include "./systematics/interpolated_par_selfveto.h"
#include "./systematics/interpolated_par_cosmicray.h"
#include "./systematics/interpolated_par_hadronicinteraction.h"

#include "astro_model_base.h"
#include "./models/astro_model_single_plaw.h" // change if other model desired
#include "./models/astro_model_plaw_singlep.h" // change if other model desired

#include "boost/python/numeric.hpp"
#include "boost/python/extract.hpp"

namespace NuFit 
{

	class analysis
	{
		public:
			analysis();
			~analysis();
			void create();
			void change_astro_model(NuFit::astro_model_base *astro);
			void get_histograms(std::string outfile, std::map<std::string, double> &pars); // writes histograms to file
			std::vector<TH3D*> get_hist_mcsum(std::map<std::string, double> &pars); // returns pointer to mcsum histogram - relevant for toyMC class
			void set_verbosity(bool flag);

			double get_likelihood(const double *pars); // includes factor of 2 for wilk's 
			double get_likelihood_gof(const double *pars); // includes factor of 2 for wilk's	
			double get_likelihood_abs(const double *pars); // includes factor of 2 for wilk's
			double get_lnprob(boost::python::numeric::array pars); // positive logllh to be used from python	
		        void get_par_names(std::vector<std::string> &names);	
			void get_par_names(std::map<std::string, unsigned int> &names);	
                        unsigned int get_npars();
			unsigned int get_n_llh_evals();
			void reset_n_llh_evals();
			std::vector<std::string> get_analysis_names();

			void set_hist(std::string name, TH3D *hist);
			void cache_data_hists();
			void restore_data_hists();

			void reset_auxillary_data();
			void update_auxillary_data(std::map<std::string, double> &pars);	
			

	
				
		private:
			unsigned int n_llh_evals;
		        model_base *model;
			bool verbose_llh;	
	};
}

#endif
