/* author: Hans Niederhausen
 * date: 03/14/2017
 *
 * wrapper class to interface ROOT's Minuit2 with analysis details and stats tasks
 *
 */

#ifndef _STATS_H
#define _STATS_H

#include "../include/analysis.h"
#include "../include/helpers.h"

#include <iostream>
#include <fstream>

// and now the ROOT Math imports

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <TROOT.h>
#include <TFile.h>

namespace NuFit 
{

	class stats
	{
		public:
			stats(NuFit::analysis &analysis_);
			~stats();
			void fit(bool werrors);
			void set_options(std::map<std::string, NuFit::helpers::par_options> opts);
			void set_tolerance(const double tol); // set ROOT's minimizer tolerance
			void get_bestpars(std::map<std::string, double> &pars);
			void change_astro_model(NuFit::astro_model_base *astro); // interface to analysis::change_astro_model
			double get_profile_llh_value(std::map<std::string, double> &point); // evaluates profile LLH at point. Parameters not specified by point will be minimized to their conditional best-fit values; Requires 	
			double get_llh_value(const std::map<std::string, double> &point);
			double get_llh_value_abs(const std::map<std::string, double> &point);
			void scan_llh(std::string outfile, std::map<std::string, NuFit::helpers::scan_options> &opts, std::map<std::string, double> &seeds); // opts specifies the grid // seeds specifies at what values the other variables should be kept;
			void scan_profile_llh(std::string outfile, std::map<std::string, NuFit::helpers::scan_options> &opts); // provides profile llh scan and stores it to file
			void scan_profile_llh_asimov(std::string outfile_scan, std::map<std::string, NuFit::helpers::scan_options> &opts, std::string outfile_asimov_hists, std::map<std::string, double> &asimov_point);
            void scan_profile_llh_asimov(std::string outfile_scan, std::map<std::string, NuFit::helpers::scan_options> &opts, std::string infile_asimov_hists);
			void set_flush_rate(unsigned int rate);

			void run_toyfits(std::string infile, std::string outfile, std::map<std::string, double> point, std::map<std::string, double> seed, unsigned int toy_id_low, unsigned int toy_id_high, bool randomize_auxilliary_data=false, std::string infile_aux = ""); // calculated profile llh fits at a given point in parameter space for all toy mc realizations
			// usually it makes sense to set the seed here to the simulated truth
                        void run_predicted_gof_bayes(std::string outfile, std::string infile, std::string infile_parameter_values, unsigned int toy_id_low, unsigned int toy_id_high);
                        double calculate_gof_toy(TFile *infile, int toy_index, std::map<std::string, double> &point, std::vector<std::string> &names);
			void get_histograms_from_samples(std::string outdir, std::map<std::string, double> simpars_fixed, std::vector<std::string> simpars_float, std::string infile_parameter_samples, unsigned int nsamples);

				
				
		private:
			ROOT::Math::Minimizer *min;
			ROOT::Math::Functor target_func;	
			NuFit::analysis &analysis; 

			std::vector<double> seeds;
			std::vector<double> stepsizes;
			std::vector<double> limits_low;
			std::vector<double> limits_high;
			std::vector<double> bestpars;	

			int flush_rate;

			bool increment_loop_state(unsigned int *counters, const unsigned int *bounds, int &loop_index, int &size);
			double get_profile_llh_value(std::map<std::string, double> &point, bool verbose, bool scan);
			void write_profile_llh(std::ofstream &file, const std::vector<std::map<std::string, double>> &points);
			void write_toy_fits(std::ofstream &file, const std::vector<std::map<std::string, double>> &points);
			void run_toyfit(TFile *infile, int toy_index, std::map<std::string, double> &point, std::vector<std::string> &names);
			void update_minimizer(const std::map<std::string, double> &pars);
			void update_seed(const std::map<std::string, double> &seed);
		        	
	};

}

#endif
