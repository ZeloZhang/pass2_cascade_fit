/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * class to calculate per-bin correction factors via interpolation of discrete simulations
 *
 */

#ifndef _INTERPOLATED_SYS_H
#define _INTERPOLATED_SYS_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <unordered_map>

#include "../neutrino_input.h"

// ROOT stuff
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include "Math/Interpolator.h"
namespace NuFit
{
    class root_interpolator_wrapper
    {
        public:
            ROOT::Math::Interpolator *itp;
            root_interpolator_wrapper(ROOT::Math::Interpolator *itp1){
               itp = itp1;
            }
            double Eval(double *x, double *p) {
                ////////std::cout<<"riw::Eval()"<<itp->Eval(500)<<std::endl;
                return itp->Eval(x[0]);
            }
    };

	class interpolated_sys
	{
		public:
			interpolated_sys(std::string parameter_name, std::string analysis_name, std::string nu_flavor, std::vector<double> bins_x, std::vector<double> bins_y, std::vector<double> bins_z);	
			void add_simulated_point(double value, std::string infile, bool is_baseline=false); // be sure to add points in ascending order!
			void create_correction_functions(bool use_interpolation=false);
			//void create_correction_functions();
			double get_efficiency_correction(const double &x, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz);
			double get_relative_correction_error(const double &x, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz);

		private:	
			void fill_hists();
			void create_fit(const double &baseline_val, const double &baseline_err, std::vector<TH3D *> &hists, const unsigned int &k, const unsigned int &l, const unsigned int &m, std::vector<std::vector<std::vector<bool>>> &success, std::string component, bool use_interpolation=false);
			std::vector<double> binsx; // log E
			std::vector<double> binsy; // cos zenith
			std::vector<double> binsz; // ra

			
		
			std::unordered_map<std::string, std::vector<std::vector<std::vector<TF1>>>> correction_functions; // for faster access to fit functions (by selection for flavor)
			std::unordered_map<std::string, std::vector<std::vector<std::vector<TFitResultPtr>>>> correction_fit_results; // for faster access to correction fit results and confidence intervals (by selection for flavor)
			std::unordered_map<std::string, std::vector<std::vector<std::vector<double>>>> correction_errors; // for faster access to correction error (by selection for flavor)
			// total number of nested vector elements are NbinsX * NbinsY & NbinsZ
			std::vector<std::vector<std::vector<bool>>> success_conv; // len is number of bins
			std::vector<std::vector<std::vector<bool>>> success_prompt; // to keep track of fit status in each bin
			std::vector<std::vector<std::vector<bool>>> success_astro;

			std::vector<double> values; // len is number of datasets
			std::vector<bool> baseline; // len is number of datasets
			double baseline_value;
		
			std::vector<NuFit::neutrino_input> sim_data;
			std::string par_name; // name of parameter
			std::string name; // name of event selection that is being parametrized			
			std::string flavor;

			bool verbose;

	};
}

#endif


