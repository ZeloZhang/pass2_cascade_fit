/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * container class to read histograms and perform parametric bootstrap re-simulation
 *
 */


#ifndef _TOY_MC_H
#define _TOY_MC_H

#include <TROOT.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TMath.h>

#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <random>

#include "../analysis.h"

namespace NuFit 
{

	class toymc 
	{
		public: 
			toymc(NuFit::analysis &analysis_, double seed);	
			void run_simulation(std::string outfile, std::map<std::string, double> simpars_, unsigned int nsamples);
			void run_simulation_waux(std::string outfile, std::string outfile_waux,std::map<std::string, double> simpars_, unsigned int nsamples);
			void run_simulation_bayesian(std::string outfile, std::map<std::string, double> simpars_fixed, std::vector<std::string> simpars_float, std::string infile_parameter_samples, unsigned int nsamples);
			

		private:	
			std::vector<TH3D *> draw_sample(std::map<std::string, double> &pars);
			NuFit::analysis &analysis;
			std::mt19937 rand;
                        double generate_normal_rv(const double &mean, const double &sigma);
                        double generate_truncated_normal_rv(const double &mean, const double &sigma);
                        void generate_truncated_bivnorm_rv(double means[2], double rvars[2]);
	};

}

#endif
