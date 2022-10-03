/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * class to calculate per-bin correction factors via interpolation of discrete simulations
 *
 */

#ifndef _INTERPOLATED_PAR_H
#define _INTERPOLATED_PAR_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <unordered_map>
#include <algorithm>

#include "interpolated_par.h"
#include "interpolated_sys.h"
#include "../helpers.h"


namespace NuFit
{
	class interpolated_par
	{
		public:
			interpolated_par(std::string par_name_, std::vector<std::string> analysis_names_);
			virtual ~interpolated_par();
			double get_efficiency_correction(const double &x, const std::string &analysis_name, const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz);
			double get_efficiency_correction_error(const double &x, const std::string &analysis_name, const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz);
			double get_efficiency_correction_muon(const double &x, const std::string &analysis_name, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz);

		protected:
			std::string par_name;	
			std::vector<std::string> analysis_names;
			std::unordered_map<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>> binfits;
			// first keys are strings: analyses
			// second keys are flavors: NuE, NuMu, NuTau	

		        std::vector<std::string> components; // stores what components are adjusted	

			std::string basedir;
			std::string basedir_mlb;
	};
}

#endif


