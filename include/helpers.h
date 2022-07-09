/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * various helper classes 
 *
 */


#ifndef _HELPERS_H
#define _HELPERS_H

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>

namespace NuFit 
{

	namespace helpers
	{
		class par_options
		{
			public:
				par_options();
				par_options(std::string name_, double seed_, double stepsize_, double limit_low_, double limit_high_);
				void set_options(std::string name_, double seed_, double stepsize_, double limit_low_, double limit_high_);
				std::string name;
				double seed;
				double stepsize;
				double limit_low;
				double limit_high;
		};

                class scan_options
                {
                        public:
                                scan_options();
                                scan_options(std::string name_, unsigned int nsteps_, double range_low_, double range_high_);
                                void set_options(std::string name_, unsigned int nsteps_, double range_low_, double range_high_);
                                std::string name; 
                                unsigned int nsteps;
                                double range_low;
                                double range_high;
                };

		enum neutrino_type {
			// uses pdg encoding
			NuE = 12,
			NuEBar = -12,
			NuMu = 14,
			NuMuBar = -14,
			NuTau = 16,
			NuTauBar = -16

			// if used by systematics only NuE, NuMu, NuTau are used
			// but astro fluxes can also depend on Nu*Bar
		};

	}

}

#endif
