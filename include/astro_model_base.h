/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * this class manages the astrophysical model.
 * it is purely virtual. derive model and specify flux in derived class!
 */


#ifndef _ASTRO_MODEL_BASE_H
#define _ASTRO_MODEL_BASE_H

#include "helpers.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>

#include <TMath.h>

namespace NuFit
{

	class astro_model_base
	{
		public:	
			astro_model_base(); 
			virtual ~astro_model_base() { }; // this is an abstract base class. need virtual destructor
			
			virtual double get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const = 0;
			void get_par_names(std::vector<std::string> &names);
			void get_par_names(std::map<std::string, unsigned int> &map);
			unsigned int get_npars() const;
			std::string get_model_name() const;
	
		protected:
			std::string model_name;
			std::vector<std::string> par_names; // keep track of parameter names
			std::map<std::string, unsigned int> parameters;
			unsigned int npars;
	
	};

}

#endif
