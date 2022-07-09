/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * this class manages the single powerlaw astro physical flux model
 * 
 */
#ifndef _ASTRO_MODEL_PADOVANI_H
#define _ASTRO_MODEL_PADOVANI_H

#include "../helpers.h"
#include "../astro_model_base.h"
#include <TSpline.h>
//#include <vector>
//#include <string>
//#include <iostream>
//#include <fstream>
//#include <stdio.h>
//#include <sstream>

namespace NuFit
{

	class astro_model_padovani: public astro_model_base {
		public:
			astro_model_padovani();
			~astro_model_padovani() { };
			double get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const;

		//private:
			double blazar(double &energy) const;
			TSpline *tspline;
			double f0;
			double E0;
	
	};

}

#endif

