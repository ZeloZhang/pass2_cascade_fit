/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * this class manages the double powerlaw astro physical flux model
 * 
 */

#ifndef _ASTRO_MODEL_DOUBLE_PLAW_H
#define _ASTRO_MODEL_DOUBLE_PLAW_H

#include "../helpers.h"
#include "../astro_model_base.h"

namespace NuFit
{

	class astro_model_double_plaw: public astro_model_base {
		public:
			astro_model_double_plaw();
			~astro_model_double_plaw() { };
			double get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const;
	
	};

}

#endif

