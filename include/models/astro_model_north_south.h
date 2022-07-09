/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * this class manages the single powerlaw astro physical flux model
 * 
 */

#ifndef _ASTRO_MODEL_NORTH_SOUTH_H
#define _ASTRO_MODEL_NORTH_SOUTH_H

#include "../helpers.h"
#include "../astro_model_base.h"

namespace NuFit
{

	class astro_model_north_south: public astro_model_base {
		public:
			astro_model_north_south();
			~astro_model_north_south() { };
			double get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const;
	
	};

}

#endif

