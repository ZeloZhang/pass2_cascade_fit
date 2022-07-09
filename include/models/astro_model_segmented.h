/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * this class manages the single powerlaw astro physical flux model
 * 
 */

#ifndef _ASTRO_MODEL_SEGMENTED_H
#define _ASTRO_MODEL_SEGMENTED_H

#include "../helpers.h"
#include "../astro_model_base.h"

#include <TH1.h>
#include <TROOT.h>
#include <TMath.h>

namespace NuFit
{

	class astro_model_segmented: public astro_model_base {
		public:
			astro_model_segmented();
			~astro_model_segmented() { flux_hist -> Delete(); };
			double get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const;

		private:
			std::vector<double> flux_bins;
			int nbins;
			TH1D* flux_hist;
	
	};

}

#endif

