/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * class to calculate per-bin correction factors via interpolation of discrete simulations
 *
 */

#ifndef _INTERPOLATED_PAR_COSMICRAY_H
#define _INTERPOLATED_PAR_COSMICRAY_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <map>

#include "interpolated_par.h"
#include "../helpers.h"
#include "../hists.h"

namespace NuFit
{
	class interpolated_par_cosmicray: public interpolated_par
	{
		public:
			interpolated_par_cosmicray(std::string par_name_, std::vector<std::string> analysis_names_, const std::map<std::string, NuFit::hists*> &map_analyses);
			~interpolated_par_cosmicray();		
	};
}

#endif


