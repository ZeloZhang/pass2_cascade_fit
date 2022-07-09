/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * container class to read and store simulated neutrino events
 *
 */

#ifndef _NEUTRINO_INPUT_H
#define _NEUTRINO_INPUT_H

#include "helpers.h"

#include <TROOT.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>

namespace NuFit
{

	class neutrino_input
	{
		public: 
			neutrino_input(std::string name, std::vector<double> &bins_x, std::vector<double> &bins_y, std::vector<double> &bins_z);
			TH3D conv; // histogram of observables. changes during fitting.
			TH3D prompt; // histogram of observables. changes during fitting.
			TH3D astro; // histogram of observables. changes during fitting.
	
			std::vector<unsigned int> run;
			std::vector<unsigned int> event;
			std::vector<double> energy_prim;
			std::vector<double> coszenith_prim;
			std::vector<double> ra_prim;
			std::vector<double> logenergy_rec;
			std::vector<double> coszenith_rec;
			std::vector<double> ra_rec;
			std::vector<double> conv_weight;
			std::vector<double> prompt_weight;
			std::vector<double> astro_weight;
			std::vector<NuFit::helpers::neutrino_type> ptype;
	
			void read(std::string infile);
			unsigned int get_size();
 			unsigned int get_nbinsx() const;
                        unsigned int get_nbinsy() const ;
                        unsigned int get_nbinsz() const;
			void clear(); // careful. this clears all std::vectors
	
		private:
			neutrino_input();
			unsigned int size;

                        // copy the bins
                        std::vector<double> binsx; // log-energy
                        std::vector<double> binsy; // cos-zenith
                        std::vector<double> binsz; // right ascension

                        unsigned int nbinsx;
                        unsigned int nbinsy;
                        unsigned int nbinsz;
	
	};

}

#endif
