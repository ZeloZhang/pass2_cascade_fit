/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * container class to read and store simulated muon events
 *
 */

#ifndef _MUON_INPUT_H
#define _MUON_INPUT_H

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

namespace NuFit 
{

	class muon_input
	{
		public: 
			muon_input(std::string name, std::vector<double> &bins_x, std::vector<double> &bins_y, std::vector<double> &bins_z);
			TH3D hist; // histogram of observables. changes during fitting.
			TH3D hist_orig; // histogram of observables. reflects nominal input simulation.
			TH3D sigma; // histogram of observables. keep MC uncertainty information 
			TH3D efficiency_correction; // histogram of observables. efficiency correction from systematic uncertainties. 
	
			std::vector<unsigned int> run;
			std::vector<unsigned int> event;
			std::vector<double> energy_prim;
			std::vector<double> coszenith_prim;
			std::vector<double> ra_prim;
			std::vector<double> logenergy_rec;
			std::vector<double> coszenith_rec;
			std::vector<double> ra_rec;
			std::vector<double> muon_weight;
			std::vector<double> muon_weight_iter;
			std::vector<double> k;
			std::vector<double> m;
			std::vector<double> l;
	
			void read(std::string &infile);
			unsigned int get_size();
	
		private:
			muon_input();
			unsigned int size;
	
	};

}

#endif
