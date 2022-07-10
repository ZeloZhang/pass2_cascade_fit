/* author: Hans Niederhausen
 * date: 11/14/2016
 *
 * this class will parse and make available the data used during fitting.
 *
 */

#ifndef _HISTS_H
#define _HISTS_H

#include <TROOT.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TMath.h>

#include "neutrino_input.h"
#include "muon_input.h"
#include "data_input.h"

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>

namespace NuFit 
{

	class hists 
	{
		public:	
		        hists(std::string name, std::vector<double> userbins_x, std::vector<double> userbins_y, std::vector<double> userbins_z);
	
			// these histograms go into the LLH template fit
			TH3D atm_conv;
			TH3D atm_conv_orig;
			TH3D atm_prompt;
			TH3D atm_prompt_orig;
			TH3D astro;
			TH3D mcsum;
			TH3D sigmasq;
			TH3D gof;
			TH3D neglogl;
	
			// these classes hold histograms per flavor. they also hold the data that has been read from file
			NuFit::neutrino_input nue;
			NuFit::neutrino_input numu;
			NuFit::neutrino_input nutau;
			NuFit::muon_input muon;
			NuFit::data_input data;
	
			std::string name;
	
			unsigned int get_nbinsx() const;
			std::vector<double> get_binsx();	
			unsigned int get_nbinsy() const;
			std::vector<double> get_binsy();
			unsigned int get_nbinsz() const;
			std::vector<double> get_binsz();

			double get_bincontent(const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz);
			void set_bincontent(const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz, const double &value);

			double get_bincontent_muon(const unsigned int &binx, const unsigned int &biny, const unsigned int &binz);
			void set_bincontent_muon(const unsigned int &binx, const unsigned int &biny, const unsigned int &binz, const double &value);
	
			void read(std::string f_nue, std::string f_numu, std::string f_nutau, std::string f_muon, std::string f_data);	
	
		private: 
			hists();
	
	                // copy the bins
	                std::vector<double> binsx; // log-energy
	                std::vector<double> binsy; // cos-zenith
	                std::vector<double> binsz; // right ascension
	
	                unsigned int nbinsx;
	                unsigned int nbinsy;
	                unsigned int nbinsz;

			/*std::map<std::string, TH3D*> nue_components;
			std::map<std::string, TH3D*> numu_components;
			std::map<std::string, TH3D*> nutau_components;
			std::map<std::string, std::map<std::string, TH3D *>> nuhists;*/
	
		};

}

#endif
