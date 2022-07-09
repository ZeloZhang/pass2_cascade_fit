/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * container class to read and maintain exp. data events
 *
 */


#ifndef _DATA_INPUT_H
#define _DATA_INPUT_H

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

	class data_input
	{
		public: 
			data_input(std::string name, std::vector<double> &bins_x, std::vector<double> &bins_y, std::vector<double> &bins_z);
			TH3D hist; // data histogram. fixed.
	
			std::vector<unsigned int> run;
			std::vector<unsigned int> event;
			std::vector<double> logenergy_rec;
			std::vector<double> coszenith_rec;
			std::vector<double> ra_rec;
	
			void read(std::string &infile);
			unsigned int get_size();
	
		private:
			data_input();
			unsigned int size;
	
	};

}

#endif
