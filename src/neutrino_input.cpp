#include "../include/neutrino_input.h"

NuFit::neutrino_input::neutrino_input(std::string name, std::vector<double> &bins_x, std::vector<double> &bins_y, std::vector<double> &bins_z)
	: conv((name+std::string("_conv")).c_str(), (name+std::string("_conv")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])),
	  prompt((name+std::string("_prompt")).c_str(), (name+std::string("_prompt")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])),
	  astro((name+std::string("_astro")).c_str(), (name+std::string("_astro")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])), 
	  sigma((name+std::string("_sigma")).c_str(), (name+std::string("_sigma")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])),
	  conv_efficiency_correction((name+std::string("_conv_efficiency_correction")).c_str(), (name+std::string("_conv_efficiency_correction")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])),
	  prompt_efficiency_correction((name+std::string("_prompt_efficiency_correction")).c_str(), (name+std::string("_prompt_efficiency_correction")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])),
	  astro_efficiency_correction((name+std::string("_astro_efficiency_correction")).c_str(), (name+std::string("_astro_efficiency_correction")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0]))
{ 
    binsx = bins_x;
    nbinsx = bins_x.size()-1;

    binsy = bins_y;
    nbinsy = bins_y.size()-1;

    binsz = bins_z;
    nbinsz = bins_z.size()-1;

    for (unsigned int k=0; k<nbinsx; ++k)
    {
        for (unsigned int l=0; k<nbinsy; ++k)
        {
            for (unsigned int m=0; k<nbinsz; ++k)
            {
                conv_efficiency_correction.SetBinContent(k+1,l+1,m+1,1);
                prompt_efficiency_correction.SetBinContent(k+1,l+1,m+1,1);
                astro_efficiency_correction.SetBinContent(k+1,l+1,m+1,1);
            }
        }
    }
}

void NuFit::neutrino_input::read(std::string infile) 
{
	// file infile needs to be ASCII containing the following 11 columns
	// run_id, event_id, Enu, theta_nu, azimuth_nu, Erec, theta_rec, azimuth_rec, conv_weight, prompt_weight, astro_weight 

	int ncols = 12;
	 
	std::string line;	
	std::ifstream thisfile(infile);
	//if(!thisfile) std::cout << "FATAL! file not found." << std::endl;
	if (thisfile.is_open()) 
	{
                 // first line contains variable headers
                 getline(thisfile, line);
		 while(getline(thisfile, line))
		 {
			 std::stringstream ss(line);
			 double value[ncols];
			 int col=0;
			 while (ss >> value[col]) col++;  	

			 if( col!=ncols ) {
			 	std::cout << "FATAL! unexpected number of columns!" << std::endl;
				std::cout << "... exiting" << std::endl;

                         	//for(unsigned int i=0; i<ncols; ++i)
				//                    std::cout << value[i] << " ";
				exit(1);
			 }

		     run.push_back((unsigned int) value[0]);
			 event.push_back((unsigned int) value[1]);
			 ptype.push_back((NuFit::helpers::neutrino_type) value[2]);
			 energy_prim.push_back(value[3]);
			 coszenith_prim.push_back(value[4]);
			 ra_prim.push_back(value[5]);
			 logenergy_rec.push_back(value[6]);
		     coszenith_rec.push_back(value[7]);
			 ra_rec.push_back(value[8]);
			 conv_weight.push_back(value[9]);
			 conv_weight_iter.push_back(0);
			 prompt_weight.push_back(value[10]);
			 prompt_weight_iter.push_back(0);
			 astro_weight.push_back(value[11]);	
			 astro_weight_iter.push_back(0);	
             m.push_back(sigma.GetXaxis()->FindBin(value[8]));
             l.push_back(sigma.GetXaxis()->FindBin(value[7]));
             k.push_back(sigma.GetXaxis()->FindBin(value[6]));
		 }
	}
	else 
	{
		std::cout << "FATAL! unable to open file: " << infile << std::endl; 
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	size = run.size();
	return;
}

unsigned int NuFit::neutrino_input::get_size() {
	return size;
}

unsigned int NuFit::neutrino_input::get_nbinsx() const {
        return nbinsx;
}

unsigned int NuFit::neutrino_input::get_nbinsy() const {
        return nbinsy;
}

unsigned int NuFit::neutrino_input::get_nbinsz() const {
        return nbinsz;
}

void NuFit::neutrino_input::clear()
{
	run.clear();
	event.clear();
	energy_prim.clear();
	coszenith_prim.clear();
	ra_prim.clear();
	logenergy_rec.clear();
	coszenith_rec.clear();
	ra_rec.clear();
	conv_weight.clear();
	prompt_weight.clear();
	astro_weight.clear();
	ptype.clear();

	return;

}
