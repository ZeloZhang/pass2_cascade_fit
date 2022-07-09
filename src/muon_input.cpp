#include "../include/muon_input.h"

NuFit::muon_input::muon_input(std::string name, std::vector<double> &bins_x, std::vector<double> &bins_y, std::vector<double> &bins_z)
	: hist((name+std::string("")).c_str(), (name+std::string("")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])) { }
	 
	

void NuFit::muon_input::read(std::string &infile) 
{
	// file infile needs to be ASCII containing the following 9 columns
	// run_id, event_id, Enu, theta_nu, azimuth_nu, Erec, theta_rec, azimuth_rec, muon_weight 

	int ncols = 9;
	 
	std::string line;	
	std::ifstream thisfile(infile);
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
                                exit(1);
                         }

		         run.push_back((unsigned int) value[0]);
			 event.push_back((unsigned int) value[1]);
			 energy_prim.push_back(value[2]);
			 coszenith_prim.push_back(value[3]);
			 ra_prim.push_back(value[4]);
			 logenergy_rec.push_back(value[5]);
		         coszenith_rec.push_back(value[6]);
			 ra_rec.push_back(value[7]);
			 muon_weight.push_back(value[8]);	
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

unsigned int NuFit::muon_input::get_size() {
	        return size;
}
