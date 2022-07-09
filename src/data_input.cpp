#include "../include/data_input.h"

NuFit::data_input::data_input(std::string name, std::vector<double> &bins_x, std::vector<double> &bins_y, std::vector<double> &bins_z)
	: hist((name+std::string("")).c_str(), (name+std::string("")).c_str(), bins_x.size()-1, &(bins_x[0]), bins_y.size()-1, &(bins_y[0]), bins_z.size()-1, &(bins_z[0])) { }
	 
	

void NuFit::data_input::read(std::string &infile) 
{
	// file infile needs to be ASCII containing the following 9 columns
	// run_id, event_id, Erec, theta_rec, azimuth_rec 

	int ncols = 5;
	 
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
			 logenergy_rec.push_back(value[2]);
		         coszenith_rec.push_back(value[3]);
			 ra_rec.push_back(value[4]);	
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

unsigned int NuFit::data_input::get_size() {
	        return size;
}
