#include "../../include/models/astro_model_padovani.h"

NuFit::astro_model_padovani::astro_model_padovani()
{


	// model name
	model_name = std::string("single_powerlaw");

	// specify names and ordering (indices) of model parameters
	std::string par0("astro_norm");
	std::string par1("astro_index");
	std::string par2("astro_Yng");
	std::string par3("astro_cutoff");
	par_names.push_back(par0);
	par_names.push_back(par1);
	par_names.push_back(par2);
	par_names.push_back(par3);

	// store the ordering in a map
	for (unsigned int i=0; i<par_names.size(); ++i)
		parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));
	
	npars = parameters.size();

	// read padovani model
        std::string infile=std::string("/data/user/hmniederhausen/analysis/ultimate/cascade_fit/final_fitcode_gitcopy/src/models/padovani_final_flux.dat");

        std::vector<double> x;
        std::vector<double> y;
        int ncols = 2;

        std::string line;
        std::fstream thisfile(infile);
        if (thisfile.is_open())
        {
	
	
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


                         x.push_back(value[0]);
                         y.push_back(value[1]);



                 }

        } 

		 
        else
        {
                std::cout << "FATAL! unable to open file: " << infile << std::endl;
                std::cout << "... exiting" << std::endl;
                exit(1);
        }

	tspline = new TSpline3("name", &(x[0]), &(y[0]), x.size());
	E0 = 1.e5;
	f0 = TMath::Power(10, tspline->Eval(TMath::Log10(E0))) * TMath::Power(E0, -2);

	return;
} 

double NuFit::astro_model_padovani::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above
	
	//return pars[0] * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[1]) + pars[2] / 0.8 * blazar(energy);
	double tf = pars[0] * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * pars[1]) * TMath::Exp(-energy / TMath::Power(10, pars[3])) + pars[2] / 0.8 * blazar(energy);
	if (ptype == NuFit::helpers::NuE)
		return 2.0 * (1. - 0.22) * tf;
	else if (ptype == NuFit::helpers::NuEBar)
		return 2.0 * 0.22 * tf;
	else
		return tf;
}

double NuFit::astro_model_padovani::blazar(double &energy) const {
	if (energy>E0) {
		//double val = tspline.Eval(TMath::Log10(energy));
		return TMath::Power(energy, -2) * TMath::Power(10, tspline->Eval(TMath::Log10(energy)));
		}
	else
		return f0 * TMath::Power(energy/E0, -1);
}
