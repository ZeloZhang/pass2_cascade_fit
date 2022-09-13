#include "../../include/bootstrap/toymc.h"

NuFit::toymc::toymc(NuFit::analysis &analysis_, double seed) : analysis(analysis_), rand(seed) { }

std::vector<TH3D *> NuFit::toymc::draw_sample(std::map<std::string, double> &pars) 
{
	std::vector<TH3D *> hists_mcsum = analysis.get_hist_mcsum(pars);
	std::vector<TH3D *> hists_mcsum_orig;
	for (unsigned int i=0; i<hists_mcsum.size(); ++i)
	{
		TH3D *hist = hists_mcsum[i];
        // store original mcsum
        TH3D *hist_clone = (TH3D*) hist -> Clone();
        std::string name = std::string(hist_clone->GetName(), 0, 1000);
        hist_clone->SetName((std::string(name)+std::string("_")+std::string("orig")).c_str());
        hists_mcsum_orig.push_back(hist_clone);
		for (int k=0; k<hist->GetNbinsX(); ++k)
		{
			for (int l=0; l<hist->GetNbinsY(); ++l)
			{
				for (int m=0; m<hist->GetNbinsZ(); ++m)
				{
					double mean = hist->GetBinContent(k+1, l+1, m+1);
					std::poisson_distribution<int> distribution(mean);
					int data = distribution(rand);	
					hist->SetBinContent(k+1, l+1, m+1, (double)data);
					hist->SetBinError(k+1, l+1, m+1, TMath::Sqrt(data));
				}
			}
		}
	}
    std::vector<TH3D *> hists;
    for(unsigned int i=0; i<hists_mcsum.size(); ++i)
    {
        hists.push_back(hists_mcsum[i]);
        hists.push_back(hists_mcsum_orig[i]);
	}
	return hists;
}

void NuFit::toymc::run_simulation(std::string outfile, std::map<std::string, double> simpars, unsigned int nsamples)
{
	std::cout << "\n";
	std::cout << "... sampling " << nsamples << " toymc datasets for fixed set of model parameters." << std::endl;
	TFile fout(outfile.c_str(), "RECREATE");

	for (unsigned int i=0; i<nsamples; ++i) 
	{
		std::vector<TH3D*> hists = draw_sample(simpars);
		for (unsigned int j=0; j<hists.size(); ++j)
		{
			// ROOT uses C-style character arrays
			TH3D *hist = hists[j];
			std::string name = std::string(hist->GetName(), 0, 1000);
			hist->SetName((std::string(name)+std::string("_")+std::to_string(i)).c_str());
			hist->Write();
			delete hists[j];
		}
	}
	fout.Close();
	std::cout << "... done" << std::endl;
	return;
}

void NuFit::toymc::run_simulation_waux(std::string outfile, std::string outfile_waux, std::map<std::string, double> simpars, unsigned int nsamples) {
    std::cout << "\n";
    std::cout << "... sampling " << nsamples << " toymc datasets for fixed set of model parameters." << std::endl;
    TFile fout(outfile.c_str(), "RECREATE");
	
	std::ofstream file;
	file.open(outfile_waux);
	file.precision(5);
	file.setf(std::ios_base::scientific);

	file << "dom_efficiency" << "\t" << "delta_cr" << "\t" << "absorption" << "\t" << "scattering" << "\t" << "holeicep0" << "\t" << "holeicep1" << "\t" << "selfveto" << "\t" << "hadronicinteraction" << "\n";	
	//file << "dom_efficiency" << "\t" << "delta_cr" << "\t" << "absorption" << "\t" << "scattering" << "\t" <<  "muon_norm_mlb" << "\n";
    //file << "dom_efficiency" << "\t" << "delta_cr" << "\n";

	double means_simulated[2] = {simpars["absorption"], simpars["scattering"]};
	double domeff_simulated = simpars["dom_efficiency"];
	double deltacr_simulated = simpars["delta_cr"];
	double holeicep0_simulated = simpars["holeicep0"];
	double holeicep1_simulated = simpars["holeicep1"];
    double selfveto_simulated = simpars["selfveto"];
    double hadronicinteraction_simulated = simpars["hadronicinteraction"];

	double prior_domeff_sigma = 0.1;
	double prior_deltacr_sigma = 0.05; 
	double prior_holeicep0_sigma = 1;
	double prior_holeicep1_sigma = 0.2;
    double prior_selfveto_sigma = 1000;
    double prior_hadronicinteraction_sigma = 2;
    std::cout<<simpars.size();
    for (unsigned int i=0; i<nsamples; ++i)
    {
        std::vector<TH3D*> hists = draw_sample(simpars);
        for (unsigned int j=0; j<hists.size(); ++j)
        {
            // ROOT uses C-style character arrays
            TH3D *hist = hists[j];
            std::string name = std::string(hist->GetName(), 0, 1000);
            hist->SetName((std::string(name)+std::string("_")+std::to_string(i)).c_str());
            hist->Write();
            delete hists[j];
        }

	// update auxiliary data
	double prior_domeff_mean_current = generate_truncated_normal_rv(domeff_simulated, prior_domeff_sigma);
	double prior_deltacr_mean_current = generate_normal_rv(deltacr_simulated, prior_deltacr_sigma);
	double prior_holeicep0_mean_current = generate_truncated_normal_rv(holeicep0_simulated, prior_holeicep0_sigma);
	double prior_holeicep1_mean_current = generate_truncated_normal_rv(holeicep1_simulated, prior_holeicep1_sigma);
    double prior_selfveto_mean_current = generate_truncated_normal_rv(selfveto_simulated, prior_selfveto_sigma);
    double prior_hadronicinteraction_mean_current = generate_truncated_normal_rv(hadronicinteraction_simulated, prior_hadronicinteraction_sigma);

	double rvars[2] = {0.0, 0.0};
	generate_truncated_bivnorm_rv(means_simulated, rvars);
	double prior_absorption_current = rvars[0];
	double prior_scattering_current = rvars[1];

	file << prior_domeff_mean_current << "\t" << prior_deltacr_mean_current << "\t" << prior_absorption_current << "\t" << prior_scattering_current << "\t" << prior_holeicep0_mean_current << "\t" << prior_holeicep1_mean_current << "\t" << prior_selfveto_mean_current << "\t" << prior_hadronicinteraction_mean_current << "\n";
    }

    fout.Close();
	file.close();
    std::cout << "... done" << std::endl;
    return;	
}

void NuFit::toymc::run_simulation_bayesian(std::string outfile, std::map<std::string, double> simpars_fixed, std::vector<std::string> simpars_float, std::string infile_parameter_samples, unsigned int nsamples)
{
	std::cout << "\n";
	std::cout << "... reading parameter samples (corresponding to prior or posterior distribution)" << std::endl;
	// check if we have prameters for all models
	if( (simpars_fixed.size() + simpars_float.size()) != analysis.get_npars() )
	{
		// we have unspecified model parameters. exit.
		std::cout << "FATAL !!! user input specified " << simpars_fixed.size() + simpars_float.size() << " parameters. but model expects " << analysis.get_npars() << " parameters." << std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	// create default point container
	std::vector<std::map<std::string, double>> points;
	std::map<std::string, double> point;
	
	// add the fix points
	for (std::map<std::string, double>::const_iterator it = simpars_fixed.begin(); it != simpars_fixed.end(); ++it)
		point.insert(std::pair<std::string, double>(it->first, it->second));

	// add the floating points
	for (unsigned int i=0; i<simpars_float.size(); ++i)
		point.insert(std::pair<std::string, double>(simpars_float[i], 0.0));


	// check if infile exists
	std::string line;
	std::fstream thisfile(infile_parameter_samples);
	unsigned int ncols = 0;
	if (thisfile.is_open())
	{
		getline(thisfile, line); // assume first line has column names

		std::istringstream ss(line);
		std::string col_name;

		std::map<std::string, int> column_index;

		while (ss >> col_name) 
		{			
			column_index.insert(std::pair<std::string, int>(col_name, ncols));
			ncols++;
		}

		// make sure the columns exist in point
		for(unsigned int i=0; i<simpars_float.size(); ++i)	
		{
			if(column_index.find(simpars_float[i]) == column_index.end())
			{
				// this is logical error. file header does not contain requested variable
				std::cout << "FATAL !!! The file does not contain a column with variable " << simpars_float[i] << std::endl;
				std::cout << "... exiting" << std::endl;
				exit(1);	
			}
		}
	
		std::vector<double> current_parameters;
		current_parameters.resize(ncols);

		unsigned int counter = 0;
		while(getline(thisfile, line))
		{
			std::stringstream ss(line);
			unsigned int col=0;
			while (ss >> current_parameters[col]) col++;

			if( col!=ncols ) {
                                std::cout << "FATAL! unexpected number of columns!" << std::endl;
                                std::cout << "... exiting" << std::endl;
                                exit(1);
                        }	

			// update point
			for (unsigned int i=0; i<simpars_float.size(); ++i)
				point[simpars_float[i]] = current_parameters[column_index[simpars_float[i]]];
			
			points.push_back(point);

			counter++;
			if(counter >= nsamples)
				break;
		}
		
	}
	else
        {
                std::cout << "FATAL! unable to open file: " << infile_parameter_samples << std::endl;
                std::cout << "... exiting" << std::endl;
                exit(1);
        }

	// lets crosscheck wether we parsed correctly.
	//for (unsigned int i=0; i<points.size(); ++i) 
	//{
	//	std::map<std::string, double> &tpoint = points[i];
	//	for(std::map<std::string, double>::const_iterator it=tpoint.begin(); it!=tpoint.end(); ++it)
	//		std::cout << it->first << " " << it->second << ", ";
	//	std::cout << "\n";
	//}
	//std::cout << std::endl;
	
	// run toy generation
	std::cout << "\n";
	std::cout << "... sampling " << nsamples << " toymc datasets for fixed set of model parameters." << std::endl;
	TFile fout(outfile.c_str(), "RECREATE");
	for (unsigned int i=0; i<points.size(); ++i)
	{
		std::vector<TH3D*> hists = draw_sample(points[i]);
			for (unsigned int j=0; j<hists.size(); ++j)
		{
			// ROOT uses C-style character arrays
			TH3D *hist = hists[j];
			std::string name = std::string(hist->GetName(), 0, 1000);
			hist->SetName((std::string(name)+std::string("_")+std::to_string(i)).c_str());
			hist->Write();
			delete hists[j];
		}	
	}
	fout.Close();
	std::cout << "... done" << std::endl;
	return;
	
}

double NuFit::toymc::generate_truncated_normal_rv(const double &mean, const double &sigma) 
{
	// 0-truncated normal random variable
	double rv = -1.0;
	while (rv<=0.0)
		rv = generate_normal_rv(mean, sigma);
	
	return rv;
}

void NuFit::toymc::generate_truncated_bivnorm_rv(double means[2], double rvars[2])
{
	// 0-truncated bivariate normal random variable
        
	// Spice paper ellipse
	// double rho=-0.1;
	double sigma=0.07;
	
	double asin_rho = -0.1001674211615598; // asin(rho)

	rvars[0]=-1.0;
	rvars[1]=-1.0;

	while ((rvars[0]<=0.0)||(rvars[1]<=0.0)) {
	// Box Muller Transforms to get (0,0) centered RV
		std::uniform_real_distribution<double> unif(0.0, 1.0);
		double u1 = unif(rand);
		double u2 = unif(rand);
		
		double v1 = sigma * TMath::Sqrt(-2.0 * TMath::Log(u2)) * TMath::Cos(TMath::Pi() * 2.0 * u1);
		double v2 = sigma * TMath::Sqrt(-2.0 * TMath::Log(u2)) * TMath::Sin(TMath::Pi() * 2.0 * u1 + asin_rho);
		rvars[0] = v1 + means[0]; // move to mean
		rvars[1] = v2 + means[1];
	}
	
	return;
}

double NuFit::toymc::generate_normal_rv(const double &mean, const double &sigma) 
{
	// normal random variable
	std::normal_distribution<double> norm_rv(mean, sigma);
	return norm_rv(rand);
}
