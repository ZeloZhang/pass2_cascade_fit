#include "../include/stats.h"

NuFit::stats::stats(NuFit::analysis &analysis_):
	target_func(&analysis_, &NuFit::analysis::get_likelihood, analysis_.get_npars()),
	analysis( analysis_ )
{
	min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); // supports all minimizers available via ROOT (including ROOT's interface to GSL minimizers). 
	//min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "SteepestDescent"); // supports all minimizers available via ROOT (including ROOT's interface to GSL minimizers). 

	min->SetFunction(target_func); // likelihood function without terms that have no parameter dependence	
	min->SetMaxFunctionCalls(1000000); // ROOT default
	min->SetMaxIterations(100000); // ROOT default
	min->SetTolerance(0.001); // ROOT default, can be changed using NuFit::stats::set_tolerance()

	flush_rate = 20;

	return;
} 

NuFit::stats::~stats() 
{
	delete min;
}

void NuFit::stats::set_options(std::map<std::string, NuFit::helpers::par_options> opts) 
{
	seeds.clear();
	stepsizes.clear();
	limits_low.clear();
	limits_high.clear();

	unsigned int npars = analysis.get_npars();

	// check whether we have received options for all parameter
	if(opts.size()!=npars) {
		std::cout << "!!! FATAL: received unexpected number of parameters: " << opts.size() << " instead of: " << npars << std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	// get model parameters and set user options
	std::vector<std::string> names;
	analysis.get_par_names(names);

	for(unsigned int i=0; i<npars; ++i) {
		// check wether parameter name exists in user options
		if (opts.count(names[i])==0) {
			std::cout << "!!! FATAL: parameter " << names[i] << " not found. Did you specify it correctly?" << std::endl;
			std::cout << "!!! only found the following parameters:" << std::endl;
			for (std::map<std::string, NuFit::helpers::par_options>::iterator it = opts.begin(); it != opts.end(); ++it )
				std::cout << it->first << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}

		// name exists
		seeds.push_back(opts.at(names[i]).seed);
		stepsizes.push_back(opts.at(names[i]).stepsize);
		limits_low.push_back(opts.at(names[i]).limit_low);
		limits_high.push_back(opts.at(names[i]).limit_high);	

	}
	return;
}

void NuFit::stats::set_tolerance(const double tol)
{
	min->SetTolerance(tol);
	return;
}

void NuFit::stats::fit(bool werrors)
{
	unsigned int npars = analysis.get_npars();
	min->Clear();
	bestpars.clear();

	std::vector<std::string> par_names;
	analysis.get_par_names(par_names);

	for (unsigned int i=0; i<npars; ++i)
		min->SetLimitedVariable(i, par_names[i].c_str(), seeds[i], stepsizes[i], limits_low[i], limits_high[i]);


	std::cout << std::endl;
	std::cout << "... commencing numerical minimization (" << npars << " parameters)" << std::endl;
	std::cout << std::endl;

	// time for full likelihood evaluation
	std::clock_t start;
	double duration;
	
	start = std::clock();
	min->Minimize();

	const double *xs = min->X();

	std::cout << "... done" << std::endl;
	std::cout << std::endl;
	std::cout << "RESULTS (fitstatus: " << min->Status() << ", edm: "<< min->Edm() << ")" << std::endl;

	// set seed to best fit (for subsequent minimizations)
	
	seeds.clear();
	for (unsigned int i=0; i<npars; ++i) {
                std::cout << par_names[i] << " = " << xs[i] << std::endl;
		bestpars.push_back(xs[i]);
		seeds.push_back(xs[i]);	
	}


	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	std::cout << "likelihood value (gof): " << analysis.get_likelihood_gof(xs) << std::endl;
	std::cout << "likelihood value (abs): " << analysis.get_likelihood_abs(xs) << std::endl;
	std::cout << "time to find minimium and evaluate exact likelihood: " << duration << "s" << std::endl;
	std::cout << "number of llh evaluations: " << analysis.get_n_llh_evals() << std::endl;
	//std::cout << "saturated gof statistic: " << analysis.get_saturated_gof(xs) << std::endl;
	std::cout << std::endl;

	
	if (werrors) {
		std::cout << "... calculating profile LLH approximate confidence intervals via ROOT interface" << std::endl;
		double errs_low[npars];
		double errs_up[npars];

		for(unsigned int i=0; i<npars; ++i) {
			errs_low[i]=0;
			errs_up[i]=0;
			min->GetMinosError(i, errs_low[i], errs_up[i], 0);
			std::cout << "... " << par_names[i] << " done" << std::endl;
		}

		std::cout << "... done" << std::endl;
		std::cout << std::endl;

		std::cout << "RESULTS" << std::endl;	
		for(unsigned int i=0; i<npars; ++i) 
		 	std::cout << par_names[i] << " = " << bestpars[i] << " (" << errs_low[i] << " +" << errs_up[i] << ")" << std::endl;

	}	
	
	return;
}

void NuFit::stats::get_bestpars(std::map<std::string, double> &pars)
{
	// check if number of parameters requested matches model parameters
	unsigned int npars = analysis.get_npars();
	if (pars.size() != npars) {
		std::cout << "!!! received npars=" << pars.size() << " but model has: " << npars <<" parameters" << std::endl;
		std::cout << "exiting ..." << std::endl;
		exit(1);
	}

	std::map<std::string, unsigned int> parameters;
	analysis.get_par_names(parameters);

	for (std::map<std::string, unsigned int>::iterator it=parameters.begin(); it != parameters.end(); ++it) 
	{
		// check if parameter exists
		if(pars.count(it->first)==0) 
		{
			std::cout << "!!! can not find model parameter: " << it -> first << std::endl;
			std::cout << "!!! the following parameters are found: " << std::endl;
			for (std::map<std::string, double>::iterator it2=pars.begin(); it2 != pars.end(); ++it2)
				std::cout << it2->first << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}
		pars.at(it->first)=bestpars.at(it->second);
	}
	
	return;
}

void NuFit::stats::change_astro_model(NuFit::astro_model_base *astro) 
{
	analysis.change_astro_model(astro);
	// need to update ROOT functor to update total number of parameter dimensions
	target_func = ROOT::Math::Functor(&analysis, &NuFit::analysis::get_likelihood, analysis.get_npars());
	min->SetFunction(target_func);
}

double NuFit::stats::get_profile_llh_value(std::map<std::string, double> &point) {
	return  NuFit::stats::get_profile_llh_value(point, true, false);
}

double NuFit::stats::get_profile_llh_value(std::map<std::string, double> &point, bool verbose, bool scan) 
{
	if(verbose) 
	{
		unsigned int npars = analysis.get_npars();
		// this is a single profile_llh evaluation. not a scan!
		// be more verbose and do more error checking on user input
		// also since this is not a scan, we need to setup the minimizer here
		
		std::cout << "... requested profile llh" << std::endl;	
		min->Clear();
		if( npars<=point.size() ) 
		{
			// user input can not be profiled
			std::cout << "!!! FATAL: can not calculate profile likelihood for dim(point) >= number of model parameters. For dim(point) == number of model parameters use stats::get_llh_value(std::map<std::string, double> &point) instead!" << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}

		std::vector<std::string> par_names;
		analysis.get_par_names(par_names);
		std::map<std::string, unsigned int> parameters;
		analysis.get_par_names(parameters);

		// figure out what parameters need to be fixed to a point
		std::vector<unsigned int> indices;
		for (std::map<std::string, double>::iterator it=point.begin(); it!=point.end(); ++it)
		{
			if(parameters.count(it->first)==0) 
			{
				// point specified by user concerns a parameter that does not exist :-(
				std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
				std::cout << "!!! the model has the following parameters:" << std::endl;
				for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
					std::cout << it2->first << std::endl;
				std::cout << "... exiting" << std::endl;
				exit(1);
			}

			// user provided a valid point
			indices.push_back(parameters.at(it->first));
		}

		for (unsigned int i=0; i<npars; ++i) 
		{
			// fix requested variables	
			if(std::find(indices.begin(), indices.end(), i) != indices.end())
			{
				std::cout << "... treating variable: " << par_names[i] << " with index " << i << " as constant." << std::endl;
				min->SetFixedVariable(i, par_names[i].c_str(), point.at(par_names[i]));
			}
			else
				min->SetLimitedVariable(i, par_names[i].c_str(), seeds[i], stepsizes[i], limits_low[i], limits_high[i]);
				
		}	

		std::cout << std::endl;
		std::cout << "... commencing numerical minimization (" << npars-point.size() << " parameters)" << std::endl;
		std::cout << std::endl;
	}

	// this be executed for both cases:
	// single point evaluation by user as well as single point evaluation as part of a llh scan

	min->Minimize();

	const double *xs = min->X();
	double profile_llh = analysis.get_likelihood_gof(xs);
	
	unsigned int npars = analysis.get_npars();
	std::vector<std::string> par_names;
	analysis.get_par_names(par_names);

	if(verbose) 
	{
		// this is not a scan
		// print some information for user 
		std::cout << "... done" << std::endl;
		std::cout << std::endl;
		std::cout << "RESULTS (fitstatus: " << min->Status() << ", edm: "<< min->Edm() << ")" << std::endl;

		for (unsigned int i=0; i<npars; ++i) 
			std::cout << par_names[i] << " = " << xs[i] << std::endl;

	}

	else
	{
		// this is a scan. or toyfit.  need to update containers:
		for (unsigned int i=0; i<npars; ++i)	
			point.at(par_names[i])=xs[i]; 	
	
		point.at("fit_status")=min->Status();
		point.at("llh")=profile_llh;	

		// set next seed point
		if(scan)
		{
		// the next scan point will be in neighborhood of this solution.
                // set seed to current best fit for next iteration, also update point with current solution
			if(min->Status()==0) 
			{
				for (unsigned int i=0; i<npars; ++i)
					seeds[i]=xs[i];
			}
			else
			{	
				std::cout << " ... fit failed. seed next scan point from global best fit." << std::endl;
				for (unsigned int i=0; i<npars; ++i)
					seeds[i]=bestpars[i];
			}
		}

	}

	return profile_llh;
}

void NuFit::stats::scan_llh(std::string outfile, std::map<std::string, NuFit::helpers::scan_options> &opts, std::map<std::string, double> &seeds) {
	std::cout << std::endl;
	std::cout << "... requested simple llh scan." << std::endl;

	// check if file exists
	if (std::ifstream(outfile))
	{
		std::cout << "!!! FATAL: File " << outfile << " already exists. Refuse to overwrite!" << std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	// check if can write to file
	std::ofstream file(outfile);
	file.precision(5);
	file.setf(std::ios_base::scientific);
	if (!file)
	{
		std::cout << "!!! FATAL: unable to write to file " << outfile << std::endl;
		exit(1);	
	}	

	unsigned int npars = analysis.get_npars();
	if( npars != (opts.size()+seeds.size()) ) 
	{
		std::cout << "!!! FATAL: can not calculate profile likelihood for number of model parameters != number of scan parameters + number of fixed parameters." << std::endl;
		std::cout << "number of model parameters: " << npars << " number of scan parameters: " << opts.size() << " number of fixed parameters: " << seeds.size() << std::endl; 
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	std::vector<std::string> par_names;
	analysis.get_par_names(par_names);
	std::map<std::string, unsigned int> parameters;
	analysis.get_par_names(parameters);

	for(unsigned int i=0; i<par_names.size(); ++i)
		file << par_names[i] << "\t";
	file << "-2lnL\tfit_status" << std::endl;

	std::vector<std::vector<double>> grid_axes;
	std::vector<std::string> axis_names;
	unsigned int ncalls=1;

	for (std::map<std::string, NuFit::helpers::scan_options>::iterator it=opts.begin(); it!=opts.end(); ++it)
	{
		if(parameters.count(it->first)==0) 
		{
			// parameters specified by user is not a model parameter:-(
			std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
			std::cout << "!!! the model has the following parameters:" << std::endl;
			for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
				std::cout << it2->first << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}	

                unsigned int nsteps = (it->second).nsteps;
                double range_low = (it->second).range_low;
                double range_high = (it->second).range_high;
                double stepsize = (range_high - range_low) / (nsteps-1);

		// create axis
		std::vector<double> axis;
		for(unsigned int i=0; i<nsteps; ++i)
			axis.push_back(range_low + i * stepsize);

		grid_axes.push_back(axis);
		axis_names.push_back(it->first);
		ncalls*=axis.size();
	}

	// inform user about number of llh calculations implied by his/her grid
	std::cout << "... scaning in " << grid_axes.size() << " dimensions. specified grid with " << ncalls << " points." << std::endl;
	std::cout << "(number of grid points = number of llh evaluations)" << std::endl;
	std::cout << std::endl;

	// create hypothesis
	double these_pars[npars];

        // create scan point - this will be updated to solution after each LLH evaluation
        std::map<std::string, double> point;
        for(unsigned int i=0; i<par_names.size(); ++i)
                point.insert(std::pair<std::string, double>(par_names[i], 0.0));

	// set first the parameter values that are to be kept fixed
	for (std::map<std::string, double>::const_iterator it=seeds.begin(); it!=seeds.end(); ++it) 
	{
		if(parameters.count(it->first)==0)
		{
			std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
			std::cout << "!!! the model has the following parameters:" << std::endl;
			for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
				std::cout << it2->first << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}

		unsigned int par_index = parameters[it->first];
		these_pars[par_index] = it->second;
		point[it->first] = it->second;
	}

	// and here goes the complicated part. need to loop over a d-dimensional dynamic grid (size unknown at compile time);
	// algorithm adapted from https://www.codeproject.com/Tips/1004156/Iterating-Over-Dynamic-Number-of-Nested-Loops

	unsigned int size = grid_axes.size();
	int sizet = grid_axes.size();
	unsigned int counters[size];
	unsigned int bounds[size];

	int loop_index = 0;
	int n_evals = 0;

	for(unsigned int i=0; i<size; ++i)
	{
		counters[i]=0;
		bounds[i]=grid_axes[i].size();
	}	

	// also store likelihood
	point.insert(std::pair<std::string, double>(std::string("llh"), 0.0));
	// store dummmy fit status (this just evaluating, so no fitting is actually being done
	point.insert(std::pair<std::string, double>(std::string("fit_status"), -1.0));

	// write to file every 50th grid point. need to cache intermediate results.
	std::vector<std::map<std::string, double>> points;
	points.reserve(flush_rate); 

	do
	{	
		++n_evals;
		for(unsigned int i=0; i<size; ++i) 
		{	
			// update current parameter array	
			these_pars[parameters.at(axis_names[i])] = grid_axes[i][counters[i]];
			point[axis_names[i]] = grid_axes[i][counters[i]];
		}

		
		double llh = analysis.get_likelihood_gof(these_pars);
		point.at("llh")=llh;

		points.push_back(point);
		if(n_evals%flush_rate==0) {
			// flush result to file
			write_profile_llh(file, points);
			std::cout << "... " << n_evals << "/" << ncalls << " grid points done." << std::endl;
			points.clear();
		}
			
	} while (increment_loop_state(counters, bounds, loop_index, sizet));

	// need to flush remaining cached results to file
	write_profile_llh(file, points);
	file.close();
	std::cout << "... finished llh scan." << std::endl;
	return;
	
}


void NuFit::stats::scan_profile_llh(std::string outfile, std::map<std::string, NuFit::helpers::scan_options> &opts) 
{
	std::cout << std::endl;
	std::cout << "... requested profile llh scan." << std::endl;
	// check if file exists
	if (std::ifstream(outfile))
	{
		std::cout << "!!! FATAL: File " << outfile << " already exists. Refuse to overwrite!" << std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	// check if can write to file
	std::ofstream file(outfile);
	file.precision(5);
	file.setf(std::ios_base::scientific);
	if (!file)
	{
		std::cout << "!!! FATAL: unable to write to file " << outfile << std::endl;
		exit(1);	
	}

	unsigned int npars = analysis.get_npars();

	if( npars<=opts.size() ) 
	// can profile only if number of scan paramters < model parameters
	{
		std::cout << "!!! FATAL: can not calculate profile likelihood for number of scan parameters >= number of model parameters." << std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	std::vector<std::string> par_names;
	analysis.get_par_names(par_names);
	std::map<std::string, unsigned int> parameters;
	analysis.get_par_names(parameters);	

	// write header information to file
	for(unsigned int i=0; i<par_names.size(); ++i)
		file << par_names[i] << "\t";
	file << "-2lnL\tfit_status" << std::endl;

	// figure out what parameters need to be scanned	
	std::vector<unsigned int> indices;

	// a d-dimensional profile-llh requires a grid. construct grid axes from user input
	std::vector<std::vector<double>> grid_axes;
	std::vector<std::string> axis_names;
	std::map<std::string, double> grid_point;
	unsigned int ncalls=1;

	for (std::map<std::string, NuFit::helpers::scan_options>::iterator it=opts.begin(); it!=opts.end(); ++it)
	{
		if(parameters.count(it->first)==0) 
		{
			// parameters specified by user is not a model parameter:-(
			std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
			std::cout << "!!! the model has the following parameters:" << std::endl;
			for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
				std::cout << it2->first << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}
		// user provided a valid parameter
		indices.push_back(parameters.at(it->first));	

        unsigned int nsteps = (it->second).nsteps;
        double range_low = (it->second).range_low;
        double range_high = (it->second).range_high;

		double stepsize = 0;
		if (nsteps>1){
            stepsize = (range_high - range_low) / (nsteps-1);
        }

		// create axis
		std::vector<double> axis;
		for(unsigned int i=0; i<nsteps; ++i)
			axis.push_back(range_low + i * stepsize);

		grid_axes.push_back(axis);
		axis_names.push_back(it->first);
		grid_point.insert(std::pair<std::string, double>(it->first, 0.0)); // actual value will be overwritten in loop
		ncalls*=axis.size();
	}

	// inform user about number of profile llh calculations implied by his/her grid
	std::cout << "... scaning in " << grid_axes.size() << " dimensions. specified grid with " << ncalls << " points." << std::endl;
	std::cout << "(number of grid points = number of profile-llh evaluations)" << std::endl;
	std::cout << std::endl;

	// set up minimizer
	min->Clear();
        for (unsigned int i=0; i<npars; ++i)
        {
                // fix requested variables
                if(std::find(indices.begin(), indices.end(), i) != indices.end())
                {
                        std::cout << "... treating variable: " << par_names[i] << " with index " << i << " as constant." << std::endl;
                        min->SetFixedVariable(i, par_names[i].c_str(), seeds[i]); // temporary seed
                }
                else
                        min->SetLimitedVariable(i, par_names[i].c_str(), seeds[i], stepsizes[i], limits_low[i], limits_high[i]);
        }

	// and here goes the complicated part. need to loop over a d-dimensional dynamic grid (size unknown at compile time);
	// algorithm adapted from https://www.codeproject.com/Tips/1004156/Iterating-Over-Dynamic-Number-of-Nested-Loops

	unsigned int size = grid_axes.size();
	int sizet = grid_axes.size();
	unsigned int counters[size];
	unsigned int bounds[size];

	int loop_index = 0;
	int n_evals = 0;

	for(unsigned int i=0; i<size; ++i)
	{
		counters[i]=0;
		bounds[i]=grid_axes[i].size();
	}

	// create scan point - this will be updated to solution after each minimization
	std::map<std::string, double> point;
	for(unsigned int i=0; i<par_names.size(); ++i)
		point.insert(std::pair<std::string, double>(par_names[i], 0.0));

	// also store likelihood
	point.insert(std::pair<std::string, double>(std::string("llh"), 0.0));
	point.insert(std::pair<std::string, double>(std::string("fit_status"), -1.0));

	// write to file every 50th grid point. need to cache intermediate results.
	std::vector<std::map<std::string, double>> points;
	points.reserve(flush_rate); // to hold 50 elements

	do
	{	
		++n_evals;
		update_minimizer(grid_point);
		for(unsigned int i=0; i<size; ++i) 
		{
			// update minimizer state, need to find variable index
			min->SetVariableValue(parameters.at(axis_names[i]), grid_axes[i][counters[i]]);

		}
		get_profile_llh_value(point, false, true); // since this a scan profile_llh is already stored in point
		
		// some prints
		//for (std::map<std::string, double>::iterator it=point.begin(); it!=point.end(); ++it)
		//	std::cout << it->first << " " << it->second << " ";	
		//std::cout << std::endl;

		// cache result
		points.push_back(point);
		if(n_evals%flush_rate==0) {
			// flush result to file
			write_profile_llh(file, points);
			std::cout << "... " << n_evals << "/" << ncalls << " grid points done." << std::endl;
			points.clear();
		}

	} while (increment_loop_state(counters, bounds, loop_index, sizet));

	// need to flush remaining cached results to file
	write_profile_llh(file, points);
	file.close();
	std::cout << "... finished profile llh scan." << std::endl;
	return;

}

void NuFit::stats::scan_profile_llh_asimov(std::string outfile_scan, std::map<std::string, NuFit::helpers::scan_options> &opts, std::string infile_asimov_hists)
{

    // save original data hists
    analysis.cache_data_hists();

    // let's read histograms back from file and exchange for original data histograms
    TFile rootfile(infile_asimov_hists.c_str(), "READ");

    // need to read mcsum for all event selections and inject to analyses
    std::vector<std::string> names = analysis.get_analysis_names();
    for (unsigned int j=0; j<names.size(); ++j)
    {
        // read from file
        std::string histname = names[j] + std::string("_mcsum");
        TH3D* datahist = (TH3D*) rootfile.Get(histname.c_str());
        datahist->SetName((names[j]+std::string("_data")).c_str());
        analysis.set_hist(names[j], datahist);
        delete datahist;
    }
    rootfile.Close();

    // fit
    fit(false);
    // run profile scan
    scan_profile_llh(outfile_scan, opts);

    // copy original data histograms back (remove toy data from analysis)
    analysis.restore_data_hists();
    std::cout << std::endl;
    return;
}
void NuFit::stats::scan_profile_llh_asimov(std::string outfile_scan, std::map<std::string, NuFit::helpers::scan_options> &opts, std::string outfile_asimov_hists, std::map<std::string, double> &asimov_point)
{
	std::cout << std::endl;
	std::cout << "... requested profile llh scan on the asimov dataset." << std::endl;

	// save original data hists 
	analysis.cache_data_hists();

	// get asimov dataset for user provided asimov_point to file
	analysis.get_histograms(outfile_asimov_hists, asimov_point);

	// let's read histograms back from file and exchange for original data histograms
	TFile rootfile(outfile_asimov_hists.c_str(), "READ");

	// need to read mcsum for all event selections and inject to analyses
	std::vector<std::string> names = analysis.get_analysis_names();
	for (unsigned int j=0; j<names.size(); ++j)
	{
		// read from file
		std::string histname = names[j] + std::string("_mcsum");
		TH3D* datahist = (TH3D*) rootfile.Get(histname.c_str());
		datahist->SetName((names[j]+std::string("_data")).c_str());
		analysis.set_hist(names[j], datahist);
		delete datahist;
	}
	rootfile.Close();

	// overwrite rootfile with new datahist (after injection)
	analysis.get_histograms(outfile_asimov_hists, asimov_point);

	// run profile scan
	scan_profile_llh(outfile_scan, opts);
	
	// copy original data histograms back (remove toy data from analysis)
	analysis.restore_data_hists();
	std::cout << std::endl;
	return;
}

bool NuFit::stats::increment_loop_state(unsigned int *counters, const unsigned int *bounds, int &loop_index, int &size) {
	if (loop_index >= size)
		return false;
	
	counters[loop_index]+=1;

	bool result = true;
	
	if (counters[loop_index] >= bounds[loop_index])
	{
		counters[loop_index]=0;
		++loop_index;
		result = increment_loop_state(counters, bounds, loop_index, size);
		--loop_index;
	}

	return result;
	
}

void NuFit::stats::write_profile_llh(std::ofstream &file, const std::vector<std::map<std::string, double>> &points) 
{
	std::vector<std::string> par_names;
	analysis.get_par_names(par_names);

	std::string llh("llh");
	std::string fit_status("fit_status");
	
	for(unsigned int i=0; i<points.size(); ++i) 
	{
		for(unsigned int j=0; j<par_names.size(); ++j) 
		{
			file << points.at(i).at(par_names[j]) << "\t";
		}

		file << points.at(i).at(llh) << "\t";
		file << (int)points.at(i).at(fit_status) << std::endl;
	}

	return;

}

void NuFit::stats::write_toy_fits(std::ofstream &file, const std::vector<std::map<std::string, double>> &points)
{
        std::vector<std::string> par_names;
        analysis.get_par_names(par_names);

        std::string llh("llh");
        std::string fit_status("fit_status");
	std::string toy_index("toy_index");

        for(unsigned int i=0; i<points.size(); ++i)
        {
		file << (int)points.at(i).at(toy_index) << "\t";

                for(unsigned int j=0; j<par_names.size(); ++j)
                {
                        file << points.at(i).at(par_names[j]) << "\t";
                }

                file << points.at(i).at(llh) << "\t";
                file << (int)points.at(i).at(fit_status) << std::endl;
        }

        return;

}

void NuFit::stats::set_flush_rate(unsigned int rate) 
{
	flush_rate = rate;
	return;
}

double NuFit::stats::get_llh_value_abs(const std::map<std::string, double> &point)
{
        unsigned int npars = analysis.get_npars();
        if( npars!=point.size() )
        {
                std::cout << "!!! FATAL: can not calculate profile likelihood for dim(point) != number pf model parameters" << std::endl;
                std::cout << "... exiting" << std::endl;
                exit(1);
        }

        // need to order parameters to pass to likelihood
        std::map<std::string, unsigned int> parameters;
        analysis.get_par_names(parameters);

        double pars[npars];

        for (std::map<std::string, double>::const_iterator it=point.begin(); it!=point.end(); ++it)
        {
                if(parameters.count(it->first)==0)
                {
                        std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
                        std::cout << "!!! the model has the following parameters:" << std::endl;
                        for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
                                std::cout << it2->first << std::endl;
                        std::cout << "... exiting" << std::endl;
                        exit(1);
                }
                // parameter exists
                unsigned int par_index = parameters[it->first];
                pars[par_index]=it->second;
        }

        return analysis.get_likelihood_abs(pars);
}

double NuFit::stats::get_llh_value(const std::map<std::string, double> &point)
{
	unsigned int npars = analysis.get_npars();
	if( npars!=point.size() ) 
	{
		std::cout << "!!! FATAL: can not calculate profile likelihood for dim(point) != number pf model parameters" << std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	// need to order parameters to pass to likelihood
	std::map<std::string, unsigned int> parameters;
	analysis.get_par_names(parameters);

	double pars[npars];

	for (std::map<std::string, double>::const_iterator it=point.begin(); it!=point.end(); ++it)
	{
		if(parameters.count(it->first)==0)
		{
			std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
			std::cout << "!!! the model has the following parameters:" << std::endl;
			for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
				std::cout << it2->first << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}
		// parameter exists
		unsigned int par_index = parameters[it->first];
		pars[par_index]=it->second;
	}

	return analysis.get_likelihood_gof(pars);
}

void NuFit::stats::run_toyfits(std::string infile, std::string outfile, std::map<std::string, double> pars, std::map<std::string, double> seed_vals, unsigned int toy_id_low, unsigned int toy_id_high, bool randomize_auxilliary_data, std::string infile_aux)
{
	std::cout << std::endl;
	std::cout << "... requested fits on simulated toy-data." << std::endl;

	if(toy_id_high <= toy_id_low)
	{
		std::cout << "!!! FATAL:  toy_id_high <= toy_id_low. can't perform fits." << std::endl;
		std::cout << "... exiting" << std::endl;
	}

	unsigned int nfits = toy_id_high - toy_id_low;

	// check if file exists
	if (std::ifstream(outfile))
	{
		std::cout << "!!! FATAL: File " << outfile << " already exists. Refuse to overwrite!" << std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}	

	// check if can write to file
	std::ofstream file(outfile);
	file.precision(5);
	file.setf(std::ios_base::scientific);
	if (!file)
	{
		std::cout << "!!! FATAL: unable to write to file " << outfile << std::endl;
		exit(1);	
	}	

	// how many parameters?
	unsigned int npars_requested = pars.size();
	unsigned int npars_model = analysis.get_npars();

	if(npars_requested > npars_model)
	{
		std::cout << "!!! FATAL: received more parameters (" << npars_requested << ") than specified by model (" << npars_model << ")" <<std::endl;
		std::cout << "... exiting" << std::endl;
		exit(1);
	}

	// check if all specified parameters are actually in the model

    std::cout << "... start fitting " << nfits << " artificial toy datasets." << std::endl;

    // need analysis names to obtain analysis names
    std::vector<std::string> names = analysis.get_analysis_names();

    // save original data hists 
	analysis.cache_data_hists();

	// need to set up minimizer
	update_minimizer(pars);

	TFile rootfile(infile.c_str(), "READ");

	std::vector<std::string> par_names;
	analysis.get_par_names(par_names);
	std::map<std::string, double> point;
	for(unsigned int i=0; i<par_names.size(); ++i){
		point.insert(std::pair<std::string, double>(par_names[i], 0.0));
    }
	 
	point.insert(std::pair<std::string, double>("llh", 0.0));
	point.insert(std::pair<std::string, double>("fit_status", -1.0));
	point.insert(std::pair<std::string, double>("toy_index", 0.0));

	// write to file every 50th grid point. need to cache intermediate results.
	std::vector<std::map<std::string, double>> points;
	points.reserve(flush_rate);

    // write header
    file << "toy_id\t";
    for(unsigned int i=0; i<par_names.size(); ++i){
        file << par_names[i] << "\t";
    }
    file << "-2lnL\tfit_status" << std::endl;

	// if randomize auxilliary data, need to read randomized experiments
	std::vector< std::map<std::string, double> > aux_data;

	if(randomize_auxilliary_data) {
		std::string line_aux;	
		std::fstream thisfile_aux(infile_aux);
		if (thisfile_aux.is_open()) {
			// first line contains variable headers
			getline(thisfile_aux, line_aux);
			std::stringstream ss(line_aux);
			std::string var;
			std::vector<std::string> var_names;
			while ( ss >> var ) {
				var_names.push_back(var);	
			}

			std::map<std::string, double> aux_pars;
			for(unsigned int i=0; i<var_names.size(); ++i)
				aux_pars.insert(std::pair<std::string, double>(var_names[i],0.0));

			// now read data
			while(getline(thisfile_aux, line_aux)) {
				std::stringstream ss1(line_aux);
				double value[var_names.size()];
				unsigned int col=0;
				while (ss1 >> value[col]) col++; 
			 	if( col!=var_names.size() ) {
			 	       	std::cout << "FATAL! unexpected number of columns!" << std::endl;
			 	       	std::cout << "... exiting" << std::endl;
                         		//for(unsigned int i=0; i<ncols; ++i)
			 	       	//                    std::cout << value[i] << " ";
			 	       	exit(1);
			 	}	

				for(unsigned int i=0; i<var_names.size(); ++i)	
				{
					aux_pars[var_names[i]]=value[i];
					//std::cout << var_names[i] << ": " << value[i] << " ";
				}
				std::cout << std::endl;
				aux_data.push_back(aux_pars);
			}				
		thisfile_aux.close();				
		}
		else
		{
			std::cout << "FATAL! unable to open file: " << infile_aux << std::endl; 
			std::cout << "... exiting" << std::endl;
			exit(1);
		}
		
	}

	for (unsigned int i=toy_id_low; i<toy_id_high; ++i)
	{	
		if(randomize_auxilliary_data) {
			// get new auxilliary best fits
			analysis.update_auxillary_data(aux_data[i]);
		}
		update_seed(seed_vals);
		point["toy_index"]=(double)i;
		run_toyfit(&rootfile, i, point, names);	
		points.push_back(point);

		if(i%flush_rate==0) {
			// flush result to file
			write_toy_fits(file, points);
			std::cout << "... " << i+1-toy_id_low << "/" << nfits << " fits done." << std::endl;
			points.clear();
		}

		// need to set up minimizer again
		update_minimizer(pars);
	
	}

	write_toy_fits(file, points);
	points.clear();

	rootfile.Close();
	file.close();

	// copy original data histograms back (remove toy data from analysis)
	analysis.restore_data_hists();

	// reset auxillary best fits
	analysis.reset_auxillary_data();
		

	return;
}

void NuFit::stats::run_toyfit(TFile *rootfile, int toy_index, std::map<std::string, double> &point, std::vector<std::string> &names)
{
	for (unsigned int j=0; j<names.size(); ++j)
	{
		// read from file
		std::string histname = names[j] + std::string("_mcsum_") + std::to_string(toy_index);
		TH3D* datahist = (TH3D*) rootfile->Get(histname.c_str());
		// and replace original datafile of this analysis
		// (inject toy data from file to analysis)
		analysis.set_hist(names[j], datahist);
		delete datahist;
	}

	// do we need to fit or simply calculate likelihoods
	double llh = 0.0;
	if(point.size() == analysis.get_npars())
	{
		// does not overwrite point with solution
		llh = get_llh_value(point);
		point.at("llh") = llh;
		point.at("fit_status") = 0.0; // No fit was run, so fit status is good :-)
	}
	else
		// overwrites point with solution	
		get_profile_llh_value(point, false, false); // treat as if this was a scan

	return;
}

void NuFit::stats::update_minimizer(const std::map<std::string, double> &pars)
{
	std::vector<std::string> par_names;
	analysis.get_par_names(par_names);
	std::map<std::string, unsigned int> parameters;
	analysis.get_par_names(parameters);

	// figure out what parameters need to be constant
	std::vector<unsigned int> indices;

	for (std::map<std::string, double>::const_iterator it=pars.begin(); it!=pars.end(); ++it)
	{
		if(parameters.count(it->first)==0) 
		{
			// parameters specified by user is not a model parameter:-(
			std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
			std::cout << "!!! the model has the following parameters:" << std::endl;
			for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
				std::cout << it2->first << std::endl;
			std::cout << "... exiting" << std::endl;
			exit(1);
		}
		indices.push_back(parameters.at(it->first));
	}

	min->Clear();
    for (unsigned int i=0; i<analysis.get_npars(); ++i)
    {
        // fix requested variables
        if(std::find(indices.begin(), indices.end(), i) != indices.end())
        { 
            min->SetFixedVariable(i, par_names[i].c_str(), pars.at(par_names[i])); 
            std::cout << "... treating variable: " << par_names[i] << " with index " << i << " as constant." << std::endl;
        }
        else
            min->SetLimitedVariable(i, par_names[i].c_str(), seeds[i], stepsizes[i], limits_low[i], limits_high[i]);
    }	
	return;
}

void NuFit::stats::update_seed(const std::map<std::string, double> &seed_vals)
{
        std::vector<std::string> par_names;
        analysis.get_par_names(par_names);
        std::map<std::string, unsigned int> parameters;
        analysis.get_par_names(parameters);

        // figure out what parameters need to be constant
        std::vector<unsigned int> indices;

        for (std::map<std::string, double>::const_iterator it=seed_vals.begin(); it!=seed_vals.end(); ++it)
        {
                if(parameters.count(it->first)==0)
                {
                        // parameters specified by user is not a model parameter:-(
                        std::cout << "!!! FATAL: parameter " << it->first << " not found. Are you sure this parameter is in the model?" << std::endl;
                        std::cout << "!!! the model has the following parameters:" << std::endl;
                        for (std::map<std::string, unsigned int>::iterator it2 = parameters.begin(); it2 != parameters.end(); ++it2 )
                                std::cout << it2->first << std::endl;
                        std::cout << "... exiting" << std::endl;
                        exit(1);
                }
                indices.push_back(parameters.at(it->first));
        }

	for (unsigned int i=0; i<analysis.get_npars(); ++i)
        {
                // update seeded variables
                if(std::find(indices.begin(), indices.end(), i) != indices.end())
                {
                        seeds[i] = seed_vals.at(par_names[i]);
                }
	}
	return;
}


double NuFit::stats::calculate_gof_toy(TFile *rootfile, int toy_index, std::map<std::string, double> &point, std::vector<std::string> &names)
{

        analysis.cache_data_hists();

        for (unsigned int j=0; j<names.size(); ++j)
        {
                // read from file
                std::string histname = names[j] + std::string("_mcsum_") + std::to_string(toy_index);
                TH3D* datahist = (TH3D*) rootfile->Get(histname.c_str());
                // and replace original datafile of this analysis
                // (inject toy data from file to analysis)
                analysis.set_hist(names[j], datahist);
                delete datahist;
        }

        double llh_gof = get_llh_value(point);

        analysis.restore_data_hists();

        return llh_gof;
}

void NuFit::stats::run_predicted_gof_bayes(std::string outfile, std::string infile, std::string infile_parameter_values, unsigned int toy_id_low, unsigned int toy_id_high)
{
        std::cout << std::endl;
        std::cout << "... requested simulation of posterior predicted deviance using toy-data." << std::endl;

        if(toy_id_high <= toy_id_low)
        {
                std::cout << "!!! FATAL:  toy_id_high <= toy_id_low. can't perform fits." << std::endl;
                std::cout << "... exiting" << std::endl;
        }

        // check if file exists
        if (std::ifstream(outfile))
        {
                std::cout << "!!! FATAL: File " << outfile << " already exists. Refuse to overwrite!" << std::endl;
                std::cout << "... exiting" << std::endl;
                exit(1);
        }

        //  check if can write to file
        std::ofstream file(outfile);
        file.precision(5);
        file.setf(std::ios_base::scientific);

        if (!file)
        {
                std::cout << "!!! FATAL: unable to write to file " << outfile << std::endl;
                exit(1);
        }

        // read infile with parameters
        std::string line;
        std::ifstream infile_par_vals(infile_parameter_values);
        std::map<std::string, double> point;
        std::vector<std::map<std::string, double>> points;
        std::vector<std::string> column_names;
        unsigned int ncols = 0;
        if (infile_par_vals.is_open())
        {
                getline(infile_par_vals, line); // assume first line has column names
                std::istringstream ss(line);
                std::string col_name;
                while (ss >> col_name)
                {
                        column_names.push_back(col_name);
                        point.insert(std::pair<std::string, double>(col_name, 0.0));
                        ncols++;
                }

                // check if number of parameters matches model expectations.
                // also do the names match?
                unsigned int npars_model = analysis.get_npars();
                if(column_names.size() != npars_model)
                {
                        std::cout << "FATAL! number of parameters required by model (" << npars_model << " parameters) does not match file content (" << column_names.size() << " parameters)" << std::endl;
                        std::cout << "... exiting" << std::endl;
                        exit(1);
                }

                std::vector<std::string> par_names;
                analysis.get_par_names(par_names);
                for(unsigned int i=0;i<column_names.size(); ++i) {
                        std::string par_name = column_names[i];
                        if(std::find(par_names.begin(), par_names.end(), par_name) == par_names.end())
                        {
                                // cat find par_name in model
                                std::cout << "FATAL! can't find parameter " << par_name << " in model parameters! " << std::endl;
                                std::cout << "... exiting" << std::endl;
                                exit(1);
                        }
                }
                // now read all the model parameters            
                std::vector<double> current_parameters;
                current_parameters.resize(ncols);
                while(getline(infile_par_vals, line))
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
                        for (unsigned int i=0; i<ncols; ++i)
                                point[column_names[i]]=current_parameters[i];

                        points.push_back(point);

                }

                /*
+               for(unsigned int i=0; i<points.size(); ++i)
+               {
+                       for (std::map<std::string, double>::const_iterator it = points[i].begin(); it != points[i].end(); ++it) std::cout << it->first << " " << it->second << " ";
+                       std::cout << "\n";
+               }
+               std::cout << std::endl;
+               */

                // read toy data from file
                std::vector<std::string> names = analysis.get_analysis_names();
                TFile rootfile(infile.c_str(), "READ");

                // and write resulting samples to outfile
                file << "Lsat(obs, theta)" << " " << "Lsat(rep, theta)";
                for (unsigned int i=0; i<column_names.size(); ++i)
                        file << " " << column_names[i];
                file << "\n";

                unsigned int ncalc=0;
                unsigned int ncalc_tot=toy_id_high - toy_id_low;
                for (unsigned int i=toy_id_low; i<toy_id_high; ++i)
                {
                        file << calculate_gof_toy(&rootfile, i, points[i], names) << " " << get_llh_value(points[i]);
                        for (unsigned int j=0; j<column_names.size(); ++j)
                                file << " " << points[i][column_names[j]];
                        file << "\n";
                        ncalc++;
                        if(ncalc%100==0)
                                std::cout << "... " << ncalc * 1.0 / ncalc_tot * 100 << "\% done. " << std::endl;
                }

                file << std::endl;

                rootfile.Close();
                file.close();




        }
        else
        {
                 std::cout << "FATAL! unable to open required file:\n" << infile_parameter_values << std::endl;
                 std::cout << "... exiting" << std::endl;
                 exit(1);
        }

        //outfile.close();
        infile_par_vals.close();

}

void NuFit::stats::get_histograms_from_samples(std::string outdir, std::map<std::string, double> simpars_fixed, std::vector<std::string> simpars_float, std::string infile_parameter_samples, unsigned int nsamples)
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

        for (unsigned int i=0; i<points.size(); ++i)
        {
                std::string toutfile=outdir+std::string("histograms_sample")+std::to_string(i)+std::string(".root");
                analysis.get_histograms(toutfile, points[i]);
        }

        return;
}




