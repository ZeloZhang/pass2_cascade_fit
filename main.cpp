#include <cstdlib>
#include <string>
#include <map>
#include <TMath.h>
#include "include/stats.h"
#include "include/analysis.h"
#include "include/helpers.h"
#include "include/models/astro_model_single_plaw_wcutoff.h"
#include "include/bootstrap/toymc.h"



int main(int argc, char **argv)
{

	std::string outdir=std::string("/data/user/zzhang1/fit_pass2_sys/output/");

	if (argc!=5) {
		std::cout << "expect four integer arguments: jobid[1,2] and njobs[1,2]. quitting!" << std::endl;
		return 1;
	}

	// get arguments relevant for profile scan
	int jobid1 = std::stoi(argv[1]);
	int njobs1 = std::stoi(argv[2]);
	int jobid2 = std::stoi(argv[3]);
	int njobs2 = std::stoi(argv[4]);
    //int random_seed = jobid;
//        int njobs = 1;

	/** create analysis. */
	NuFit::analysis wrapper;

	/** need to specifically create the analysis. */
	/** this function contains all the analysis specific details and can be edited in ./src/analysis.cpp */

	wrapper.create(); 

	/** get occasional prints of minimization progress from inside likelihood function */
	//wrapper.set_verbosity(true);

	/** use stats class to analyze NuFit::analysis */
	NuFit::stats min(wrapper);
	/** option to change model (NuFit::analysis uses single powerlaw as default) */
	//NuFit::astro_model_single_plaw_wcutoff *new_astro_model = new NuFit::astro_model_single_plaw_wcutoff();
	//min.change_astro_model(new_astro_model);

	
	/** specify seed, stepsize, limits for each parameter */
	/** parameter names need to match specification in model class, but parameter ordering does not matter! */
	/** e.g. ./src/astro_model_single_plaw.cpp + ./src/model_base.cpp */
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */
	NuFit::helpers::par_options astro_norm("astro_norm", 1.58, 0.001, 0., 10.);
	NuFit::helpers::par_options astro_index("astro_index", 2.53, 0.01, 0, 10.0);
	//NuFit::helpers::par_options energy_cut("energy_cut", 6, 0.01, 3.0, 7.0);	
    NuFit::helpers::par_options muon_norm("muon_norm", 1.5, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.1, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 0.0, 0.1, 0.0, 100.0);
    /*
	NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.01, -0.5, 0.5);
    NuFit::helpers::par_options dom_eff("dom_efficiency", 1.0, 0.01, 0.0, 10.0);
    NuFit::helpers::par_options scattering("scattering", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options absorption("absorption", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options holeicep0("holeicep0", 0.0, 0.05, -2.0, 1.0);
	NuFit::helpers::par_options holeicep1("holeicep1", 0.0, 0.01, -0.2, 0.2);
	NuFit::helpers::par_options selfveto("selfveto", 1000, 10, 50.0, 3000);
	//NuFit::helpers::par_options cosmicray("cosmicray", 0.0, 0.01, -2, 3);
	NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 0.0, 0.01, -2, 3);
    */

	/** .. package everything */
	std::map<std::string, NuFit::helpers::par_options> options;
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm.name, astro_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index.name, astro_index));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(energy_cut.name, energy_cut));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(conv_norm.name, conv_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm.name, muon_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(prompt_norm.name, prompt_norm));	
    /*
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(cr_index.name, cr_index));
    options.insert(std::pair<std::string, NuFit::helpers::par_options>(dom_eff.name, dom_eff));
    options.insert(std::pair<std::string, NuFit::helpers::par_options>(scattering.name, scattering));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(absorption.name, absorption));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep0.name, holeicep0));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep1.name, holeicep1));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(selfveto.name, selfveto));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(cosmicray.name, cosmicray));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(hadronicinteraction.name, hadronicinteraction));
    */
	/** .. and send to minimizer */
	min.set_options(options);

	/** write seed histograms */
	std::map<std::string, double> pars;
	pars.insert(std::pair<std::string, double>(muon_norm.name, muon_norm.seed));
	pars.insert(std::pair<std::string, double>(conv_norm.name, conv_norm.seed));
	pars.insert(std::pair<std::string, double>(prompt_norm.name, prompt_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_norm.name, astro_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_index.name, astro_index.seed));
	//pars.insert(std::pair<std::string, double>(energy_cut.name, energy_cut.seed));
    /*
	pars.insert(std::pair<std::string, double>(cr_index.name, cr_index.seed));
    pars.insert(std::pair<std::string, double>(dom_eff.name, dom_eff.seed));
    pars.insert(std::pair<std::string, double>(scattering.name, scattering.seed));
	pars.insert(std::pair<std::string, double>(absorption.name, absorption.seed));
	pars.insert(std::pair<std::string, double>(holeicep0.name, holeicep0.seed));
	pars.insert(std::pair<std::string, double>(holeicep1.name, holeicep1.seed));
	pars.insert(std::pair<std::string, double>(selfveto.name, selfveto.seed));
	//pars.insert(std::pair<std::string, double>(cosmicray.name, cosmicray.seed));
	pars.insert(std::pair<std::string, double>(hadronicinteraction.name, hadronicinteraction.seed));
    */
    min.set_tolerance(1);

    min.fit(false); // ** if true -> get profile LLH errors after minimization from ROOT Minuit2
    //min.fit(true);

    std::string outfile("./output/hists_seed.root");
	wrapper.get_histograms(outfile,pars);
    min.get_bestpars(pars);
	outfile = std::string("./output/hists_fit.root");
	wrapper.get_histograms(outfile, pars);

    /*
    int nsteps_selfveto = 29;
    double selfveto_min = 100;
    double selfveto_max = 2900;
    double ds = (selfveto_max-selfveto_min) / (nsteps_selfveto-1);
    if (nsteps_selfveto % njobs1 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job1 = nsteps_selfveto / njobs1;
    double selfveto_min_tj = selfveto_min + jobid1 * nsteps_job1 * ds;
    double selfveto_max_tj = selfveto_min_tj + (nsteps_job1-1) * ds;

    //scan
    NuFit::helpers::scan_options s_selfveto("selfveto", nsteps_job1, selfveto_min_tj, selfveto_max_tj);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_selfveto.name, s_selfveto));

	min.scan_profile_llh(outdir+std::string("outfile_profile_llh_2d_part_")+std::to_string(jobid1)+std::string(".txt"), scan_profile);
    
    */ 
    //min.get_bestpars(pars);

    /*
    int nsteps_index = 13;
    double index_min = 2.2;
    double index_max = 2.8;
    double ds = (index_max-index_min) / (nsteps_index-1);
    if (nsteps_index % njobs1 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job1 = nsteps_index / njobs1;
    double index_min_tj = index_min + jobid1 * nsteps_job1 * ds;
    double index_max_tj = index_min_tj + (nsteps_job1-1) * ds;

    int nsteps_norm = 16;
    double norm_min = 1.0;
    double norm_max = 2.5;
    ds = (norm_max-norm_min)/(nsteps_norm-1);
    if (nsteps_norm % njobs2 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job2 = nsteps_norm / njobs2;
    double norm_min_tj = norm_min + jobid2 * nsteps_job2 * ds;
    double norm_max_tj = norm_min_tj + (nsteps_job2-1) * ds;

    //scan
    NuFit::helpers::scan_options s_astro_norm("astro_norm", nsteps_job2, norm_min_tj, norm_max_tj);
    NuFit::helpers::scan_options s_astro_index("astro_index", nsteps_job1, index_min_tj, index_max_tj);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm.name, s_astro_norm));
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_index.name, s_astro_index));

    std::map<std::string, double> asimov_point;
    asimov_point.insert(std::pair<std::string, double>(muon_norm.name, 1.5));
    asimov_point.insert(std::pair<std::string, double>(conv_norm.name, 1.0689));
    asimov_point.insert(std::pair<std::string, double>(prompt_norm.name, 0));
    asimov_point.insert(std::pair<std::string, double>(astro_norm.name, 1.66));
    asimov_point.insert(std::pair<std::string, double>(astro_index.name, 2.53));
    asimov_point.insert(std::pair<std::string, double>(dom_eff.name, 1.0));
    asimov_point.insert(std::pair<std::string, double>(scattering.name, 1.0));
    asimov_point.insert(std::pair<std::string, double>(absorption.name, 1.0));
    asimov_point.insert(std::pair<std::string, double>(holeicep0.name, 0.0));
    asimov_point.insert(std::pair<std::string, double>(holeicep1.name, 0.0));
    asimov_point.insert(std::pair<std::string, double>(cr_index.name, 0.0));
    //asimov_point.insert(std::pair<std::string, double>(cosmicray.name, 0.0));
    asimov_point.insert(std::pair<std::string, double>(hadronicinteraction.name, 0.0));
    //min.scan_profile_llh_asimov(outdir+std::string("outfile_profile_llh_2d_part_")+std::to_string(jobid1)+std::string("_")+std::to_string(jobid2)+std::string(".txt"), scan_profile, outdir+std::string("injected_hist_")+std::to_string(jobid1)+std::string("_")+std::to_string(jobid2)+std::string(".root"), asimov_point);
    min.scan_profile_llh_asimov(outdir+std::string("outfile_profile_llh_2d_part_")+std::to_string(jobid1)+std::string("_")+std::to_string(jobid2)+std::string(".txt"), scan_profile, outdir+std::string("injected_hist")+std::string(".root"));
    //min.get_bestpars(pars);
    */

    /** write best fit histograms */
	//std::string outfile = outdir+std::string("hists_fit.root");
    //wrapper.get_histograms(outfile, pars);
	
	/** example of how to do a profile llh scan */
    /*
	int nsteps_index = 5;	
	double index_min = 2.3;
	double index_max = 2.6;

	double ds = (index_max-index_min) / (nsteps_index-1);

	if (nsteps_index % njobs != 0) {
		std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
		return 1;
	}

	int nsteps_job = nsteps_index / njobs;
        */
	//double index_min_tj = index_min + jobid * nsteps_job * ds;
	//double index_max_tj = index_min_tj + (nsteps_job-1) * ds;
	
        //scan
        /*
	NuFit::helpers::scan_options s_astro_norm("astro_norm", 20, 1.5, 1.7);
	//NuFit::helpers::scan_options s_astro_index("astro_index", nsteps_job, index_min_tj, index_max_tj);
	NuFit::helpers::scan_options s_astro_index("astro_index", 20, 2.3, 2.6);

	//std::cout << ds << " " << index_min_tj << " " << index_max_tj << " " << nsteps_job << std::endl;
	
	std::map<std::string, NuFit::helpers::scan_options> scan;
	scan.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm.name, s_astro_norm));
	scan.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_index.name, s_astro_index));
	

	min.set_flush_rate(1); // set frequency of file dumps of scan results during scan

	//min.scan_profile_llh(outdir+std::string("outfile_profile_llh_part")+std::to_string(jobid)+std::string(".txt"), scan);
	//min.scan_profile_llh(outdir+std::string("outfile_profile_llh_2d_part")+std::string(".txt"), scan);
        
                        
        std::map<std::string, double> seed;
        seed.insert(std::pair<std::string, double>(muon_norm.name,1.66672));
        seed.insert(std::pair<std::string, double>(muon_norm1.name,1));
        seed.insert(std::pair<std::string, double>(conv_norm.name,1.03438));
        seed.insert(std::pair<std::string, double>(prompt_norm.name,5.02903e-5));
	min.scan_llh(outdir+std::string("outfile_llh_2d_part")+std::string(".txt"),scan,seed);
        */
	
	//scan_profile
	/*
	NuFit::helpers::scan_options s_astro_norm_profile("astro_norm", 10, 1.5, 1.7);
	NuFit::helpers::scan_options s_astro_index_profile("astro_index", 10, 2.3, 2.6);

	std::map<std::string, NuFit::helpers::scan_options> scan_profile;
	scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm_profile.name, s_astro_norm_profile));
	scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_index_profile.name, s_astro_index_profile));
	
	min.scan_profile_llh(outdir+std::string("outfile_profile_llh_2d_part")+std::string(".txt"), scan_profile);
        */
 /*       
        //toymc
        std::string outfile = outdir+std::string("hists_fit.root");
        wrapper.get_histograms(outfile, pars);


        pars[muon_norm.name]=1.37842;
        pars[conv_norm.name]=1.0;
        pars[prompt_norm.name]=1.0;
        pars[astro_norm.name]=1.0;
        pars[astro_index.name]=2.0;
        
        
        pars[dom_eff.name]=1.1003;
        //pars[scattering.name]=1.07412;
        //pars[absorption.name]=0.989128;
        pars[cr_index.name]=0.0548299;
        //pars[holeice.name]=1.74587;
        //pars[muon_norm1.name]=0.7497; // mlb
        
      
        
        //double random_seed = 1;
        NuFit::toymc sim(wrapper, random_seed);

        std::string outfile_toy = outdir+std::string("hists_toymc_job")+std::to_string(random_seed)+std::string(".root");
        std::string outfile_aux = outdir+std::string("outfile_aux_job")+std::to_string(random_seed)+std::string(".txt");

        unsigned int nsamples = 5;
        std::cout<< pars.size();
        sim.run_simulation_waux(outfile_toy, outfile_aux, pars, nsamples);
        //sim.run_simulation(outfile_toy,pars,nsamples);

        std::map<std::string, double> point2;
        //point2.insert(std::pair<std::string, double>("astro_norm", 0.9));
        //point2.insert(std::pair<std::string, double>("astro_index", 2.5));
        min.set_flush_rate(1);
        
        
        std::string outfile_toyfits = outdir+std::string("outfile_toy_from_bestfit_job")+std::to_string(random_seed)+std::string("_promptonly.txt");
        pars[prompt_norm.name]=0.0;
        pars[astro_norm.name]=2.0;
        pars[astro_index.name]=3.0;
        std::cout<<pars.size();
        min.run_toyfits(outfile_toy, outfile_toyfits, point2, pars, 0, nsamples, true, outfile_aux);

        //std::string outfile_toyfits_constrain = outdir+std::string("outfile_toy_from_bestfit_job")+std::to_string(random_seed)+std::string("_promptonly_constrain.txt");
        //std::map<std::string, double> point3;
        //point3.insert(std::pair<std::string, double>("astro_norm", 0.0));
        //min.run_toyfits(outfile_toy, outfile_toyfits_constrain, point3, pars, 0, nsamples, true, outfile_aux);
   */     
        /*
        pars[prompt_norm.name]=0.0;
        pars[astro_norm.name]=1.0;
        pars[astro_index.name]=3.0;
        std::string outfile_toyfits = outdir+std::string("outfile_toy_from_bestfit_job")+std::to_string(random_seed)+std::string("_promptonly_diffseed.txt");
        min.run_toyfits(outfile_toy, outfile_toyfits, point2, pars, 0, nsamples, false, outfile_aux);

        std::string outfile_toyfits_constrain = outdir+std::string("outfile_toy_from_bestfit_job")+std::to_string(random_seed)+std::string("_promptonly_constrain.txt");
        std::map<std::string, double> point3;
        point3.insert(std::pair<std::string, double>("astro_norm", 0.0));
        min.run_toyfits(outfile_toy, outfile_toyfits_constrain, point3, pars, 0, nsamples, false, outfile_aux);
        */
        /*
        std::map<std::string, double> point3;
        min.set_flush_rate(1);
        std::string outfile_toyfits2 = outdir+std::string("outfile_toy_from_bestfit_job")+std::to_string(random_seed)+std::string("_wastro.txt");
        min.run_toyfits(outfile_toy, outfile_toyfits2, point3, pars, 0, nsamples, true, outfile_aux);
        */
       
        return 0;
}
