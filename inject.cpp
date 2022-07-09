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
	NuFit::helpers::par_options astro_norm("astro_norm", 1.58, 0.01, 0., 10.);
	NuFit::helpers::par_options astro_index("astro_index", 2.53, 0.01, 0, 10.0);
	//NuFit::helpers::par_options energy_cut("energy_cut", 6, 0.01, 3.0, 7.0);	
    NuFit::helpers::par_options muon_norm("muon_norm", 1.5, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.1, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 0.0, 0.1, 0.0, 100.0);
	//NuFit::helpers::par_options prompt_norm("prompt_norm", 0.0, 0.0, 0.0, 100.0);
	NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.01, -0.5, 0.5);
    NuFit::helpers::par_options dom_eff("dom_efficiency", 1.0, 0.01, 0.0, 10.0);
    NuFit::helpers::par_options scattering("scattering", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options absorption("absorption", 1.0, 0.01, 0.0, 10.0);
    NuFit::helpers::par_options holeicep0("holeicep0", 0.0, 0.01, -1.0, 2.0);
	//NuFit::helpers::par_options holeicep0("holeicep0", 0.0, 0.00, -1.0, 2.0);
	NuFit::helpers::par_options holeicep1("holeicep1", 0.0, 0.01, -0.2, 0.2);
	//NuFit::helpers::par_options holeicep1("holeicep1", 0.0, 0.00, -0.2, 0.2);
    NuFit::helpers::par_options selfveto("selfveto", 750, 50, 5.0, 3000);
	//NuFit::helpers::par_options selfveto("selfveto", 1000, 0, 5.0, 3000);
	//NuFit::helpers::par_options oldsv("oldsv", 700, 50, 100, 5000);
    //NuFit::helpers::par_options cosmicray("cosmicray", 0, 0.01, -1, 2);
    //NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 0, 0.01, -1, 2);

	/** .. package everything */
	std::map<std::string, NuFit::helpers::par_options> options;
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm.name, astro_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index.name, astro_index));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(energy_cut.name, energy_cut));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(conv_norm.name, conv_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm.name, muon_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(prompt_norm.name, prompt_norm));	
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(cr_index.name, cr_index));
    options.insert(std::pair<std::string, NuFit::helpers::par_options>(dom_eff.name, dom_eff));
    options.insert(std::pair<std::string, NuFit::helpers::par_options>(scattering.name, scattering));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(absorption.name, absorption));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep0.name, holeicep0));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep1.name, holeicep1));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(selfveto.name, selfveto));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(oldsv.name, oldsv));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(cosmicray.name, cosmicray));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(hadronicinteraction.name, hadronicinteraction));
	/** .. and send to minimizer */
	min.set_options(options);

	/** write seed histograms */
    /*
	std::map<std::string, double> pars;
	pars.insert(std::pair<std::string, double>(muon_norm.name, muon_norm.seed));
	pars.insert(std::pair<std::string, double>(conv_norm.name, conv_norm.seed));
	pars.insert(std::pair<std::string, double>(prompt_norm.name, prompt_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_norm.name, astro_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_index.name, astro_index.seed));
	//pars.insert(std::pair<std::string, double>(energy_cut.name, energy_cut.seed));
	pars.insert(std::pair<std::string, double>(cr_index.name, cr_index.seed));
    pars.insert(std::pair<std::string, double>(dom_eff.name, dom_eff.seed));
    pars.insert(std::pair<std::string, double>(scattering.name, scattering.seed));
	pars.insert(std::pair<std::string, double>(absorption.name, absorption.seed));
	pars.insert(std::pair<std::string, double>(holeicep0.name, holeicep0.seed));
	pars.insert(std::pair<std::string, double>(holeicep1.name, holeicep1.seed));
	pars.insert(std::pair<std::string, double>(selfveto.name, selfveto.seed));
	//pars.insert(std::pair<std::string, double>(oldsv.name, oldsv.seed));
    min.set_tolerance(1);
    */
    int nsteps = 29;
    double s_min = 100;
    double s_max = 2900;

    //scan
    NuFit::helpers::scan_options s("oldsv", nsteps, s_min, s_max);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s.name, s));

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
    asimov_point.insert(std::pair<std::string, double>(selfveto.name, 1000.0));
    //asimov_point.insert(std::pair<std::string, double>(oldsv.name, 1000.0));
    //asimov_point.insert(std::pair<std::string, double>(cosmicray.name, 1.0));
    //asimov_point.insert(std::pair<std::string, double>(hadronicinteraction.name, 1.0));

    min.scan_profile_llh_asimov(outdir+std::string("outfile_profile_llh_2d")+std::string(".txt"), scan_profile, outdir+std::string("injected_hist")+std::string(".root"), asimov_point);

       
        return 0;
}
