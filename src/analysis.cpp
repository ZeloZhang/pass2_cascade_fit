#include "../include/analysis.h"

NuFit::analysis::analysis() {
    verbose_llh = false;    
    n_llh_evals = 0;
}

NuFit::analysis::~analysis() {
    delete model;
}

void NuFit::analysis::create()
{
    std::cout << "... start wrapping analyses" << std::endl;
    std::cout << std::endl;

    // here goes the analysis code (parsing of data and mc)

    // define cascade signal sample
    std::vector<double> binsx; // logE
    for (unsigned int i=0; i<23; ++i)    
    //for (unsigned int i=4; i<23; ++i)    
    //for (unsigned int i=0; i<8; ++i)    
    //for (unsigned int i=4; i<8; ++i)    
    //for (unsigned int i=4; i<11; ++i)    
    //for (unsigned int i=7; i<23; ++i)    
    //for (unsigned int i=0; i<11; ++i)    
        binsx.push_back(2.6 + 0.2 * i);

    std::vector<double> binsy; // cosz (2 bins: Northern Sky, Southern Sky)
    //for (unsigned int i=0; i<11; ++i)
    //    binsy.push_back(-1.0+0.2*i);

    binsy.push_back(-2);
    binsy.push_back(0.2);
    binsy.push_back(0.6);
    binsy.push_back(2);

    std::vector<double> binsz; // ra (1 bin)
    binsz.push_back(-10.);
    binsz.push_back(10.);

    std::string defdir("/data/user/zzhang1/fit_pass2_sys/input/");
    std::string dir_baseline=defdir+std::string("baseline/");

    NuFit::hists analysis1("cascade_all", binsx, binsy, binsz); 
    std::string name_nue("nue_cascade.txt");
    std::string name_numu("numu_cascade.txt");
    std::string name_nutau("nutau_cascade.txt");
    std::string name_mu("mgun_cascade.txt"); 
    std::string name_data("data/data_cascade.txt");

    analysis1.read(dir_baseline+name_nue, dir_baseline+name_numu, dir_baseline+name_nutau, dir_baseline+name_mu, defdir+name_data);

    /*
    std::vector<double> binsx12; // logE
    for (unsigned int i=11; i<23; ++i)    
        binsx12.push_back(2.6 + 0.2 * i);
    std::vector<double> binsy12; // cosz (2 bins: Northern Sky, Southern Sky)
    binsy12.push_back(-1.0);
    binsy12.push_back(0.2);
    binsy12.push_back(0.6);
    binsy12.push_back(1.0);

    NuFit::hists analysis12("cascade_all_he1", binsx12, binsy12, binsz); 
    analysis12.read(dir_baseline+name_nue, dir_baseline+name_numu, dir_baseline+name_nutau, dir_baseline+name_mu, defdir+name_data);
    */

    // define muon background sample
    std::vector<double> binsx2; // logE
    binsx2.push_back(2.6);
    //binsx2.push_back(3.2);
    //binsx2.push_back(4);
    binsx2.push_back(4.778);

    std::vector<double> binsy2; // cosz (1 bin: all-sky)
    binsy2.push_back(-2.0);
    binsy2.push_back(2.0);

    std::vector<double> binsz2; // ra (1 bin)
    binsz2.push_back(-10.);
    binsz2.push_back(10.);

    NuFit::hists analysis2("muon", binsx2, binsy2, binsz2); 
    std::string name_nue2("nue_muon.txt");
    std::string name_numu2("numu_muon.txt");
    std::string name_nutau2("nutau_muon.txt");
    std::string name_mu2("mgun_muon.txt");
    std::string name_data2("data/data_muon.txt");
    analysis2.read(dir_baseline+name_nue2, dir_baseline+name_numu2, dir_baseline+name_nutau2, dir_baseline+name_mu2, defdir+name_data2);
        
    // define numu control sample
    std::vector<double> binsx3;
    for (unsigned int i=0; i<12; ++i)
    //for (unsigned int i=4; i<11; ++i)
         binsx3.push_back(2.6 + 0.2 * i);

    NuFit::hists analysis3("hybrid", binsx3, binsy2, binsz);
    std::string name_nue3("nue_hybrid.txt");
    std::string name_numu3("numu_hybrid.txt");
    std::string name_nutau3("nutau_hybrid.txt");
    std::string name_mu3("mgun_hybrid.txt");
    std::string name_data3("data/data_hybrid.txt");
    analysis3.read(dir_baseline+name_nue3, dir_baseline+name_numu3, dir_baseline+name_nutau3, dir_baseline+name_mu3, defdir+name_data3);

    std::vector<NuFit::hists> analyses;
    analyses.push_back(analysis1);
    //analyses.push_back(analysis12);
    analyses.push_back(analysis2);
    analyses.push_back(analysis3);

    // create astro model. then create base model
    NuFit::astro_model_single_plaw *astro = new astro_model_single_plaw();
    //NuFit::astro_model_plaw_singlep *astro = new astro_model_plaw_singlep();
    //base_model needs to know the input data as well as the astro model
    //NuFit::model_base *mymodel = new model_base(analyses, astro);
    //model = mymodel;
    
    // comment out code below if you don't need systematics        
    std::vector<std::string> analysis_names;
    std::map<std::string, NuFit::hists*> map_analyses;
    for (unsigned int i=0; i<analyses.size(); ++i) 
    {
        analysis_names.push_back(analyses[i].name);
        map_analyses.insert(std::pair<std::string, NuFit::hists*>(analyses[i].name, &(analyses[i])));
    }

    std::map<std::string, NuFit::interpolated_par *> systematics;
    
    // start with dom efficiency
    std::string sysname_eff("dom_efficiency");
    NuFit::interpolated_par_domeff *domeff = new interpolated_par_domeff(sysname_eff, analysis_names, map_analyses);

    // scattering
    std::string sysname_scat("scattering");
    NuFit::interpolated_par_scattering *scattering = new interpolated_par_scattering(sysname_scat, analysis_names, map_analyses);

    // absorption
    std::string sysname_abs("absorption");
    NuFit::interpolated_par_absorption *absorption = new interpolated_par_absorption(sysname_abs, analysis_names, map_analyses);
        
    // holeicep0
    std::string sysname_p0("holeicep0");
    NuFit::interpolated_par_p0 *holeicep0 = new interpolated_par_p0(sysname_p0, analysis_names, map_analyses);

    // holeicep1
    std::string sysname_p1("holeicep1");
    NuFit::interpolated_par_p1 *holeicep1 = new interpolated_par_p1(sysname_p1, analysis_names, map_analyses);

    // selfveto
    std::string sysname_sv("selfveto");
    NuFit::interpolated_par_selfveto *selfveto = new interpolated_par_selfveto(sysname_sv, analysis_names, map_analyses);

    
    // hadronic interaction
    std::string sysname_hadronicinteraction("hadronicinteraction");
    NuFit::interpolated_par_hadronicinteraction *hadronicinteraction = new interpolated_par_hadronicinteraction(sysname_hadronicinteraction, analysis_names, map_analyses);

    // cosmicray
    //std::string sysname_cosmicray("cosmicray");
    //NuFit::interpolated_par_cosmicray *cosmicray = new interpolated_par_cosmicray(sysname_cosmicray, analysis_names, map_analyses);

    // all objects in this map will be destroyed by model_base_sys later on
    systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_eff, domeff));
    systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_scat, scattering));
    systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_abs, absorption));
    systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_p0, holeicep0));
    systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_p1, holeicep1));
    systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_sv, selfveto));
    systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_hadronicinteraction, hadronicinteraction));
    //systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_cosmicray, cosmicray));
    
    NuFit::model_base_sys *mymodel = new model_base_sys(analyses, astro, systematics, 10);
    model = mymodel; 

    // end of analysis code
    std::cout << std::endl;
    std::cout << "... wrapping done" << std::endl; 
}

// the functions below interface to the model class

void NuFit::analysis::change_astro_model(NuFit::astro_model_base *astro)
{
    model->change_astro_model(astro);
    return;
}

double NuFit::analysis::get_likelihood(const double *pars)
{
    // this is being minimized when stats class is used (interface to ROOT Minuit2)
    double neglogl = model->likelihood_say(pars);
    //double neglogl = model->likelihood(pars);
    ++n_llh_evals;

    // add verbosity if requested
    if(verbose_llh) {
        if(n_llh_evals%50==0 || n_llh_evals==1) {
            std::cout << "... evaluated llh " << n_llh_evals << " times already." << std::endl;
            std::cout << "f( ";
            for(unsigned int i=0; i<get_npars(); ++i)
                std::cout << pars[i] << " ";
            std::cout << ") = " << neglogl << std::endl;
            std::cout << std::endl;
        }
        
    }

    // factor of 2 for Wilk's theorem
    return 2.0 * neglogl;
}

double NuFit::analysis::get_likelihood_gof(const double *pars)
{
    // factor of 2 for Wilk's theorem
    return 2.0 * model->likelihood_gof_say(pars);
    //return 2.0 * model->likelihood_gof(pars);
}

double NuFit::analysis::get_likelihood_abs(const double *pars)
{
    // factor of 2 for Wilk's theorem
    return 2.0 * model->likelihood_abs_say(pars);
    //return 2.0 * model->likelihood_abs(pars);
}

/*
double NuFit::analysis::get_lnprob(boost::python::numeric::array pars_)
{
    // assume that pars is of length npars
    unsigned int npars = get_npars();    
    double pars[npars];
    for (unsigned int i=0; i<npars; ++i)
        pars[i] = boost::python::extract<double>(pars_[i]);

    // make log-likelihood positive
    return (-1.0) * model->likelihood_abs_say(pars); 
    //return (-1.0) * model->likelihood_abs(pars); 
}
*/

unsigned int NuFit::analysis::get_n_llh_evals() 
{
    return n_llh_evals;
}

unsigned int NuFit::analysis::get_npars()
{
    return model->get_npars();
}

void NuFit::analysis::get_par_names(std::vector<std::string> &names)
{
    model->get_par_names(names);
    return;
}

void NuFit::analysis::get_par_names(std::map<std::string, unsigned int> &names)
{
    model->get_par_names(names);
    return;
}

void NuFit::analysis::get_histograms(std::string outfile, std::map<std::string, double> &pars)
{
    model->get_histograms(outfile, pars);
    return;
}

std::vector<TH3D*> NuFit::analysis::get_hist_mcsum(std::map<std::string, double> &pars)
{ 
    return model->get_hist_mcsum(pars);
}

void NuFit::analysis::reset_n_llh_evals() 
{
    n_llh_evals=0;
    return;
}

void NuFit::analysis::set_verbosity(bool flag=false) {
    verbose_llh=flag;
    return;
}

std::vector<std::string> NuFit::analysis::get_analysis_names()
{
    return model->get_analysis_names();
}

void NuFit::analysis::set_hist(std::string analysis, TH3D* hist)
{
    model -> set_hist(analysis, hist);
    return;
}

void NuFit::analysis::cache_data_hists()
{
    model -> cache_data_hists();
    return;
}

void NuFit::analysis::restore_data_hists()
{
    model -> restore_data_hists();
    return;
}

void NuFit::analysis::update_auxillary_data(std::map<std::string, double> &pars)
{
    model -> update_auxillary_data(pars);
    return;
}

void NuFit::analysis::reset_auxillary_data()
{
    model -> reset_auxillary_data();
    return;
}

/*
#include <boost/python.hpp>
using namespace boost::python;
BOOST_PYTHON_MODULE(analysis)
{
    // expose functions create
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
        class_<NuFit::analysis>("analysis", init<>())    
                .def("create", &NuFit::analysis::create)
        .def("get_lnprob", &NuFit::analysis::get_lnprob)
        ;
}
*/
