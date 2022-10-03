#include "../../include/systematics/interpolated_par.h"

NuFit::interpolated_par::interpolated_par(std::string par_name_, std::vector<std::string> analysis_names_)
{
    par_name = par_name_;
    analysis_names = analysis_names_;
    basedir=std::string("/data/user/zzhang1/fit_pass2_sys/input/");

    // note that by default only those flavors will be adjusted that have corresponding simulation added
    // see derived class constructors (e.g. interpolated_par_domeff.cpp

    // define what components and flavors are effected by this systematics parameter
    // if derived implementation does not effect e.g. astro (like self veto)
    // then simply don't add astro in the derived constructor and astro corrections will always return 1.0 (skipped)
    components.push_back(std::string("Conv"));
    components.push_back(std::string("Prompt"));
    components.push_back(std::string("Astro"));
    return;
}

NuFit::interpolated_par::~interpolated_par() 
{ 
    for(std::unordered_map<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>>::iterator it=binfits.begin(); it!=binfits.end(); ++it)
    {
        for(std::unordered_map<std::string, NuFit::interpolated_sys *>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); )
        {
            // clean up memory
            (it->second).erase(it2++);
        }
    }
    std::cout << std::endl;
    std::cout << "... cleaned base class for interpolated parameters from memory." << std::endl;
}

double NuFit::interpolated_par::get_efficiency_correction(const double &x, const std::string &analysis_name, const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz)
{
    // flavor: NuE, NuMu, NuTau
    // component: Conv, Prompt, Astro, Muon
    
    // per bin fits for analysis
    std::unordered_map<std::string, NuFit::interpolated_sys *> &analysis = binfits[analysis_name];

    // first find flavor for analysis
    std::unordered_map<std::string, NuFit::interpolated_sys *>::const_iterator it = analysis.find(flavor);

    // if not found then this systematic does not effect this flavor 
    if (it == analysis.end()) return 1.0;
        // eventually need to work in effect on muon background
    else 
    {
        // need to find correction factor for component 
        // check if component is actually effected by this systematic
        if (std::find(components.begin(), components.end(), component) != components.end()) return analysis[flavor]->get_efficiency_correction(x, component, binx, biny, binz);
        else return 1.0;
    }
    
    return 1.0; // should never be reached
}

double NuFit::interpolated_par::get_efficiency_correction_error(const double &x, const std::string &analysis_name, const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz)
{
    // flavor: NuE, NuMu, NuTau
    // component: Conv, Prompt, Astro, Muon
    
    // per bin fits for analysis
    std::unordered_map<std::string, NuFit::interpolated_sys *> &analysis = binfits[analysis_name];

    // first find flavor for analysis
    std::unordered_map<std::string, NuFit::interpolated_sys *>::const_iterator it = analysis.find(flavor);

    // if not found then this systematic does not effect this flavor 
    if (it == analysis.end()) return 0.0;
        // eventually need to work in effect on muon background
    else 
    {
        // need to find correction factor for component 
        // check if component is actually effected by this systematic
        if (std::find(components.begin(), components.end(), component) != components.end()) return analysis[flavor]->get_efficiency_correction_error(x, component, binx, biny, binz);
        else return 0.0;
    }
    return 0.0; // should never be reached
}

double NuFit::interpolated_par::get_efficiency_correction_muon(const double &x, const std::string &analysis_name, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz)
{
    // this will currently map onto neutrino systematics (since we don't have muon background systematics simulation. will assume a soft spectrum for the computation
    // --> model after NuE Conventional
    
    // per bin fits for analysis
    std::unordered_map<std::string, NuFit::interpolated_sys *> &analysis = binfits[analysis_name];
    
    // check whether interpolated par effects muons
    std::unordered_map<std::string, NuFit::interpolated_sys *>::const_iterator it = analysis.find("Muon");
    
    // if not found then this systematic does not effect this flavor
    /*
        if (it == analysis.end()) return 1.0;
                // eventually need to work in effect on muon background
        else
        {
                 // analysis["Muon"] is of type NuFit::neutrino_input (since we are using neutrino systematics)
                return analysis["Muon"]->get_efficiency_correction(x, "Conv", binx, biny, binz); 
        }
    */
    return 1.0; // should never be reached
}


