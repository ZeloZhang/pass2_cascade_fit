#include "../../include/systematics/interpolated_par_scattering.h"

NuFit::interpolated_par_scattering::interpolated_par_scattering(std::string par_name_, std::vector<std::string> analysis_names_, const std::map<std::string, NuFit::hists*> &map_analyses): interpolated_par(par_name_, analysis_names_)

{
    // need to know what analyses are implemented -> check /src/analysis.cpp

    std::string dir=basedir+std::string("sys/txt/");
    std::string baselinedir=basedir+std::string("baseline/");

    std::string dir_muon=basedir+std::string("sys/txt/");
    std::string baselinedir_muon=basedir+std::string("baseline/");

    // cascade event selection first
    std::string cascade("cascade_all");
    if (std::find(analysis_names.begin(), analysis_names.end(), cascade) != analysis_names.end())
    {    
        std::vector<double> binsx = map_analyses.at(cascade)->get_binsx();
        std::vector<double> binsy = map_analyses.at(cascade)->get_binsy();
        std::vector<double> binsz = map_analyses.at(cascade)->get_binsz();

        NuFit::interpolated_sys *nue_scattering_cascade = new interpolated_sys(par_name, cascade, "NuE", binsx, binsy, binsz);
        nue_scattering_cascade -> add_simulated_point(1.0, baselinedir + std::string("nue_cascade.txt"), true); // this is baseline hist
        nue_scattering_cascade -> add_simulated_point(1.1, dir + std::string("s+.10/nue_cascade.txt"));
        nue_scattering_cascade -> create_correction_functions();

        NuFit::interpolated_sys *numu_scattering_cascade = new interpolated_sys(par_name, cascade, "NuMu", binsx, binsy, binsz);  
        numu_scattering_cascade -> add_simulated_point(1.0, baselinedir + std::string("numu_cascade.txt"), true); // this is baseline hist
        numu_scattering_cascade -> add_simulated_point(1.1, dir + std::string("s+.10/numu_cascade.txt"));
        numu_scattering_cascade -> create_correction_functions();

        NuFit::interpolated_sys *nutau_scattering_cascade = new interpolated_sys(par_name, cascade, "NuTau", binsx, binsy, binsz);
        nutau_scattering_cascade -> add_simulated_point(1.0, baselinedir + std::string("nutau_cascade.txt"), true); // this is baseline hist
        nutau_scattering_cascade -> add_simulated_point(1.1, dir + std::string("s+.10/nutau_cascade.txt"));
        nutau_scattering_cascade -> create_correction_functions();

        // map from flavor to correction
        std::unordered_map<std::string, NuFit::interpolated_sys *> cascade_corrections;
        cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuE", nue_scattering_cascade));
        cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuMu", numu_scattering_cascade));
        cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuTau", nutau_scattering_cascade));
        //cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("Muon", nue_scattering_cascade));

        // and now add
        binfits.insert(std::pair<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>>(cascade, cascade_corrections));
    }


    // muon event selection
    std::string muon("muon");
    if (std::find(analysis_names.begin(), analysis_names.end(), muon) != analysis_names.end())
    {    
        std::vector<double> binsx = map_analyses.at(muon)->get_binsx();
        std::vector<double> binsy = map_analyses.at(muon)->get_binsy();
        std::vector<double> binsz = map_analyses.at(muon)->get_binsz();

        NuFit::interpolated_sys *nue_scattering_muon = new interpolated_sys(par_name, muon, "NuE", binsx, binsy, binsz);
        nue_scattering_muon -> add_simulated_point(1.0, baselinedir_muon + std::string("nue_muon.txt"), true); // this is baseline hist
        nue_scattering_muon -> add_simulated_point(1.1, dir_muon + std::string("s+.10/nue_muon.txt"));
        nue_scattering_muon -> create_correction_functions();

        NuFit::interpolated_sys *numu_scattering_muon = new interpolated_sys(par_name, muon, "NuMu", binsx, binsy, binsz); 
        numu_scattering_muon -> add_simulated_point(1.0, baselinedir_muon + std::string("numu_muon.txt"), true); // this is baseline hist
        numu_scattering_muon -> add_simulated_point(1.1, dir_muon + std::string("s+.10/numu_muon.txt"));
        numu_scattering_muon -> create_correction_functions();

        NuFit::interpolated_sys *nutau_scattering_muon = new interpolated_sys(par_name, muon, "NuTau", binsx, binsy, binsz);
        nutau_scattering_muon -> add_simulated_point(1.0, baselinedir_muon + std::string("nutau_muon.txt"), true); // this is baseline hist
        nutau_scattering_muon -> add_simulated_point(1.1, dir_muon + std::string("s+.10/nutau_muon.txt"));
        nutau_scattering_muon -> create_correction_functions();

        // map from flavor to correction
        std::unordered_map<std::string, NuFit::interpolated_sys *> muon_corrections;
        muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuE", nue_scattering_muon));
        muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuMu", numu_scattering_muon));
        muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuTau", nutau_scattering_muon));
        //muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("Muon", nue_scattering_muon));

        // store 
        binfits.insert(std::pair<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>>(muon, muon_corrections));
    }


    // hybrid event selection
    std::string hybrid("hybrid");
    if (std::find(analysis_names.begin(), analysis_names.end(), hybrid) != analysis_names.end())
    {
        std::vector<double> binsx = map_analyses.at(hybrid)->get_binsx();
        std::vector<double> binsy = map_analyses.at(hybrid)->get_binsy();
        std::vector<double> binsz = map_analyses.at(hybrid)->get_binsz();    

        NuFit::interpolated_sys *nue_scattering_hybrid = new interpolated_sys(par_name, hybrid, "NuE", binsx, binsy, binsz);
        nue_scattering_hybrid -> add_simulated_point(1.0, baselinedir+std::string("nue_hybrid.txt"), true); // this is baseline hist
        nue_scattering_hybrid -> add_simulated_point(1.1, dir+std::string("s+.10/nue_hybrid.txt"));
        nue_scattering_hybrid -> create_correction_functions();


        NuFit::interpolated_sys *numu_scattering_hybrid = new interpolated_sys(par_name, hybrid, "NuMu", binsx, binsy, binsz); 
        numu_scattering_hybrid -> add_simulated_point(1.0, baselinedir+std::string("numu_hybrid.txt"), true); // this is baseline hist
        numu_scattering_hybrid -> add_simulated_point(1.1, dir+std::string("s+.10/numu_hybrid.txt"));
        numu_scattering_hybrid -> create_correction_functions();

        NuFit::interpolated_sys *nutau_scattering_hybrid = new interpolated_sys(par_name, hybrid, "NuTau", binsx, binsy, binsz);
        nutau_scattering_hybrid -> add_simulated_point(1.0, baselinedir+std::string("nutau_hybrid.txt"), true); // this is baseline hist
        nutau_scattering_hybrid -> add_simulated_point(1.1, dir+std::string("s+.10/nutau_hybrid.txt"));
        nutau_scattering_hybrid -> create_correction_functions();

        // map from flavor to correction
        std::unordered_map<std::string, NuFit::interpolated_sys *> hybrid_corrections;
        hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuE", nue_scattering_hybrid));
        hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuMu", numu_scattering_hybrid));
        hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuTau", nutau_scattering_hybrid));
        //hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("Muon", nue_scattering_hybrid));

        // store
        binfits.insert(std::pair<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>>(hybrid, hybrid_corrections));
    }

    // specify what histograms are effected by this systematic
    components.clear();
    components.push_back(std::string("Conv"));
    components.push_back(std::string("Prompt"));
    components.push_back(std::string("Astro"));

    return;
}

NuFit::interpolated_par_scattering::~interpolated_par_scattering() { }


