#include "../../include/systematics/interpolated_sys.h"

NuFit::interpolated_sys::interpolated_sys(std::string parameter_name, std::string analysis_name, std::string nu_flavor, std::vector<double> bins_x, std::vector<double> bins_y, std::vector<double> bins_z)
{
    par_name = parameter_name;
    name = analysis_name;
    flavor = nu_flavor;

    binsx = bins_x;
    binsy = bins_y;
    binsz = bins_z;

    verbose = false;
}

void NuFit::interpolated_sys::add_simulated_point(double value, std::string infile, bool is_baseline)
{
    // if this is not the first point, check the this point is larger than the last one
    if(values.size()>0) 
    {
        if(value <= values[values.size()-1]) 
        {
            std::cout << "!!! FATAL: systematics datasets have to be added in strictly ascending order." << std::endl;
            std::cout << "received value: " << value << " but last value: " << values[values.size()-1] << std::endl;
            std::cout << "... exiting" << std::endl;
            exit(1);
        }
    }

    std::cout << std::endl;
    std::cout << "... reading " << par_name << " systematics for analysis " << name << " and flavor " << flavor << std::endl;
    std::cout << "... from file: " << infile << std::endl;
    std::cout << "... simulated " << par_name << " value: " << value << std::endl;

    std::string histname=par_name+std::string("_")+name+std::string("_")+flavor+std::string("_")+std::to_string(value);
    NuFit::neutrino_input simdata(histname, binsx, binsy, binsz);

    // fill systematics histograms
    simdata.read(infile);
    sim_data.push_back(simdata);
    values.push_back(value);

    // remember whether this discrete point is baseline
    baseline.push_back(is_baseline);

    std::cout << "... done" << std::endl;
    std::cout << std::endl;    
}

//void NuFit::interpolated_sys::create_correction_functions()
//{
//    return NuFit::interpolated_sys::create_correction_functions(false);
//}

void NuFit::interpolated_sys::create_correction_functions(bool use_interpolation)
{
    // loop over bins and add fit functions to matrix of fit functions
    // first fill histograms
    
    fill_hists();

    // find baseline hist:
    unsigned int baseline_index = 0;
    bool found = false;
    for (unsigned int i=0; i<baseline.size(); ++i)
    {
        if(baseline[i]) 
        {
            baseline_index = i;
            found = true;
            baseline_value = values[i];
            break;
        }
    }

    if (!found) 
    {
        // user did not add a baseline dataset
        std::cout << "!!! FATAL: no dataset is marked as baseline. User needs to add a baseline dataset." << std::endl;
        std::cout << "... exiting" << std::endl;
        exit(1);
    }

    // now create vector of fit functions
    NuFit::neutrino_input &dataset = sim_data.at(baseline_index);

    std::vector<std::vector<std::vector<TF1>>> vec;

    correction_functions.insert(std::pair<std::string, std::vector<std::vector<std::vector<TF1>>>> ("Conv", vec));
    correction_functions.insert(std::pair<std::string, std::vector<std::vector<std::vector<TF1>>>> ("Prompt", vec));
    correction_functions.insert(std::pair<std::string, std::vector<std::vector<std::vector<TF1>>>> ("astro", vec));

    correction_functions["Conv"].resize(dataset.get_nbinsx(), std::vector<std::vector<TF1>>(dataset.get_nbinsy(), std::vector<TF1>(dataset.get_nbinsz())));
    correction_functions["Prompt"].resize(dataset.get_nbinsx(), std::vector<std::vector<TF1>>(dataset.get_nbinsy(), std::vector<TF1>(dataset.get_nbinsz())));
    correction_functions["Astro"].resize(dataset.get_nbinsx(), std::vector<std::vector<TF1>>(dataset.get_nbinsy(), std::vector<TF1>(dataset.get_nbinsz())));

    std::vector<std::vector<std::vector<TFitResultPtr>>> vec2;
    correction_fit_results.insert(std::pair<std::string, std::vector<std::vector<std::vector<TFitResultPtr>>>> ("Conv", vec2));
    correction_fit_results.insert(std::pair<std::string, std::vector<std::vector<std::vector<TFitResultPtr>>>> ("Prompt", vec2));
    correction_fit_results.insert(std::pair<std::string, std::vector<std::vector<std::vector<TFitResultPtr>>>> ("astro", vec2));

    correction_fit_results["Conv"].resize(dataset.get_nbinsx(), std::vector<std::vector<TFitResultPtr>>(dataset.get_nbinsy(), std::vector<TFitResultPtr>(dataset.get_nbinsz())));
    correction_fit_results["Prompt"].resize(dataset.get_nbinsx(), std::vector<std::vector<TFitResultPtr>>(dataset.get_nbinsy(), std::vector<TFitResultPtr>(dataset.get_nbinsz())));
    correction_fit_results["Astro"].resize(dataset.get_nbinsx(), std::vector<std::vector<TFitResultPtr>>(dataset.get_nbinsy(), std::vector<TFitResultPtr>(dataset.get_nbinsz())));

    std::vector<std::vector<std::vector<double>>> vec3;
    correction_errors.insert(std::pair<std::string, std::vector<std::vector<std::vector<double>>>> ("Conv", vec3));
    correction_errors.insert(std::pair<std::string, std::vector<std::vector<std::vector<double>>>> ("Prompt", vec3));
    correction_errors.insert(std::pair<std::string, std::vector<std::vector<std::vector<double>>>> ("astro", vec3));

    correction_errors["Conv"].resize(dataset.get_nbinsx(), std::vector<std::vector<double>>(dataset.get_nbinsy(), std::vector<double>(dataset.get_nbinsz())));
    correction_errors["Prompt"].resize(dataset.get_nbinsx(), std::vector<std::vector<double>>(dataset.get_nbinsy(), std::vector<double>(dataset.get_nbinsz())));
    correction_errors["Astro"].resize(dataset.get_nbinsx(), std::vector<std::vector<double>>(dataset.get_nbinsy(), std::vector<double>(dataset.get_nbinsz())));

    success_conv.resize(dataset.get_nbinsx(), std::vector<std::vector<bool>>(dataset.get_nbinsy(), std::vector<bool>(dataset.get_nbinsz())));
    success_prompt.resize(dataset.get_nbinsx(), std::vector<std::vector<bool>>(dataset.get_nbinsy(), std::vector<bool>(dataset.get_nbinsz())));
    success_astro.resize(dataset.get_nbinsx(), std::vector<std::vector<bool>>(dataset.get_nbinsy(), std::vector<bool>(dataset.get_nbinsz())));

    unsigned int size = sim_data.size();
    std::vector<TH3D *> hists_conv;
    hists_conv.reserve(size);
    std::vector<TH3D *> hists_prompt;
    hists_prompt.reserve(size);
    std::vector<TH3D *> hists_astro;
    hists_astro.reserve(size);
        
    for (unsigned int k=0; k<dataset.get_nbinsx(); ++k) 
    {
        for (unsigned int l=0; l<dataset.get_nbinsy(); ++l) 
        {
            for (unsigned int m=0; m<dataset.get_nbinsz(); ++m) 
            {
                hists_conv.clear();
                hists_prompt.clear();
                hists_astro.clear();

                double baseline_expectation_conv = dataset.conv.GetBinContent(k+1, l+1, m+1);
                double baseline_error_conv = dataset.conv.GetBinError(k+1, l+1, m+1);

                //std::cout << name << " ratio conv (" << k << "," << l << "," << m << ")=" << sim_data[size-1].conv.GetBinContent(k+1,l+1,m+1) / baseline_expectation_conv << std::endl;

                double baseline_expectation_prompt = dataset.prompt.GetBinContent(k+1, l+1, m+1);
                double baseline_error_prompt = dataset.prompt.GetBinError(k+1, l+1, m+1);

                //std::cout << name << " ratio prompt (" << k << "," << l << "," << m << ")=" << sim_data[size-1].prompt.GetBinContent(k+1,l+1,m+1) / baseline_expectation_prompt << std::endl;

                double baseline_expectation_astro = dataset.astro.GetBinContent(k+1, l+1, m+1);
                double baseline_error_astro = dataset.astro.GetBinError(k+1, l+1, m+1);

                //std::cout << name << " ratio astro (" << k << "," << l << "," << m << ")=" << sim_data[size-1].astro.GetBinContent(k+1,l+1,m+1) / baseline_expectation_astro << std::endl;

                for(unsigned int n=0; n<size; ++n)
                {
                    hists_conv.push_back(&(sim_data[n].conv));
                    hists_prompt.push_back(&(sim_data[n].prompt));
                    hists_astro.push_back(&(sim_data[n].astro));
                }

                // create interpolations for conventional, prompt, conv
                create_fit(baseline_expectation_conv, baseline_error_conv, hists_conv, k, l ,m, success_conv, "Conv", use_interpolation);
                create_fit(baseline_expectation_prompt, baseline_error_prompt, hists_prompt, k, l ,m, success_prompt, "Prompt", use_interpolation);
                create_fit(baseline_expectation_astro, baseline_error_astro, hists_astro, k, l ,m, success_astro, "Astro", use_interpolation);
            }
        }
    }

    // clear all the input arrays. we only need to keep the correction functions
    for(unsigned int i=0; i<size; ++i) sim_data[i].clear();    
}

void NuFit::interpolated_sys::create_fit(const double &baseline_val, const double &baseline_err, std::vector<TH3D *> &hists, const unsigned int &k, const unsigned int &l, const unsigned int &m, std::vector<std::vector<std::vector<bool>>> &success, std::string component, bool use_interpolation)
{
    // subtract one, since we only fit the systematics datasets. 
    // the interpolation is forced to go through baseline point
    // thus the baseline point is not another degree of freedom

    std::vector<double> ratios; // y-values of fit  
    //ratios.reserve(hists.size()-1);
    std::vector<double> ratios_err;
    //ratios_err.reserve(hists.size()-1);
    std::vector<double> vals; // x-values of fit
    //vals.reserve(hists.size()-1);
    std::vector<double> vals_err;
    //vals_err.reserve(hists.size()-1);

    // setup TF1
    std::string fitname(name+"_fit_"+std::to_string(k)+std::to_string(l)+std::to_string(m));
    //std::string formula("[0]+[1]*(x-"+std::to_string(baseline_value)+")");
    std::string formula("1.0+[0]*(x-"+std::to_string(baseline_value)+")");
    TF1 fitfunc(fitname.c_str(), formula.c_str(), -10., 10.);


    bool valid_point = true;

    if(baseline_val==0)
    {
        valid_point = false;
        if( ((flavor.compare("NuTau") != 0) || ((flavor.compare("NuTau") == 0)&&(component.compare("Astro")==0))) && verbose )
        {
            std::cout << "!!! empty baseline expectation in histogram bin " << k+1 << " " << l+1 << " " << m+1 << " -> cant compute ratio." <<std::endl;
            std::cout << "... ignoring " << par_name << " systematics for flavor " << flavor << " and component "<< component << " in " << name << " for this bin." << std::endl;
            std::cout << std::endl;
        }

        // disable systematics parametrization for this bin
        //fitfunc.SetParameter(0, 1.0);
        //fitfunc.SetParameter(1, 0.0);
        fitfunc.SetParameter(0, 0.0);
        correction_functions[component][k][l][m]=fitfunc;        
    }
    else 
    {
        for (unsigned int n=0; n<hists.size(); ++n)
        {    
            if(baseline[n]==false)
            {
                double bincontent = hists[n]->GetBinContent(k+1,l+1,m+1);
                double binerror = hists[n]->GetBinError(k+1,l+1,m+1);
                double ratio = bincontent / baseline_val;
                ratios.push_back(ratio);
                ratios_err.push_back(ratio * TMath::Sqrt(TMath::Power(binerror/bincontent,2)+TMath::Power(baseline_err/baseline_val,2)));
                vals.push_back(values[n]);
                vals_err.push_back(0.0);
            
                // need to check whether result is reasonable
                // skip per-bin parametrization if uncertainty on ratio > 10%
                // or if result is nan

                /*
                if (ratios_err[n] / ratios[n] > 0.1)
                {
                    valid_point = false;
                    std::cout << "!!! insufficient statistics in histogram bin " << k+1 << " " << l+1 << " " << m+1 << " -> error on ratio > 10%. (current parameter point: " << values[n] << ")" << std::endl;
                    std::cout << "... ignoring " << par_name << " systematics for flavor " << flavor << " and component "<< component << " in " << name << " for this bin." << std::endl;
                    std::cout << std::endl;
                }
                */
            }
        }
    
        if (valid_point)
        {
            if (use_interpolation)
            {
                TF1 fitfunc;
                std::vector<std::vector<double>> x{{0.9,1,1.1},{0.9,1,1.1}};
                std::vector<std::vector<double>> y{{0.9,1,1.1},{1.8,2,2.2}};
                //char name_char[10+sizeof(char)];
                //std::sprintf(name_char,"f");
                ROOT::Math::Interpolator *itp = new ROOT::Math::Interpolator(3,ROOT::Math::Interpolation::kLINEAR);
                itp->SetData(vals,ratios);
                //itp123->SetData(x[0],y[0]);
                NuFit::root_interpolator_wrapper *inter_wrapper = new NuFit::root_interpolator_wrapper(itp);
                fitfunc = TF1("f",inter_wrapper,&NuFit::root_interpolator_wrapper::Eval,5,3000,0,"NuFit::root_interpolator_wrapper","Eval");

                correction_functions[component][k][l][m]=fitfunc;

                //std::cout<<"construction test Eval"<<correction_functions[component][k][l][m].Eval(500)<<std::endl;
            }
            else{
                // data to be interpolated/fitted
                TGraphErrors graph(ratios.size(), &(vals[0]), &(ratios[0]), &(vals_err[0]), &(ratios_err[0]));
                int status = graph.Fit(&fitfunc, "QS");        
                if ((int(status)==0) && ((fitfunc.GetProb() > 0.001) || (fitfunc.GetNDF()==0)))
                {
                    correction_functions[component][k][l][m]=fitfunc;
                    correction_fit_results[component][k][l][m]=status;
                }
                else 
                {
                    valid_point = false;
                }

                // some plotting functions 
                // need all points (also baseline)
                        /*
                std::vector<double> x_all;
                std::vector<double> y_all;
                for (unsigned int n=0; n<hists.size(); ++n)
                {
                    x_all.push_back(values[n]);
                    double bincontent = hists[n]->GetBinContent(k+1,l+1,m+1);
                    y_all.push_back(bincontent/baseline_val);
                }

                TGraph graph_all(x_all.size(), &(x_all[0]), &(y_all[0])); 
                            graph_all.SetMarkerColor(4);
                            graph_all.SetMarkerStyle(21);
                    
                double pval = fitfunc.GetProb();
                double par0 = fitfunc.GetParameter(0);
                //double par1 = fitfunc.GetParameter(1);
                
                std::string thisbin=name+std::string(" ")+flavor+std::string(" ")+component+std::string("):")+std::to_string(k)+std::string("-")+std::to_string(l)+std::string("-")+std::to_string(m);
                TCanvas c1("binfits", "binfits", 200, 10, 700, 500);
                c1.SetFillColor(0);
                c1.SetGrid();

                graph.SetLineColor(1);
                   graph.SetLineWidth(2);
                   graph.SetMarkerColor(2);
                   graph.SetMarkerStyle(21);
                   //graph.SetTitle((std::string("bin-fit(")+thisbin+std::string(", p=")+std::to_string(pval)+std::string(", par0=")+std::to_string(par0)+std::string(", par1=")+std::to_string(par1)).c_str());
                graph_all.SetTitle((std::string("bin-fit(")+thisbin+std::string(", p=")+std::to_string(pval)+std::string(", par0=")+std::to_string(par0)).c_str());
                   graph_all.GetXaxis()->SetTitle(par_name.c_str());
                   graph_all.GetYaxis()->SetTitle("ratio");
                graph_all.Draw("AP");
                   graph.Draw("P");        
        
                c1.SaveAs((std::string("./pdf/")+par_name+std::string("_")+name+std::string("_")+flavor+std::string("_")+component+std::string("_")+std::to_string(k)+std::string("_")+std::to_string(l)+std::string("_")+std::to_string(m)+std::string(".jpg")).c_str()); */
            }
        }
        
        if ( (!valid_point))
        {
            // invalid point    
            // disable systematics parametrization for this bin
            if (verbose) std::cout << "... ignoring " << par_name << " systematics for flavor " << flavor << " and component "<< component << " in " << name << " for bin " << k << " " << l << " " << m << std::endl;    
            //fitfunc.SetParameter(0, 1.0);
            //fitfunc.SetParameter(1, 0.0);
            fitfunc.SetParameter(0, 0.0);
            correction_functions[component][k][l][m]=fitfunc;
        }
    }
    success[k][l][m] = valid_point;
}

void NuFit::interpolated_sys::fill_hists()
{
    for (unsigned int i=0; i<sim_data.size(); ++i)
    {
        NuFit::neutrino_input &dataset = sim_data[i];
        //dataset.conv.Reset();
        dataset.conv.Sumw2();
        //dataset.prompt.Reset();
        dataset.prompt.Sumw2();
        //dataset.astro.Reset();
        dataset.astro.Sumw2();
        for (unsigned int j=0; j<dataset.get_size(); ++j) 
        {
            dataset.conv.Fill(dataset.logenergy_rec[j], dataset.coszenith_rec[j], dataset.ra_rec[j], dataset.conv_weight[j]);
            dataset.prompt.Fill(dataset.logenergy_rec[j], dataset.coszenith_rec[j], dataset.ra_rec[j], dataset.prompt_weight[j]);
            dataset.astro.Fill(dataset.logenergy_rec[j], dataset.coszenith_rec[j], dataset.ra_rec[j], dataset.astro_weight[j] * 2.2 * 1.e-18 * TMath::Power(dataset.energy_prim[j] / 1.e5, -2.5)); // assume astro for global fit
            //dataset.astro.Fill(dataset.logenergy_rec[j], dataset.coszenith_rec[j], dataset.ra_rec[j], dataset.astro_weight[j] * 2.2 * 1.e-18 * TMath::Power(dataset.energy_prim[j] / 1.e5, -2.54)); // assume astro for global fit
        }

    }

        /* diagnostics
                std::cout << name << std::endl;
                // check what's filled into histograms
                for (unsigned int i=0; i<sim_data[0].get_nbinsx(); ++i)
                {
                        std::cout << sim_data[0].astro.GetXaxis()->GetBinCenter(i+1) << " ";
                        for (unsigned int j=0; j<sim_data.size(); ++j)
                        {
                                std::cout << sim_data[j].astro.GetBinContent(i+1,1,1) << " ";
                        }
                        std::cout << "\n";
                }
                std::cout << std::endl;
        */

}

double NuFit::interpolated_sys::get_efficiency_correction(const double &x, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz) {
    // root histograms use range 1 to N to index bins
    // but correction factors are stored in array of 0 to N-1
    // if you are interested in correction for ROOT histogram bin k, l, m
    // this function needs to be queried with k-1, l-1, m-1    
    
    //return correction_functions[component][binx][biny][binz].Eval(x);

    //std::cout<<"componet"<<std::endl;
    //std::cout<<component<<std::endl;
    //std::cout<<"NuFit::interpolated_sys::get_efficiency_correction::k,l,m"<<std::endl;
    //std::cout<<binx<<" "<<biny<<" "<<binz<<std::endl;
    
    double corr = correction_functions[component][binx][biny][binz].Eval(x);

    if((corr>0.1)&&(corr<2.0))
        return corr;
    else if(corr<=0.1)
        return 0.1;
    else
        return 2.0;

}


double NuFit::interpolated_sys::get_relative_correction_error(const double &x, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz) {
    // root histograms use range 1 to N to index bins
    // but correction factors are stored in array of 0 to N-1
    // if you are interested in correction for ROOT histogram bin k, l, m
    // this function needs to be queried with k-1, l-1, m-1    
    
    double corr = correction_functions[component][binx][biny][binz].Eval(x);
    if((corr>0.1)&&(corr<2.0))
        corr = corr;
    else if(corr<=0.1)
        corr = 0.1;
    else
        corr = 2.0;

    double sigma[1];
    double xx[1];
    xx[0] = x;
    correction_fit_results[component][binx][biny][binz]->GetConfidenceIntervals(1,1,1,xx,sigma,0.68);

    double relative_correction_error = sigma[0]/corr;

    return relative_correction_error;


}
