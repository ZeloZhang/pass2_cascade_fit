#include "../include/hists.h"

NuFit::hists::hists(std::string name, std::vector<double> userbins_x, std::vector<double> userbins_y, std::vector<double> userbins_z) 
	: atm_conv((name+std::string("_all_conv")).c_str(), (name+std::string("_all_conv")).c_str(), userbins_x.size()-1, &(userbins_x[0]), userbins_y.size()-1, &(userbins_y[0]), userbins_z.size()-1, &(userbins_z[0])),
	  atm_prompt((name+std::string("_all_prompt")).c_str(), (name+std::string("_all_prompt")).c_str(), userbins_x.size()-1, &(userbins_x[0]), userbins_y.size()-1, &(userbins_y[0]), userbins_z.size()-1, &(userbins_z[0])), 
	  astro((name+std::string("_all_astro")).c_str(), (name+std::string("_all_astro")).c_str(), userbins_x.size()-1, &(userbins_x[0]), userbins_y.size()-1, &(userbins_y[0]), userbins_z.size()-1, &(userbins_z[0])),
	  mcsum((name+std::string("_mcsum")).c_str(), (name+std::string("_mcsum")).c_str(), userbins_x.size()-1, &(userbins_x[0]), userbins_y.size()-1, &(userbins_y[0]), userbins_z.size()-1, &(userbins_z[0])), 
	 

	  nue((name+std::string("_nue")), userbins_x, userbins_y, userbins_z),
	  numu((name+std::string("_numu")), userbins_x, userbins_y, userbins_z),
	  nutau((name+std::string("_nutau")), userbins_x, userbins_y, userbins_z),
          muon((name+std::string("_muon")), userbins_x, userbins_y, userbins_z), 
	  data((name+std::string("_data")), userbins_x, userbins_y, userbins_z)
           
{ 

	// keep a copy of input bins
	binsx = userbins_x;
	nbinsx = userbins_x.size()-1;

	binsy = userbins_y;
	nbinsy = userbins_y.size()-1;

	binsz = userbins_z;
	nbinsz = userbins_z.size()-1;

        (*this).name = name;

	/*
	nue_components.insert(std::pair<std::string, TH3D*>("Conv", &(nue.conv)));
	nue_components.insert(std::pair<std::string, TH3D*>("Prompt", &(nue.prompt)));
	nue_components.insert(std::pair<std::string, TH3D*>("Astro", &(nue.astro)));


        numu_components.insert(std::pair<std::string, TH3D*>("Conv", &(numu.conv)));
        numu_components.insert(std::pair<std::string, TH3D*>("Prompt", &(numu.prompt)));
        numu_components.insert(std::pair<std::string, TH3D*>("Astro", &(numu.astro)));


        nutau_components.insert(std::pair<std::string, TH3D*>("Conv", &(nutau.conv)));
        nutau_components.insert(std::pair<std::string, TH3D*>("Prompt", &(nutau.prompt)));
        nutau_components.insert(std::pair<std::string, TH3D*>("Astro", &(nutau.astro)));

	nuhists.insert(std::pair<std::string, std::map<std::string, TH3D *>>("NuE", nue_components));
	nuhists.insert(std::pair<std::string, std::map<std::string, TH3D *>>("NuMu", numu_components));
	nuhists.insert(std::pair<std::string, std::map<std::string, TH3D *>>("NuTau", nutau_components));
	*/

}

void NuFit::hists::read(std::string f_nue, std::string f_numu, std::string f_nutau, std::string f_muon, std::string f_data) {
	// read input from all neutrinos
	// and initialize data histogram
	
	std::cout << "... preparing event selection: " << name << std::endl;
	std::cout << std::endl;
	
	std::cout << "... reading nue" << std::endl;
	nue.read(f_nue);
	std::cout << "... done" << std::endl;
	std::cout << std::endl;
	std::cout << "... reading numu" << std::endl;
	numu.read(f_numu);
	std::cout << "... done" << std::endl;
	std::cout << std::endl;
	std::cout << "... reading nutau" << std::endl;
	nutau.read(f_nutau);
	std::cout << "... done" << std::endl;
	std::cout << std::endl;
	std::cout << "... reading muon" << std::endl;
	muon.read(f_muon);
	std::cout << "... done" << std::endl;
	std::cout << std::endl;
	std::cout << "... reading data" << std::endl;
	data.read(f_data);
	std::cout << "... done" << std::endl;
	std::cout << std::endl;


	

	// initialize data histograms
	std::cout << "... filling histograms" << std::endl;

	
	// nue	
	for (unsigned int i=0; i<nue.get_size(); ++i) {
		nue.conv.Fill(nue.logenergy_rec[i], nue.coszenith_rec[i], nue.ra_rec[i], nue.conv_weight[i]);
		nue.prompt.Fill(nue.logenergy_rec[i], nue.coszenith_rec[i], nue.ra_rec[i], nue.prompt_weight[i]);	
	}

	// numu
        for (unsigned int i=0; i<numu.get_size(); ++i) {
                numu.conv.Fill(numu.logenergy_rec[i], numu.coszenith_rec[i], numu.ra_rec[i], numu.conv_weight[i]);
                numu.prompt.Fill(numu.logenergy_rec[i], numu.coszenith_rec[i], numu.ra_rec[i], numu.prompt_weight[i]);
	}

	// nutau
        for (unsigned int i=0; i<nutau.get_size(); ++i) 
                nutau.prompt.Fill(nutau.logenergy_rec[i], nutau.coszenith_rec[i], nutau.ra_rec[i], nutau.prompt_weight[i]);
	



	// add nue, numu, nutau histograms
	
	atm_conv.Add(&nue.conv);	
	atm_conv.Add(&numu.conv);

	// copy the histogram
	atm_conv_orig = *((TH3D*) atm_conv.Clone());

	atm_prompt.Add(&nue.prompt);
	atm_prompt.Add(&numu.prompt);
	atm_prompt.Add(&nutau.prompt);

	// copy the histogram
	atm_prompt_orig = *((TH3D*) atm_prompt.Clone());

	//astro.Add(&nue.astro);
	//astro.Add(&numu.astro);
	//astro.Add(&nutau.astro);
	

	// muon
        for (unsigned int i=0; i<muon.get_size(); ++i)
		muon.hist.Fill(muon.logenergy_rec[i], muon.coszenith_rec[i], muon.ra_rec[i], muon.muon_weight[i]);

	// copy the histogram
	muon.hist_orig = *((TH3D*) muon.hist.Clone());
	

	// data
        for (unsigned int i=0; i<data.get_size(); ++i)
		data.hist.Fill(data.logenergy_rec[i], data.coszenith_rec[i], data.ra_rec[i]);


	std::cout << "... done" << std::endl;
	std::cout << std::endl;
}

unsigned int NuFit::hists::get_nbinsx() const {
	return nbinsx;
}

unsigned int NuFit::hists::get_nbinsy() const {
	return nbinsy;
}

unsigned int NuFit::hists::get_nbinsz() const {
	return nbinsz;
}

std::vector<double> NuFit::hists::get_binsx() {
	return binsx;
}

std::vector<double> NuFit::hists::get_binsy() {
	return binsy;
}

std::vector<double> NuFit::hists::get_binsz() {
	return binsz;
}

double NuFit::hists::get_bincontent_muon(const unsigned int &binx, const unsigned int &biny, const unsigned int &binz)
{
	// wrap around ROOTs TH3D --> bins range from i=1 to i=N
	return muon.hist.GetBinContent(binx, biny, binz);
}

void NuFit::hists::set_bincontent_muon(const unsigned int &binx, const unsigned int &biny, const unsigned int &binz, const double &value)
{
        // wrap around ROOTs TH3D --> bins range from i=1 to i=N
        muon.hist.SetBinContent(binx, biny, binz, value);
	return;
}

double NuFit::hists::get_bincontent(const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz)
{
	// wrap around ROOTs TH3D --> bins range from i=1 to i=N
		
	double bincontent = 0.0;
	if (!flavor.compare("NuE"))
	{
		//bincontent = nue_components[component]->GetBinContent(binx, biny, binz);
		if (!component.compare("Conv"))
			bincontent = nue.conv.GetBinContent(binx, biny, binz);
		else if (!component.compare("Prompt"))
			bincontent = nue.prompt.GetBinContent(binx, biny, binz);
		else if (!component.compare("Astro"))
        		bincontent = nue.astro.GetBinContent(binx, biny, binz);
	}
	else if (!flavor.compare("NuMu"))
	{
		//bincontent = numu_components[component]->GetBinContent(binx, biny, binz);
	        if (!component.compare("Conv"))
                        bincontent = numu.conv.GetBinContent(binx, biny, binz);
                else if (!component.compare("Prompt"))
                        bincontent = numu.prompt.GetBinContent(binx, biny, binz);
                else if (!component.compare("Astro"))
                        bincontent = numu.astro.GetBinContent(binx, biny, binz);	
	}
	else if (!flavor.compare("NuTau"))
	{
		//bincontent = nutau_components[component]->GetBinContent(binx, biny, binz);
                if (!component.compare("Conv"))
                        bincontent = nutau.conv.GetBinContent(binx, biny, binz);
                else if (!component.compare("Prompt"))
                        bincontent = nutau.prompt.GetBinContent(binx, biny, binz);
                else if (!component.compare("Astro"))
                        bincontent = nutau.astro.GetBinContent(binx, biny, binz);
	}
	return bincontent;
	
	//return nuhists[flavor][component]->GetBinContent(binx, biny, binz);
}

void NuFit::hists::set_bincontent(const std::string &flavor, const std::string &component, const unsigned int &binx, const unsigned int &biny, const unsigned int &binz, const double &value)
{
	
        if (!flavor.compare("NuE"))
        {
		//nue_components[component]->SetBinContent(binx, biny, binz, value);
                if (!component.compare("Conv"))
                        nue.conv.SetBinContent(binx, biny, binz, value);
                else if (!component.compare("Prompt"))
                        nue.prompt.SetBinContent(binx, biny, binz, value);
                else if (!component.compare("Astro"))
                        nue.astro.SetBinContent(binx, biny, binz, value);
        }
        else if (!flavor.compare("NuMu"))
        {
		//numu_components[component]->SetBinContent(binx, biny, binz, value);
                if (!component.compare("Conv"))
                        numu.conv.SetBinContent(binx, biny, binz, value);
                else if (!component.compare("Prompt"))
                        numu.prompt.SetBinContent(binx, biny, binz, value);
                else if (!component.compare("Astro"))
                         numu.astro.SetBinContent(binx, biny, binz, value);
        }
        else if (!flavor.compare("NuTau"))
        {
		//nutau_components[component]->SetBinContent(binx, biny, binz, value);
                if (!component.compare("Conv"))
                        nutau.conv.SetBinContent(binx, biny, binz, value);
                else if (!component.compare("Prompt"))
                        nutau.prompt.SetBinContent(binx, biny, binz, value);
                else if (!component.compare("Astro"))
                        nutau.astro.SetBinContent(binx, biny, binz, value);
        }	
	return;
	
	//nuhists[flavor][component]->SetBinContent(binx, biny, binz, value);
}
