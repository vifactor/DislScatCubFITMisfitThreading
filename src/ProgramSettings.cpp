/*
 * ProgramSettings.cpp
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"
using namespace NonlinearFit;
using namespace boost;

/*----- Data settings -----*/
void 
ProgramSettings::DataConfig::set(const libconfig::Setting& data,
                                CalculatorParameterMap& cmap) 
{
    if(data["Q"].isArray() && data["Q"].getLength() 
	                                               == MillerCubIndicesDimension)
	{
		Q.H = data["Q"][0];
		Q.K = data["Q"][1];
		Q.L = data["Q"][2];
	}
	else
	{
		throw ProgramSettings::Exception(data["Q"].getPath());
	}
    
    resolX = data["resolution"]["x"];
	resolZ = data["resolution"]["z"];
	
	I0 = data["I0"];
	cmap[data["I0"].getPath()] = I0;
	Ibg = data["Ibg"];
	cmap[data["Ibg"].getPath()] = Ibg;
	
    file = data["file"].c_str();
}

std::ostream&
operator<<(std::ostream& out, const ProgramSettings::DataConfig &data)
{
    out << "\tFile:\t" << data.file << std::endl;
	out << "\tReflection:\t" << data.Q << std::endl;
	out << "\tResolutions (dqx, dqz):\t" << data.resolX
			<< "\t" << data.resolZ << std::endl;
	out << "\tScale:\t" << data.I0 << std::endl;
	out << "\tBackground:\t" << data.Ibg << std::endl;
	return out;
}

/*----- Fit settings -----*/
void ProgramSettings::FitConfig::set(const libconfig::Setting& fit,
                const CalculatorParameterMap& cmap)
{
    FitParameter param;
    CalculatorParameterMap::const_iterator it;
    
	nbIterations = fit["nbIterations"];
    if(fit["parameters"].isList())
    {
        for(size_t i = 0; i < fit["parameters"].getLength(); ++i)
        {
            const libconfig::Setting& stg = fit["parameters"][i];
            param.m_Name = stg["name"].c_str();
            it = cmap.find(param.m_Name);
		    if(it != cmap.end())
		    {
		        param.m_Value = it->second;
		        param.m_Lbvalue = stg["minVal"];
		        param.m_Ubvalue = stg["maxVal"];
		        fitParameters.push_back(param);
		    }
	        else
	            std::cout << "No such cmap parameter:\t" << param.m_Name 
	                       << std::endl;
		}
    }
}

std::ostream&
operator<<(std::ostream& out, const ProgramSettings::FitConfig &fit)
{
    out << "--- Fit settings ---" << std::endl;
	out << "\tNb iterations:\t" << fit.nbIterations << std::endl;
	out << "\tNb fit parameters:\t" << fit.fitParameters.size() << std::endl;
    for(size_t i = 0; i < fit.fitParameters.size(); ++i)
    {
        const FitParameter& param = fit.fitParameters[i];
        out << "\t\t" << param.m_Name << ":\t" << param.m_Value <<
                " [" << param.m_Lbvalue << " : " << param.m_Ubvalue << "]"
                                                                 << std::endl;
    }
	return out;
}

/*----- Sample settings -----*/
void 
ProgramSettings::SampleConfig::set(const libconfig::Setting& sample,
                                 CalculatorParameterMap& cmap)
{
	/*lattice parameters*/
	a0 = sample["a0"];
	/*Poisson ratio*/
	nu = sample["nu"];
	/*Sample dimensions*/
	thickness = sample["thickness"];
	width = sample["width"];
	/*dislocation settings*/
	const libconfig::Setting &dislocations = sample["dislocations"];
	if(dislocations.exists("misfit"))
	{
	    /*find how many misfit dislocation types are present*/
	    int nb_mf_types = dislocations["misfit"].getLength();
	    misfit.resize(nb_mf_types);
	    for(int i = 0; i < nb_mf_types; ++i)
	    {
	        const libconfig::Setting & interface = 
	                    dislocations["misfit"][i];

            misfit[i].rho = interface["rho"];
            cmap[interface["rho"].getPath()] = misfit[i].rho;
		    /*burgers vector*/
	        if(interface["b"].isArray() && interface["b"].getLength()
	                                             == MillerCubIndicesDimension)
	        {
		        misfit[i].b.X = interface["b"][0];
		        misfit[i].b.Y = interface["b"][1];
		        misfit[i].b.Z = interface["b"][2];
		    }
		    else
		        throw ProgramSettings::Exception(
		                            interface["b"].getPath());
		    /*dislocation line*/
	        if(interface["l"].isArray() && interface["l"].getLength() 
	                                            == MillerCubIndicesDimension)
            {
		        misfit[i].l.X = interface["l"][0];
		        misfit[i].l.Y = interface["l"][1];
		        misfit[i].l.Z = interface["l"][2];
	        }
	        else
		        throw ProgramSettings::Exception(
		                            interface["l"].getPath());
	   }
	}

	if(dislocations.exists("threading"))
	{
	    /*find how many threading dislocation types are present*/
	    int nb_th_types = dislocations["threading"].getLength();
		threading.resize(nb_th_types);
	    for(int i = 0; i < nb_th_types; ++i)
	    {
            const libconfig::Setting & layer = 
                dislocations["threading"][i];

            threading[i].rho = layer["rho"];
            cmap[layer["rho"].getPath()] = threading[i].rho;
            threading[i].rc = layer["rc"];
            cmap[layer["rc"].getPath()] = threading[i].rc;
		    /*burgers vector*/
	        if(layer["b"].isArray() && layer["b"].getLength()
	                                             == MillerCubIndicesDimension)
	        {
		        threading[i].b.X = layer["b"][0];
		        threading[i].b.Y = layer["b"][1];
		        threading[i].b.Z = layer["b"][2];
		    }
		    else
		        throw ProgramSettings::Exception(
		                            layer["b"].getPath());
	    }
	}
}

std::ostream&
operator<<(std::ostream& out, const ProgramSettings::SampleConfig &sample)
{
    out << "--- Sample settings ---" << std::endl;
	out << "\tLattice parameters: (a0)\t" << sample.a0 << std::endl;
	out << "\tSample dimensions (thickness width):\t" 
            << sample.thickness << "\t"
			<< sample.width << std::endl;
	out << "\tPoisson ratio:\t" << sample.nu << std::endl;

    out << "\tMisfit dislocations:" << std::endl;
    for(size_t i = 0; i < sample.misfit.size(); ++i)
    {
        out << "\t\t[" << i << "]---" << std::endl;
        out << "\t\tBurgers vector:\t" << sample.misfit[i].b << std::endl;
        out << "\t\tDislocation line:\t" << sample.misfit[i].l << std::endl;
        out << "\t\tDensity :\t" << sample.misfit[i].rho << std::endl;
    }

	out << "\tThreading dislocations:" << std::endl;
	for(size_t i = 0; i < sample.threading.size(); ++i)
	{
	    out << "\t\t[" << i << "]---" << std::endl;
	    out << "\t\tBurgers vector:\t" << sample.threading[i].b << std::endl;
	    out << "\t\tDensity:\t" << sample.threading[i].rho << std::endl;
	    out << "\t\tCorrelation radius:\t" << sample.threading[i].rc << std::endl;
    }

	return out;
}

std::ostream& operator<<(std::ostream& out, const FitParameter& fparam)
{
	/*variable parameter*/
	if((fparam.m_Lbvalue == 0.0) && (fparam.m_Ubvalue == 0.0))
	{
		out << fparam.m_Value;
	}
	else /*fixed parameter*/
	{
		out << fparam.m_Value << " [" << fparam.m_Lbvalue << " - " << fparam.m_Ubvalue << "]";
	}
	return out;
}

ProgramSettings::ProgramSettings()
{
    m_resfile = "resume.txt";
}

void ProgramSettings::read(const filesystem::path& cfgdir)
{
	libconfig::Config samplecfg, datacfg, fitcfg;
	filesystem::path filename; 
	
	m_samfile = cfgdir / "sample.cfg";
	m_datfile = cfgdir / "data.cfg";
	m_fitfile = cfgdir / "fit.cfg";
	// Read the file. If there is an error, report it
	try
	{	
	    filename = m_samfile;
		samplecfg.readFile(filename.c_str());
		samplecfg.setAutoConvert(true);
		const libconfig::Setting& sampleroot = samplecfg.getRoot();
		m_sampleConfig.set(sampleroot["Sample"], m_cpMap);
		std::cout << m_sampleConfig << std::endl;
		
		filename = m_datfile;
		datacfg.readFile(filename.c_str());
		datacfg.setAutoConvert(true);
		const libconfig::Setting& dataroot = datacfg.getRoot();
		m_dataConfig.resize(dataroot["Data"].getLength());
		std::cout << "--- Data settings--- " << std::endl;
		for(size_t i = 0; i < m_dataConfig.size(); i++)
		{
		    m_dataConfig[i].set(dataroot["Data"][i], m_cpMap);
		    std::cout << m_dataConfig[i] << std::endl;
		}
		
		filename = m_fitfile;
		fitcfg.readFile(filename.c_str());
		fitcfg.setAutoConvert(true);
		const libconfig::Setting& fitroot = fitcfg.getRoot();
		m_fitConfig.set(fitroot["Fit"], m_cpMap);
		std::cout << m_fitConfig << std::endl;
		
		std::cout << "--- Calculator parameters ---" << std::endl;
        for (CalculatorParameterMap::const_iterator it=m_cpMap.begin();
                it!=m_cpMap.end(); ++it)
            std::cout << "\t" << it->first << " => " << it->second << std::endl;
        std::cout << std::endl;

	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(lexical_cast<std::string>(fioex.what()) 
		                + " in\t" + filename.native());
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				lexical_cast<std::string>(pex.what()) + " in\t" 
				+ filename.native() + ":"
						+ lexical_cast<std::string>(pex.getLine()) + " - "
						+ lexical_cast<std::string>(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				lexical_cast<std::string>(nfex.what()) + "\t" 
				+ lexical_cast<std::string>(nfex.getPath())
						+ " in\t" + filename.native());
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				lexical_cast<std::string>(tex.what()) + "\t"
				+ lexical_cast<std::string>(tex.getPath()) + " in\t"
						+ filename.native());
	}
}
