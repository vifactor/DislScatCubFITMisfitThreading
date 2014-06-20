/*
 * ProgramSettings.cpp
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"
using namespace NonlinearFit;

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
    out << "--- Data settings--- " << std::endl;
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
		misfit.rho = dislocations["misfit"]["rho"];
        cmap[dislocations["misfit"]["rho"].getPath()] = misfit.rho;
		/*burgers vector*/
	    if(dislocations["misfit"]["b"].isArray() && 
	       dislocations["misfit"]["b"].getLength() == MillerCubIndicesDimension)
	    {
		    misfit.b.X = dislocations["misfit"]["b"][0];
		    misfit.b.Y = dislocations["misfit"]["b"][1];
		    misfit.b.Z = dislocations["misfit"]["b"][2];
		}
		else
		    throw ProgramSettings::Exception(toString(
		                            dislocations["misfit"]["b"].getPath()));
		/*dislocation line*/
	    if(dislocations["misfit"]["l"].isArray() && 
	       dislocations["misfit"]["l"].getLength() 
	                                            == MillerCubIndicesDimension)
	    {
		    misfit.l.X = dislocations["misfit"]["l"][0];
		    misfit.l.Y = dislocations["misfit"]["l"][1];
		    misfit.l.Z = dislocations["misfit"]["l"][2];
	    }
	    else
		    throw ProgramSettings::Exception(
		                            dislocations["misfit"]["l"].getPath());
	}
	else
	{
		misfit.rho = 0.0;
	}

	if(dislocations.exists("threading"))
	{
		const libconfig::Setting &threading = dislocations["threading"];
		if(threading.exists("edge"))
	    {
		    threading_edge.rho = threading["edge"]["rho"];
		    cmap[threading["edge"]["rho"].getPath()] = threading_edge.rho;
		    threading_edge.rc = threading["edge"]["rc"];
		    cmap[threading["edge"]["rc"].getPath()] = threading_edge.rc;
		    threading_edge.b_edge = a0;
		    threading_edge.b_screw = 0.0;
	    }
	    else
	    {
		    threading_edge.rho = 0.0;
		    cmap[threading.getPath() + ".edge.rho"] = threading_edge.rho;
		    threading_edge.rc = 0.0;
		    cmap[threading.getPath() + ".edge.rc"] = threading_edge.rc;
	    }

	    if(threading.exists("screw"))
	    {
		    threading_screw.rho = threading["screw"]["rho"];
		    cmap[threading["screw"]["rho"].getPath()] = threading_screw.rho;
		    threading_screw.rc = threading["screw"]["rc"];
		    cmap[threading["screw"]["rc"].getPath()] = threading_screw.rc;
		    threading_screw.b_edge = 0;
		    threading_screw.b_screw = a0;
	    }
	    else
	    {
		    threading_screw.rho = 0.0;
            cmap[threading.getPath() + ".screw.rho"] = threading_screw.rho;
		    threading_screw.rc = 0.0;
            cmap[threading.getPath() + ".screw.rc"] = threading_screw.rc;
	    }

	    if(threading.exists("mixed"))
	    {
		    threading_mixed.rho = threading["mixed"]["rho"];
		    cmap[threading["mixed"]["rho"].getPath()] = threading_mixed.rho;
		    threading_mixed.rc = threading["mixed"]["rc"];
		    cmap[threading["mixed"]["rho"].getPath()] = threading_mixed.rc;
		    threading_mixed.b_edge = a0;
		    threading_mixed.b_screw = a0;
	    }
	    else
	    {
		    threading_mixed.rho = 0.0;
		    cmap[threading.getPath() + ".mixed.rho"] = threading_mixed.rho;
		    threading_mixed.rc = 0.0;
		    cmap[threading.getPath() + ".mixed.rc"] = threading_mixed.rc;
	    }
	}
}

std::ostream&
operator<<(std::ostream& out, const ProgramSettings::SampleConfig &sample)
{
    out << "--- Sample settings ---" << std::endl;
	out << "\tLattice parameters: (a0)\t" << sample.a0 << ", "
			<< std::endl;
	out << "\tSample dimensions (thickness width):\t" 
            << sample.thickness << "\t"
			<< sample.width << std::endl;
	out << "\tPoisson ratio:\t" << sample.nu << std::endl;

	out << "\tMisfit dislocations:" << std::endl;
	out << "\t\tBurgers vector:\t" << sample.misfit.b << std::endl;
	out << "\t\tDislocation line:\t" << sample.misfit.l << std::endl;
	out << "\t\tDensity :\t" << sample.misfit.rho << std::endl;

	out << "\tThreading dislocations:" << std::endl;
	out << "\t\tEdge" << std::endl;
	out << "\t\t\tBurgers vector:\t"
			<< sample.threading_edge.b_edge << std::endl;
	out << "\t\t\tDensity:\t" << sample.threading_edge.rho
			<< std::endl;
	out << "\t\t\tCorrelation radius:\t"
			<< sample.threading_edge.rc << std::endl;

	out << "\t\tScrew" << std::endl;
	out << "\t\t\tBurgers vector :\t"	<< sample.threading_screw.b_screw 
                                        << std::endl;
	out << "\t\t\tDensity:\t" << sample.threading_screw.rho << std::endl;
	out << "\t\t\tCorrelation radius:\t"
			<< sample.threading_screw.rc << std::endl;

	out << "\t\tMixed" << std::endl;
	out << "\t\t\tBurgers vector (edge):\t"
			<< sample.threading_mixed.b_edge << std::endl;
	out << "\t\t\tBurgers vector (screw):\t"
			<< sample.threading_mixed.b_screw << std::endl;
	out << "\t\t\tDensity:\t" << sample.threading_mixed.rho << std::endl;
	out << "\t\t\tCorrelation radius:\t"
			<< sample.threading_mixed.rc << std::endl;
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

void ProgramSettings::read(const boost::filesystem::path& cfgdir)
{
	libconfig::Config samplecfg, datacfg, fitcfg;
	
	m_cfgfile = cfgdir / "default.cfg";
	m_samfile = cfgdir / "sample.cfg";
	m_datfile = cfgdir / "data.cfg";
	m_fitfile = cfgdir / "fit.cfg";
	// Read the file. If there is an error, report it
	try
	{	
		samplecfg.readFile(m_samfile.c_str());
		samplecfg.setAutoConvert(true);
		const libconfig::Setting& sampleroot = samplecfg.getRoot();
		m_sampleConfig.set(sampleroot["Sample"], m_cpMap);
		std::cout << m_sampleConfig << std::endl;
		
		datacfg.readFile(m_datfile.c_str());
		datacfg.setAutoConvert(true);
		const libconfig::Setting& dataroot = datacfg.getRoot();
		m_dataConfig.set(dataroot["Data"], m_cpMap);
		std::cout << m_dataConfig << std::endl;
		
		fitcfg.readFile(m_fitfile.c_str());
		fitcfg.setAutoConvert(true);
		const libconfig::Setting& fitroot = fitcfg.getRoot();
		m_fitConfig.set(fitroot["Fit"], m_cpMap);
		std::cout << m_fitConfig << std::endl;
		
        for (CalculatorParameterMap::const_iterator it=m_cpMap.begin();
                it!=m_cpMap.end(); ++it)
            std::cout << it->first << " => " << it->second << '\n';

	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(toString(fioex.what()) + " in\t" + m_cfgfile.native());
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				toString(pex.what()) + " in\t" + m_cfgfile.native() + ":"
						+ toString(pex.getLine()) + " - "
						+ toString(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				toString(nfex.what()) + "\t" + toString(nfex.getPath())
						+ " in\t" + m_cfgfile.native());
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				toString(tex.what()) + "\t" + toString(tex.getPath()) + " in\t"
						+ m_cfgfile.native());
	}
}
