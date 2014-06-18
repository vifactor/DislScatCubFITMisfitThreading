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
		throw ProgramSettings::Exception(toString(data["Q"].getPath()));
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
		    throw ProgramSettings::Exception(toString(
		                            dislocations["misfit"]["l"].getPath()));
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
		    threading_mixed.rho = 0.0;
		    threading_mixed.rc = 0.0;
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
		    threading_screw.rc = 0.0;
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
		    threading_mixed.rc = 0.0;
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

FitParameter readFParameter(const libconfig::Setting& stg)
{
	bool isFit;
	static FitParameter fparameter;

	if(stg.isScalar())
	{
		fparameter.m_Name = stg.getPath();
		fparameter.m_Value = stg;
		fparameter.m_Lbvalue = 0.0;
		fparameter.m_Ubvalue = 0.0;
	}else if(stg.isGroup())
	{
		isFit = stg["tofit"];
		if(isFit)
		{
			fparameter.m_Name = stg.getPath();
			fparameter.m_Value = stg["value"];
			fparameter.m_Lbvalue = stg["lbvalue"];
			fparameter.m_Ubvalue = stg["ubvalue"];

			if(!stg.lookupValue("scvalue", fparameter.m_Scvalue))
			{
				fparameter.m_Scvalue = 0.0;
			}

			if(fparameter.m_Lbvalue >= fparameter.m_Ubvalue)
			{
				throw ProgramSettings::Exception("Fault boundaries " + toString(stg.getPath()));
			}
		}
		else
		{
			fparameter.m_Name = stg.getPath();
			fparameter.m_Value = stg["value"];
			fparameter.m_Lbvalue = 0.0;
			fparameter.m_Ubvalue = 0.0;
		}
	}
	else
	{
		throw ProgramSettings::Exception("Inappropriate type " + toString(stg.getPath()));
	}
	return fparameter;
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

void ProgramSettings::read(const boost::filesystem::path& cfgdir)
{
	libconfig::Config cfg;
	libconfig::Config samplecfg, datacfg, fitcfg;
	boost::filesystem::path samplecfgfile, datacfgfile, fitcfgfile;
	
	m_cfgfile = cfgdir / "default.cfg";
	samplecfgfile = cfgdir / "sample.cfg";
	datacfgfile = cfgdir / "data.cfg";
	fitcfgfile = cfgdir / "fit.cfg";
	// Read the file. If there is an error, report it
	try
	{	
		/*FIXME temporary -----*/
		ProgramSettings::SampleConfig sample;
		samplecfg.readFile(samplecfgfile.c_str());
		samplecfg.setAutoConvert(true);
		const libconfig::Setting& sampleroot = samplecfg.getRoot();
		sample.set(sampleroot["Sample"], m_cpMap);
		std::cout << sample << std::endl;
		
		ProgramSettings::DataConfig data;
		datacfg.readFile(datacfgfile.c_str());
		datacfg.setAutoConvert(true);
		const libconfig::Setting& dataroot = datacfg.getRoot();
		data.set(dataroot["Data"], m_cpMap);
		std::cout << data << std::endl;
		
		ProgramSettings::FitConfig fit;
		fitcfg.readFile(fitcfgfile.c_str());
		fitcfg.setAutoConvert(true);
		const libconfig::Setting& fitroot = fitcfg.getRoot();
		fit.set(fitroot["Fit"], m_cpMap);
		std::cout << fit << std::endl;
		
        for (CalculatorParameterMap::iterator it=m_cpMap.begin(); it!=m_cpMap.end(); ++it)
            std::cout << it->first << " => " << it->second << '\n';
		/*FIXME temporary solution -----*/

        cfg.readFile(m_cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();

		readSampleSettings(root);
		readCalculatorSettings(root);
		readEngineSettings(root);
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

void ProgramSettings::readCalculatorSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &calculator = root["Calculator"];

	/*reflection*/
	if(calculator["Q"].isArray() && calculator["Q"].getLength() 
	                                               == MillerCubIndicesDimension)
	{
		m_calculatorSettings.Q.H = calculator["Q"][0];
		m_calculatorSettings.Q.K = calculator["Q"][1];
		m_calculatorSettings.Q.L = calculator["Q"][2];
	}
	else
	{
		throw ProgramSettings::Exception(toString(calculator["Q"].getPath()));
	}


	/*z-sampling for intensity calculations*/
	m_calculatorSettings.sampling = calculator["zsampling"];

	/*X-ray wavelength*/
	m_calculatorSettings.lambda = calculator["lambda"];

	m_calculatorSettings.scale = readFParameter(calculator["scale"]);
	m_calculatorSettings.background = readFParameter(calculator["background"]);

	m_calculatorSettings.qresolX = calculator["resolution"]["x"];
	m_calculatorSettings.qresolZ = calculator["resolution"]["z"];
}

void ProgramSettings::readSampleSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &sample = root["Sample"];

	/*lattice parameters*/
	m_sampleSettings.a0 = sample["a0"];

	/*Poisson ratio*/
	m_sampleSettings.nu = sample["nu"];

	/*Sample sizes*/
	m_sampleSettings.thickness = sample["thickness"];
	m_sampleSettings.width = sample["width"];

	/*dislocation settings*/
	const libconfig::Setting &dislocations = sample["dislocations"];

	/*potential fit parameters*/
	if(dislocations.exists("misfit"))
	{
		m_sampleSettings.misfit.rho = readFParameter(dislocations["misfit"]["rho"]);
		/*burgers vector*/
	    if(dislocations["misfit"]["b"].isArray() && 
	       dislocations["misfit"]["b"].getLength() == MillerCubIndicesDimension)
	    {
		    m_sampleSettings.misfit.b.X = dislocations["misfit"]["b"][0];
		    m_sampleSettings.misfit.b.Y = dislocations["misfit"]["b"][1];
		    m_sampleSettings.misfit.b.Z = dislocations["misfit"]["b"][2];
		}
		else
		    throw ProgramSettings::Exception(toString(
		                            dislocations["misfit"]["b"].getPath()));
		/*dislocation line*/
	    if(dislocations["misfit"]["l"].isArray() && 
	       dislocations["misfit"]["l"].getLength() 
	                                            == MillerCubIndicesDimension)
	    {
		    m_sampleSettings.misfit.l.X = dislocations["misfit"]["l"][0];
		    m_sampleSettings.misfit.l.Y = dislocations["misfit"]["l"][1];
		    m_sampleSettings.misfit.l.Z = dislocations["misfit"]["l"][2];
	    }
	    else
		    throw ProgramSettings::Exception(toString(
		                            dislocations["misfit"]["l"].getPath()));
	}
	else
	{
		m_sampleSettings.misfit.rho = FitParameter(dislocations.getPath() + ".misfit.rho");
	}

	if(dislocations.exists("threading"))
	{
		const libconfig::Setting &threading = dislocations["threading"];
		readThreadingDislocations(threading);
	}
	else
	{
		/*just a trick to initialize variable parameters*/
		readThreadingDislocations(dislocations);
	}
}

void ProgramSettings::readThreadingDislocations(const libconfig::Setting& stg)
{
	if(stg.exists("edge"))
	{
		m_sampleSettings.threading_edge.rho = readFParameter(stg["edge"]["rho"]);
		m_sampleSettings.threading_edge.rc = readFParameter(stg["edge"]["rc"]);
		m_sampleSettings.threading_edge.b_edge = m_sampleSettings.a0;
		m_sampleSettings.threading_edge.b_screw = 0.0;
	}
	else
	{
		m_sampleSettings.threading_mixed.rho = FitParameter(stg.getPath() + ".edge.rho");
		m_sampleSettings.threading_mixed.rc = FitParameter(stg.getPath() + ".edge.rc");
	}

	if(stg.exists("screw"))
	{
		m_sampleSettings.threading_screw.rho = readFParameter(stg["screw"]["rho"]);
		m_sampleSettings.threading_screw.rc = readFParameter(stg["screw"]["rc"]);
		m_sampleSettings.threading_screw.b_edge = 0;
		m_sampleSettings.threading_screw.b_screw = m_sampleSettings.a0;
	}
	else
	{
		m_sampleSettings.threading_screw.rho = FitParameter(stg.getPath() + ".screw.rho");
		m_sampleSettings.threading_screw.rc = FitParameter(stg.getPath() + ".screw.rc");
	}

	if(stg.exists("mixed"))
	{
		m_sampleSettings.threading_mixed.rho = readFParameter(stg["mixed"]["rho"]);
		m_sampleSettings.threading_mixed.rc = readFParameter(stg["mixed"]["rc"]);
		m_sampleSettings.threading_mixed.b_edge = m_sampleSettings.a0;
		m_sampleSettings.threading_mixed.b_screw = m_sampleSettings.a0;
	}
	else
	{
		m_sampleSettings.threading_mixed.rho = FitParameter(stg.getPath() + ".mixed.rho");
		m_sampleSettings.threading_mixed.rc = FitParameter(stg.getPath() + ".mixed.rc");
	}
}

void ProgramSettings::readEngineSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &engine = root["Engine"];

	/*output basename*/
	m_engineSettings.outfile = engine["outfile"].c_str();
	/*data file*/
	m_engineSettings.datafile = engine["datafile"].c_str();
	/*resume file*/
	m_engineSettings.resumefile = engine["resumefile"].c_str();
	/*max nb of iterations for fit*/
	m_engineSettings.nb_iter = engine["nbIterations"];
}

void ProgramSettings::print() const
{
	printEngineSettings();
	printSampleSettings();
	printCalculatorSettings();
}

void ProgramSettings::printCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;
	std::cout << "Reflection:\t" << m_calculatorSettings.Q << std::endl;
	std::cout << "X-ray wavelength:\t" << m_calculatorSettings.lambda << std::endl;
	std::cout << "Resolutions (dqx, dqz):\t" << m_calculatorSettings.qresolX
			<< "\t" << m_calculatorSettings.qresolZ << std::endl;
	std::cout << "Intensity scale coefficient:\t" << m_calculatorSettings.scale << std::endl;
	std::cout << "Intensity background:\t" << m_calculatorSettings.background << std::endl;
	std::cout << "Z sampling:\t" << m_calculatorSettings.sampling << std::endl;
}

void ProgramSettings::printSampleSettings() const
{
	std::cout << "---Sample settings---" << std::endl;
	std::cout << "Lattice parameters: (a0)\t" << m_sampleSettings.a0 << ", "
			<< std::endl;
	std::cout << "Sample sizes (thickness width):\t" << m_sampleSettings.thickness << "\t"
			<< m_sampleSettings.width << std::endl;
	std::cout << "Poisson ratio:\t" << m_sampleSettings.nu << std::endl;

	std::cout << "Misfit dislocations:" << std::endl;
	printMisfitDislocations();

	std::cout << "Threading dislocations:" << std::endl;
	printThreadingDislocations();

}

void ProgramSettings::printEngineSettings() const
{
	std::cout << "---Engine settings---" << std::endl;

	std::cout << "Data file:\t" << m_engineSettings.datafile << std::endl;
	std::cout << "Resume file:\t" << m_engineSettings.outfile << std::endl;
	std::cout << "Output basename:\t" << m_engineSettings.outfile << std::endl;
	std::cout << "Nb fit iterations:\t" << m_engineSettings.nb_iter << std::endl;
}

void ProgramSettings::printMisfitDislocations() const
{
	std::cout << "\tBurgers vector:\t" << m_sampleSettings.misfit.b << std::endl;
	std::cout << "\tDislocation line:\t" << m_sampleSettings.misfit.l << std::endl;
	std::cout << "\tDensity :\t" << m_sampleSettings.misfit.rho << std::endl;
}

void ProgramSettings::printThreadingDislocations() const
{

	std::cout << "\tEdge" << std::endl;
	std::cout << "\t\tBurgers vector:\t"
			<< m_sampleSettings.threading_edge.b_edge << std::endl;
	std::cout << "\t\tDensity:\t" << m_sampleSettings.threading_edge.rho
			<< std::endl;
	std::cout << "\t\tCorrelation radius:\t"
			<< m_sampleSettings.threading_edge.rc << std::endl;

	std::cout << "\tScrew" << std::endl;
	std::cout << "\t\tBurgers vector :\t"
			<< m_sampleSettings.threading_screw.b_screw << std::endl;
	std::cout << "\t\tDensity:\t" << m_sampleSettings.threading_screw.rho
			<< std::endl;
	std::cout << "\t\tCorrelation radius:\t"
			<< m_sampleSettings.threading_screw.rc << std::endl;

	std::cout << "\tMixed" << std::endl;
	std::cout << "\t\tBurgers vector (edge):\t"
			<< m_sampleSettings.threading_mixed.b_edge << std::endl;
	std::cout << "\t\tBurgers vector (screw):\t"
			<< m_sampleSettings.threading_mixed.b_screw << std::endl;
	std::cout << "\t\tDensity:\t" << m_sampleSettings.threading_mixed.rho
			<< std::endl;
	std::cout << "\t\tCorrelation radius:\t"
			<< m_sampleSettings.threading_mixed.rc << std::endl;
}
