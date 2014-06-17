/*
 * ProgramSettings.cpp
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"
using namespace NonlinearFit;

/*----- Data settings -----*/
void ProgramSettings::DataConfig::set(const libconfig::Setting& data)
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
	Ibg = data["Ibg"];
	
    file = data["file"].c_str();
}

std::ostream& operator<<(std::ostream& out, const ProgramSettings::DataConfig &data)
{
    out << "---Data settings---" << std::endl;
	out << "\tReflection:\t" << data.Q << std::endl;
	out << "\tResolutions (dqx, dqz):\t" << data.resolX
			<< "\t" << data.resolZ << std::endl;
	out << "\tScale:\t" << data.I0 << std::endl;
	out << "\tBackground:\t" << data.Ibg << std::endl;
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

ProgramSettings::ProgramSettings()
{
}

ProgramSettings::~ProgramSettings()
{
}

void ProgramSettings::read(const boost::filesystem::path& cfgdir)
{
	libconfig::Config cfg;
	libconfig::Config datacfg;
	boost::filesystem::path datacfgfile;
	
	m_cfgfile = cfgdir / "default.cfg";
	datacfgfile = cfgdir / "data.cfg";
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(m_cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();
		
		/*FIXME temporary*/
		ProgramSettings::DataConfig data;
		datacfg.readFile(datacfgfile.c_str());
		datacfg.setAutoConvert(true);
		const libconfig::Setting& dataroot = datacfg.getRoot();
		data.set(dataroot["Data"]);
		std::cout << data << std::endl;

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
