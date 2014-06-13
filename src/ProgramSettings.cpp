/*
 * ProgramSettings.cpp
 *
 *  Created on: 11 черв. 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"
using namespace NonlinearFit;

Range readRange(const libconfig::Setting& stg)
{
	Range range;

	range.m_min = stg[0][0];
	range.m_max = stg[0][1];
	range.m_sampling = stg[1];

	return range;
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

std::ostream& operator<<(std::ostream& out, const Range& range)
{
	out << "[" << range.m_min << ", " << range.m_max << "]:" << range.m_sampling;
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
}

ProgramSettings::~ProgramSettings()
{
}

void ProgramSettings::read(const std::string& cfgfile)
{
	libconfig::Config cfg;

	m_cfgfile = cfgfile;
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(m_cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();

		readSampleSettings(root);
		readCalculatorSettings(root);
		readEngineSettings(root);
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(toString(fioex.what()) + " in\t" + cfgfile);
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				toString(pex.what()) + " in\t" + cfgfile + ":"
						+ toString(pex.getLine()) + " - "
						+ toString(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				toString(nfex.what()) + "\t" + toString(nfex.getPath())
						+ " in\t" + cfgfile);
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				toString(tex.what()) + "\t" + toString(tex.getPath()) + " in\t"
						+ cfgfile);
	}
}

void ProgramSettings::readCalculatorSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &calculator = root["Calculator"];

	/*reflection*/
	if(calculator["Q"].isArray() && calculator["Q"].getLength() == CalculatorSettings::HEXDIM)
	{
		for(int i = 0; i < CalculatorSettings::HEXDIM; ++i)
		{
			m_calculatorSettings.Q[i] = calculator["Q"][i];
		}
	}
	else
	{
		throw ProgramSettings::Exception(toString(calculator["Q"].getPath()));
	}
	/*check the property of hexagonal Miller indices*/
	if((m_calculatorSettings.Q[0] + m_calculatorSettings.Q[1] + m_calculatorSettings.Q[2]) != 0)
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
	m_sampleSettings.c0 = sample["c0"];

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
		m_sampleSettings.misfit.b_x = dislocations["misfit"]["b_x"];
		m_sampleSettings.misfit.b_z = dislocations["misfit"]["b_z"];
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
		m_sampleSettings.threading_screw.b_screw = m_sampleSettings.c0;
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
		m_sampleSettings.threading_mixed.b_screw = m_sampleSettings.c0;
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
	std::cout << "Reflection:\t[" << m_calculatorSettings.Q[0] << ", "
			<< m_calculatorSettings.Q[1] << ", " << m_calculatorSettings.Q[2]
			<< ", " << m_calculatorSettings.Q[3] << "]" << std::endl;
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
	std::cout << "Lattice parameters: (a0, c0)\t" << m_sampleSettings.a0 << ", "
			<< m_sampleSettings.c0 << std::endl;
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
	std::cout << "\tBurgers vector (bx):\t" << m_sampleSettings.misfit.b_x
			<< std::endl;
	std::cout << "\tBurgers vector (bz):\t" << m_sampleSettings.misfit.b_z
			<< std::endl;
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
