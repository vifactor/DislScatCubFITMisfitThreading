/*
 *
 *  Created on: 30 may 2012
 *      Author: kopp
 */

#include "Engine.h"

using namespace NonlinearFit;
using namespace Geometry;

/*----------- Handy functions --------------*/
Vector3d toMisfitFrame(const Vector3d& burgers,
                     const Vector3d& line,
                     const Vector3d& normal)
{
    Vector3d misf_dir;
    Vector3d line_dir;
    Vector3d norm_dir;
    
    line_dir = normalize(line);
    norm_dir = normalize(normal);
    misf_dir = cross(norm_dir, line_dir);
    
    /*misfit, screw and bz-edge components of the dislocation*/
    double bx = inner_prod(burgers, misf_dir);
    double by = inner_prod(burgers, line_dir);
    double bz = inner_prod(burgers, norm_dir);
    
    return Vector3d(bx, by, bz);
}

Vector3d toCoplanarFrame(const Vector3d& Q, const Vector3d& normal)
{
    Vector3d norm_dir;
    double Qx, Qy, Qz;

    norm_dir = normalize(normal);
    
    Qz = inner_prod(Q, norm_dir);
    Qy = 0.0;//in coplanar geometry Qy = 0
    Qx = norm_2(Q - Qz * norm_dir);
    
    return Vector3d(Qx, Qy, Qz);
}

double toMisfitInPlaneAngle(const Vector3d& Q, const Vector3d& burgers,
                const Vector3d& line, const Vector3d& normal)
{
    /* gives angle between in-plane component of reflection vector
     * and misfit component of Burgers vector
     */
    Vector3d misf_dir;
    Vector3d line_dir;
    Vector3d norm_dir;
    Vector3d Qpar_dir;
    Vector3d Qpar;
    double Qz;
    
    line_dir = normalize(line);
    norm_dir = normalize(normal);
    misf_dir = cross(norm_dir, line_dir);
    
    Qz = inner_prod(Q, norm_dir);
    Qpar = Q - Qz * norm_dir;
    Qpar_dir = normalize(Qpar);
    
    return acos(inner_prod(Qpar_dir, misf_dir));
}
/*----------- end handy functions -----------*/
Engine::Engine()
{
	m_calculator = NULL;
	m_fitter = NULL;
	m_programSettings = NULL;
	m_fit_calculator = NULL;
	m_sample = NULL;
}

Engine::~Engine()
{
	if(m_programSettings)
		delete m_programSettings;
	if(m_calculator)
		delete m_calculator;
	if(m_fit_calculator)
		delete m_fit_calculator;
	if(m_fitter)
		delete m_fitter;
	if(m_sample)
		delete m_sample;
	for(size_t i = 0; i < m_DataPoints.size(); ++i)
	{
		if(m_DataPoints[i].m_Argument)
			delete m_DataPoints[i].m_Argument;
	}
}

void Engine::exec(const boost::filesystem::path & workDir)
{
   
    m_WorkDir = workDir;
	m_programSettings = new ProgramSettings;
	m_programSettings->read(m_WorkDir);

	m_programSettings->print();
	init();
	setupComponents();
	doWork();
	saveResult();
}

void Engine::init()
{
	if(m_calculator)
		delete m_calculator;
	if(m_fit_calculator)
		delete m_fit_calculator;
	if(m_fitter)
		delete m_fitter;
	/*vector of fit parameters with their initial values*/
	m_fParametersInit.clear();
	/*vector of fit parameters with their final values*/
	m_fParametersFinal.clear();
	/*map of calculator parameters (potential fit parameters)*/
	m_cParameters.clear();

	m_exp_intens_vals.clear(); m_ini_intens_vals.clear(); m_fin_intens_vals.clear();
	m_qx_vals.clear(); m_qz_vals.clear();

	for(size_t i = 0; i < m_DataPoints.size(); ++i)
	{
		delete m_DataPoints[i].m_Argument;
	}
	m_DataPoints.clear();
}

void Engine::doWork()
{
	calcI(m_ini_intens_vals);
	fitI();
	printFitterInfo(m_fitter);
	calcI(m_fin_intens_vals);
	//saveSettings();
}

void Engine::calcI(std::vector<double>& vals)
{
	double I;

	m_fit_calculator->reinit(m_cParameters);

	for(size_t ipt = 0; ipt < m_DataPoints.size(); ipt++)
	{
		I = m_fit_calculator->eval(m_DataPoints[ipt].m_Argument);
		vals[ipt] = I;

		//std::cout << ipt << "\t" << I << std::endl;
	}
}

std::string Engine::getFilename() const
{
	std::string filename;
	std::string extension;

	filename = m_programSettings->getEngineSettings().outfile;
	extension = stripExtension(filename);
	filename += "_ft";
	filename += "." + extension;
	return filename;
}

void Engine::saveResult() const
{
	std::ofstream fout(getFilename().c_str());
	if(!fout)
	{
		throw Engine::Exception("Unable to open the output file:\t" + getFilename());
	}

	fout<<"#qx\tqz\tIexp\tIini\tIfin"<<std::endl;
	for(size_t ipt = 0; ipt < m_DataPoints.size(); ipt++)
	{
		fout<<m_qx_vals[ipt]<<"\t"
			<<m_qz_vals[ipt]<<"\t"
			<<m_exp_intens_vals[ipt]<<"\t"
			<<m_ini_intens_vals[ipt]<<"\t"
			<<m_fin_intens_vals[ipt]<<"\t"
				<<std::endl;
	}
	fout.close();
	saveResume();
}

void Engine::fitI()
{
	//carry out fit
	//with max "calculationSettings.modefitSettings.nbIterations" iterations
	m_fitter->fit(NonlinearFitter::fitLIN, m_programSettings->getEngineSettings().nb_iter);

	//get new values (best fit values) of fit parameters
	m_fParametersFinal = m_fitter->getFitParameters();

	//update calc parameters with new fitted values
	mergeParameters(m_cParameters, m_fParametersFinal);
}

void Engine::saveResume() const
{
	double rho_screw, rho_edge, rho_mixed,
		rc_screw, rc_edge, rc_mixed,
		M_edge, M_screw, M_mixed, rho_misfit;
	std::string filename;
	std::ofstream fout, fin;

	rc_edge = m_cParameters.find("Sample.dislocations.threading.edge.rc")->second;
	rc_screw = m_cParameters.find("Sample.dislocations.threading.screw.rc")->second;
	rc_mixed = m_cParameters.find("Sample.dislocations.threading.mixed.rc")->second;
	rho_edge = m_cParameters.find("Sample.dislocations.threading.edge.rho")->second;
	rho_screw = m_cParameters.find("Sample.dislocations.threading.screw.rho")->second;
	rho_mixed = m_cParameters.find("Sample.dislocations.threading.mixed.rho")->second;
	rho_misfit = m_cParameters.find("Sample.dislocations.misfit.rho")->second;

	/*transform rc to M = rc/rd = rc*rho^(1/2) */
	M_screw = rc_screw * sqrt(rho_screw * 1e-14);
	M_edge = rc_edge * sqrt(rho_edge * 1e-14);
	M_mixed = rc_mixed * sqrt(rho_mixed * 1e-14);

	filename = m_programSettings->getEngineSettings().resumefile;

	fin.open(filename.c_str(), std::ios::in);
	if(fin.good())
	{
		fin.close();
		fout.open(filename.c_str(), std::ios::app);
	}
	else
	{
		fin.close();
		fout.open(filename.c_str(), std::ios::out);
		fout << "#cfg\tdat\tH\tK\tL\t"<<
				"rho_edge\trc_edge\tM_edge\t" <<
				"rho_screw\trc_screw\tM_screw\t" <<
				"rho_mixed\trc_mixed\tM_mixed\trho_misfit\n";
	}
	fout << m_programSettings->getConfigfile() << "\t";
	fout << m_programSettings->getEngineSettings().datafile << "\t";
	fout << m_programSettings->getCalculatorSettings().Q.H << "\t" <<
			m_programSettings->getCalculatorSettings().Q.K << "\t" <<
			m_programSettings->getCalculatorSettings().Q.L << "\t";
			
	fout << rho_edge << "\t";
	fout << rc_edge << "\t";
	fout << M_edge << "\t";
	fout << rho_screw << "\t";
	fout << rc_screw << "\t";
	fout << M_screw << "\t";
	fout << rho_mixed << "\t";
	fout << rc_mixed << "\t";
	fout << M_mixed << "\t";
	fout << rho_misfit << "\t";
	fout << std::endl;
	fout.close();
}

void Engine::saveSettings() const
{
	libconfig::Config cfg;
	boost::filesystem::path oldcfgfile, newcfgfile;

	oldcfgfile = m_programSettings->getConfigfile();
	oldcfgfile.replace_extension(".~cfg");
	newcfgfile = m_programSettings->getConfigfile();
	/*save old configuration in a file with extension ~cfg*/
	boost::filesystem::rename(newcfgfile, oldcfgfile);

	try
	{
		// Read the configuration file. If there is an error, report it
		cfg.readFile(oldcfgfile.c_str());

		//reset configuration putting the new fitted values from fitParameterListFinal
		//for(auto fparameter : fitParameterListFinal)
		for(size_t i = 0; i < m_fParametersFinal.size(); ++i)
		{
			cfg.lookup(m_fParametersFinal[i].m_Name + ".value") = m_fParametersFinal[i].m_Value;
		}
		//save new configuration
		cfg.writeFile (newcfgfile.c_str());
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Engine::Exception("I/O error while reading file:\t" + oldcfgfile.native());
	} catch (const libconfig::ParseException &pex)
	{
		throw Engine::Exception("Parse error at " +
									toString(pex.getFile()) + ":" +
									toString(pex.getLine()) + " - " +
									toString(pex.getError()));
	}catch(const libconfig::SettingNotFoundException &nfex)
	{
		throw Engine::Exception(toString(nfex.getPath()));
	}catch(libconfig::SettingTypeException& tex)
	{
		throw Engine::Exception(toString(tex.getPath()) + "(" + toString(tex.what()) + ")");
	}
}

void Engine::setupComponents()
{
	setupCalculator();
	setupCParameters();
	readData();
	setupFitter();
}

void Engine::readData()
{
	DataReader dr;
	boost::filesystem::path datapath;
	
	datapath = m_WorkDir / m_programSettings->getEngineSettings().datafile;
	dr.readFile(datapath.c_str());

	if(!dr.good())
	{
		throw Engine::Exception("File " + datapath.native() +
				" has not been read or has no header");
	}

	if(dr.columnExist("[intensity]"))
	{
		dr.getColumn(m_exp_intens_vals, "[intensity]");
	}
	else
	{
		throw Engine::Exception("Column \"[intensity]\" has not been found in " +
				m_programSettings->getEngineSettings().datafile);
	}

	if (dr.columnExist("[qx]"))
	{
		//get points without any transformation
		dr.getColumn(m_qx_vals, "[qx]");
	}
	else
	{
		throw Engine::Exception("Column \"[qx]\" has not been found in "+
				m_programSettings->getEngineSettings().datafile);
	}

	if (dr.columnExist("[qz]"))
	{
		//get points without any transformation
		dr.getColumn(m_qz_vals, "[qz]");
	}
	else
	{
		throw Engine::Exception("Column \"[qz]\" has not been found in "+
				m_programSettings->getEngineSettings().datafile);
	}

	/*allocate arguments and residuals*/
	for(size_t i = 0; i < m_exp_intens_vals.size(); ++i)
	{
		m_DataPoints.push_back(
				NonlinearFit::DataPoint(new ANACalculatorCoplanarTripleArgument(m_qx_vals[i], m_qz_vals[i]),
						m_exp_intens_vals[i]));
	}

	m_ini_intens_vals.resize(m_exp_intens_vals.size(), 0.0);
	m_fin_intens_vals.resize(m_exp_intens_vals.size(), 0.0);

	std::cout << "Nb data residuals:\t" << m_DataPoints.size() << std::endl;
}

void Engine::setupCalculator()
{
    Vector3d Q_vec, Q;
	Vector3d b_vec, l_vec, n_vec, b;
	MillerCubIndicesTransformator transformator(
	                                m_programSettings->getSampleSettings().a0);
	double phi;
	
	l_vec = transformator.toVector3d(m_programSettings->getSampleSettings().misfit.l);
	b_vec = transformator.toVector3d(m_programSettings->getSampleSettings().misfit.b);
	n_vec = Vector3d(0, 0, 1);
	b = toMisfitFrame(b_vec, l_vec, n_vec);

	/*transform hexagonal miller indices to vector*/
	Q_vec = transformator.toVector3d(m_programSettings->getCalculatorSettings().Q);
	/*get Q in coplanar frame*/
	Q = toCoplanarFrame(Q_vec, n_vec);
	/* sign of the first index of Q_vec determines the difference between reflections
	 * like [224] and [-2-24]
	*/
	Q[0] *= GSL_SIGN (Q_vec[0]);
	
	phi = toMisfitInPlaneAngle(Q_vec, b_vec, l_vec, n_vec);
                
	try
	{
		m_sample = new ANASampleCub(m_programSettings->getSampleSettings().thickness,
				m_programSettings->getSampleSettings().width);

		/*one misfit interfaces*/
		m_sample->addMisfitInterface(
					m_programSettings->getSampleSettings().misfit.rho.m_Value * 1e-7,
					b(0), b(1), b(2),
					Q(0), Q(1), Q(2),
					phi,
					m_programSettings->getSampleSettings().nu,
					m_programSettings->getSampleSettings().thickness);
		/*add threading layers*/
		m_sample->addThreadingLayer(
					m_programSettings->getSampleSettings().threading_edge.rho.m_Value * 1e-14,
					m_programSettings->getSampleSettings().threading_edge.b_edge,
					m_programSettings->getSampleSettings().threading_edge.b_screw,
					m_programSettings->getSampleSettings().threading_edge.rc.m_Value, Q(0), Q(2),
					m_programSettings->getSampleSettings().nu);
		m_sample->addThreadingLayer(
					m_programSettings->getSampleSettings().threading_screw.rho.m_Value * 1e-14,
					m_programSettings->getSampleSettings().threading_screw.b_edge,
					m_programSettings->getSampleSettings().threading_screw.b_screw,
					m_programSettings->getSampleSettings().threading_screw.rc.m_Value, Q(0), Q(2),
					m_programSettings->getSampleSettings().nu);
		m_sample->addThreadingLayer(
					m_programSettings->getSampleSettings().threading_mixed.rho.m_Value * 1e-14,
					m_programSettings->getSampleSettings().threading_mixed.b_edge,
					m_programSettings->getSampleSettings().threading_mixed.b_screw,
					m_programSettings->getSampleSettings().threading_mixed.rc.m_Value, Q(0), Q(2),
					m_programSettings->getSampleSettings().nu);

		m_calculator = new ANACalculatorCoplanarTriple(m_sample,
				m_programSettings->getCalculatorSettings().sampling);
		m_calculator->setResolution(
				m_programSettings->getCalculatorSettings().qresolX,
				m_programSettings->getCalculatorSettings().qresolZ);

		m_fit_calculator = new FitANACalculatorCoplanarTriple(m_calculator, m_sample);
	} catch (std::exception& ex)
	{
		throw Engine::Exception(
				"Calculator allocation pb (" + toString(ex.what()) + ")");
	}
}

void Engine::addFParameter(const FitParameter & param)
{
	if (param.m_Lbvalue != param.m_Ubvalue)
	{
		m_fParametersInit.push_back(param);
	}
}

void Engine::setupCParameters()
{
	const FitParameter * param;
	/*scale*/
	param = &m_programSettings->getCalculatorSettings().scale;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*background*/
	param = &m_programSettings->getCalculatorSettings().background;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*densities*/
	/*edge*/
	param = &m_programSettings->getSampleSettings().threading_edge.rho;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*screw*/
	param = &m_programSettings->getSampleSettings().threading_screw.rho;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*mixed*/
	param = &m_programSettings->getSampleSettings().threading_mixed.rho;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*critical radii*/
	/*edge*/
	param = &m_programSettings->getSampleSettings().threading_edge.rc;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*screw*/
	param = &m_programSettings->getSampleSettings().threading_screw.rc;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*mixed*/
	param = &m_programSettings->getSampleSettings().threading_mixed.rc;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);

	/*misfit*/
	param = &m_programSettings->getSampleSettings().misfit.rho;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
}

void Engine::setupFitter()
{
	m_fitter = new PortFitter;
	std::cout << "nb FitParameters:\t" << m_fParametersInit.size() << std::endl;

	m_fitter->init(m_fit_calculator, m_fParametersInit, m_cParameters, m_DataPoints);
}

void Engine::printFitterInfo(NonlinearFit::NonlinearFitter * fitter)
{
	//size_t index;

	std::cout << "In " << fitter->getNbIterations() << " iterations ||f||^2 reduced from "
				<< fitter->getFinit() << " to " << fitter->getFfin() << std::endl;
	std::cout << "Nb Function evaluations:\t" << fitter->getNbFuncEval() << std::endl;
	std::cout << "Nb Jacobian evaluations:\t" << fitter->getNbGradEval() << std::endl;
	std::cout << "Reason to stop iterations:\t<" << fitter->getReasonToStop()
			<< "> Code:\t" << fitter->getReasonToStopId() << std::endl;

	std::cout << "Fitted parameters:" << std::endl;
	if(fitter->getName().compare("LEVMAR") == 0)
	{
		//index = 0;

		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "[" << i << "]\t" << fitter->getFitParameter(i).m_Name
					<< "\t" << fitter->getFitParameter(i).m_Value << " +/- "
					<< sqrt(fitter->getCovarMatrix(i, i)) << std::endl;
		}

		std::cout << "Correlation matrix:" << std::endl;
		std::cout << "\t";
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "["<< i <<"]" << "\t";
		}
		std::cout << std::endl;
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "["<< i <<"]" << "\t";
			for(size_t j = 0; j <= i; ++j)
			{
				std::cout
						<< fitter->getCovarMatrix(i, j)
								/ sqrt(
										fitter->getCovarMatrix(i, i)
												* fitter->getCovarMatrix(j, j))
						<< "\t";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "[" << i << "]\t" << fitter->getFitParameter(i).m_Name
					<< "\t" << fitter->getFitParameter(i).m_Value << std::endl;
		}
		std::cout << "Covariance matrix has not been calculated." << std::endl;
	}
}
