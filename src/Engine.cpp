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
	m_fitter = NULL;
	m_programSettings = NULL;
	m_fit_calculator = NULL;
	m_sample = NULL;
}

Engine::~Engine()
{
	if(m_programSettings)
		delete m_programSettings;

    init();
}

void Engine::exec(const boost::filesystem::path & workDir)
{
   
    m_WorkDir = workDir;
	m_programSettings = new ProgramSettings;
	m_programSettings->read(m_WorkDir);

	init();
	setupComponents();
	doWork();
	saveResult();
	
	//update configuration files
	updateConfigFile(m_programSettings->getSampleConfigFile());
	updateConfigFile(m_programSettings->getDataConfigFile());
}

void Engine::init()
{
	for(size_t i = 0; i < m_calculators.size(); ++i)
	{
        if(m_calculators[i])
            delete m_calculators[i];
	}
	if(m_fit_calculator)
		delete m_fit_calculator;
	if(m_fitter)
		delete m_fitter;

	/*vector of fit parameters with their final values*/
	m_fParametersFinal.clear();

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
}

void Engine::calcI(std::vector<double>& vals)
{
	double I;

	m_fit_calculator->reinit(m_programSettings->getCPMap());

	for(size_t ipt = 0; ipt < m_DataPoints.size(); ipt++)
	{
		I = m_fit_calculator->eval(m_DataPoints[ipt].m_Argument);
		vals[ipt] = I;
	}
}

boost::filesystem::path
Engine::getOutFilename(const boost::filesystem::path& infile) const
{
    boost::filesystem::path outfile;
    
    outfile = infile.filename();
	outfile.replace_extension("ft");
	return m_WorkDir / outfile;
}

void Engine::saveResult() const
{
    boost::filesystem::path filename;
    
    filename = getOutFilename(m_programSettings->getDataConfig(0).file);
	std::ofstream fout(filename.c_str());
	if(!fout)
	{
		throw Engine::Exception("Unable to open the output file:\t" +
		                    filename.native());
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
	m_fitter->fit(NonlinearFitter::fitLIN,
	    m_programSettings->getFitConfig().nbIterations);

	//get new values (best fit values) of fit parameters
	m_fParametersFinal = m_fitter->getFitParameters();
	
	/*update cpMap parameters with newly fitted values*/
	mergeParameters(m_programSettings->getCPMap(),
		m_fParametersFinal);
}

void Engine::saveResume() const
{
	boost::filesystem::path filename;
	std::ofstream fout;

	/*transform rc to M = rc/rd = rc*rho^(1/2) */
	filename = m_WorkDir / m_programSettings->getResumefile();

    fout.open(filename.c_str());
    fout << "In " << m_fitter->getNbIterations() 
        << " iterations ||f||^2 reduced from "
        << m_fitter->getFinit() << " to " 
        << m_fitter->getFfin() << std::endl;
    fout << std::endl;
    for (NonlinearFit::FitParameterList::const_iterator it=m_fParametersFinal.begin();
            it!=m_fParametersFinal.end(); ++it)
    {
        fout << it->m_Name << " => " << it->m_Value << std::endl;
    }
    fout.close();
}

void Engine::updateConfigFile(const boost::filesystem::path & cfgfile) const
{
	libconfig::Config cfg;
	boost::filesystem::path oldcfgfile;

	oldcfgfile = cfgfile;
	oldcfgfile.replace_extension(".~cfg");
	/*save old configuration in a file with extension ~cfg*/
	boost::filesystem::rename(cfgfile, oldcfgfile);

	try
	{
		// Read the configuration file. If there is an error, report it
		cfg.readFile(oldcfgfile.c_str());

		//reset configuration putting the new fitted values from m_fParametersFinal
		for(size_t i = 0; i < m_fParametersFinal.size(); ++i)
		{
		    if(cfg.exists(m_fParametersFinal[i].m_Name))    
                cfg.lookup(m_fParametersFinal[i].m_Name) = 
                                                m_fParametersFinal[i].m_Value;
		}
		//save new configuration
		cfg.writeFile (cfgfile.c_str());
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
    for(size_t id = 0; id < m_programSettings->getDataConfig().size(); ++id)
    {
        setupCalculator(id);
    }
	readData();
	setupFitter();
}

void Engine::readData()
{
	DataReader dr("\t ");
	boost::filesystem::path datapath;
	
	datapath = m_WorkDir / m_programSettings->getDataConfig(0).file;
	dr.parse(datapath.native());

	if(dr.columnExists("[intensity]"))
	{
		m_exp_intens_vals = dr.columnGet("[intensity]");
	}
	else
	{
		throw Engine::Exception("Column \"[intensity]\" has not been found in " +
				m_programSettings->getDataConfig(0).file.native());
	}

	if (dr.columnExists("[qx]"))
	{
		//get points without any transformation
		m_qx_vals = dr.columnGet("[qx]");
	}
	else
	{
		throw Engine::Exception("Column \"[qx]\" has not been found in "+
				m_programSettings->getDataConfig(0).file.native());
	}

	if (dr.columnExists("[qz]"))
	{
		//get points without any transformation
		m_qz_vals = dr.columnGet("[qz]");
	}
	else
	{
		throw Engine::Exception("Column \"[qz]\" has not been found in "+
				m_programSettings->getDataConfig(0).file.native());
	}

	/*allocate arguments and residuals*/
	for(size_t i = 0; i < m_exp_intens_vals.size(); ++i)
	{
		m_DataPoints.push_back(
				NonlinearFit::DataPoint(
				        new ANACalculatorCoplanarTripleArgument(
				            m_qx_vals[i], m_qz_vals[i],
				            0 //FIXME
				            ),
						m_exp_intens_vals[i]));
	}

	m_ini_intens_vals.resize(m_exp_intens_vals.size(), 0.0);
	m_fin_intens_vals.resize(m_exp_intens_vals.size(), 0.0);

	std::cout << "Nb data residuals:\t" << m_DataPoints.size() << std::endl;
}

void Engine::setupCalculator(size_t id)
{
    Vector3d Q_vec, Q;
	Vector3d b_vec, l_vec, n_vec, b;
	MillerCubIndicesTransformator transformator(
	                                m_programSettings->getSampleConfig().a0);
	double phi;
	
	l_vec = transformator.toVector3d(m_programSettings->getSampleConfig().misfit.l);
	b_vec = transformator.toVector3d(m_programSettings->getSampleConfig().misfit.b);
	n_vec = Vector3d(0, 0, 1);
	b = toMisfitFrame(b_vec, l_vec, n_vec);

	/*transform hexagonal miller indices to vector*/
	Q_vec = transformator.toVector3d(m_programSettings->getDataConfig(id).Q);
	/*get Q in coplanar frame*/
	Q = toCoplanarFrame(Q_vec, n_vec);
	/* sign of the first index of Q_vec determines the difference between reflections
	 * like [224] and [-2-24]
	*/
	Q[0] *= GSL_SIGN (Q_vec[0]);
	
	phi = toMisfitInPlaneAngle(Q_vec, b_vec, l_vec, n_vec);
                
	try
	{
		m_sample = new ANASampleCub(m_programSettings->getSampleConfig().thickness,
				m_programSettings->getSampleConfig().width);

		/*one misfit interfaces*/
		m_sample->addMisfitInterface(
					m_programSettings->getSampleConfig().misfit.rho * 1e-7,
					b(0), b(1), b(2),
					Q(0), Q(1), Q(2),
					phi,
					m_programSettings->getSampleConfig().nu,
					m_programSettings->getSampleConfig().thickness);
		/*add threading layers*/
		m_sample->addThreadingLayer(
					m_programSettings->getSampleConfig().threading_edge.rho * 1e-14,
					m_programSettings->getSampleConfig().threading_edge.b_edge,
					m_programSettings->getSampleConfig().threading_edge.b_screw,
					m_programSettings->getSampleConfig().threading_edge.rc, Q(0), Q(2),
					m_programSettings->getSampleConfig().nu);
		m_sample->addThreadingLayer(
					m_programSettings->getSampleConfig().threading_screw.rho * 1e-14,
					m_programSettings->getSampleConfig().threading_screw.b_edge,
					m_programSettings->getSampleConfig().threading_screw.b_screw,
					m_programSettings->getSampleConfig().threading_screw.rc, Q(0), Q(2),
					m_programSettings->getSampleConfig().nu);
		m_sample->addThreadingLayer(
					m_programSettings->getSampleConfig().threading_mixed.rho * 1e-14,
					m_programSettings->getSampleConfig().threading_mixed.b_edge,
					m_programSettings->getSampleConfig().threading_mixed.b_screw,
					m_programSettings->getSampleConfig().threading_mixed.rc, Q(0), Q(2),
					m_programSettings->getSampleConfig().nu);

		m_calculators.push_back(new ANACalculatorCoplanarTriple(m_sample, 150
				/*FIXME : make sampling a part of program settings
				 m_programSettings->getCalculatorSettings().sampling*/));
		m_calculators.back()->setResolution(
				m_programSettings->getDataConfig(id).resolX,
				m_programSettings->getDataConfig(id).resolZ);

		m_fit_calculator = new FitANACalculatorCoplanarTriple(m_calculators, m_sample);
	} catch (std::exception& ex)
	{
		throw Engine::Exception(
				"Calculator allocation pb (" + toString(ex.what()) + ")");
	}
}

void Engine::setupFitter()
{
	m_fitter = new PortFitter;

	m_fitter->init(m_fit_calculator,
                    m_programSettings->getFitConfig().fitParameters,
                    m_programSettings->getCPMap(),
                    m_DataPoints);
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
