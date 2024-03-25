/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file GenerateConfigs.h
 * \brief Header file for a ConfigGen class and its derived class.
 * Collections of point configurations and spherepackings that you can simply construct.
 * They can be used as initial configurations.
	//What can we generate?
	----- Lattices/
	** lattice packing
	** Perturbed lattices
	** Imperfect lattices with vacancies
	** Paracrystals (in progress)
	
	----- RSA/
	** Saturated RSA (sphere packings)
	----- HSF/
	** equilibrium hard-sphere fluids
	----- DPP/
	** Determinantal point processes (TODO)
	----- Jamming
	** Dense Packing via LS algorithm (TODO)
 */

#ifndef GENERATECONFIGS
#define GENERATECONFIGS

#include <omp.h>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

/* in $(core) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <etc.h>



/*Compute the number of k-points in spectral densities*/
size_t numK_points(PeriodicCellList<Empty> & FundamentalCell, std::vector<GeometryVector> & ks);

template<typename T> inline void Echo(std::istream & ifile, std::ostream & ofile, T & variable) {
	ifile >> variable;	ofile << variable << "\n";
}


/*
	Callable types to generate configurations
	Instructions to the Constructors
	Ask parameters in the following orders
	dimension -> intensive parameters -> extensive parameters
*/

/** \brief  Abstract class to define a point process.
 */
class ConfigGen
{
protected:
	double rho;		//< number density
	size_t N,	//< (expected) number of particles
		numThreads;	//< the number of threads that are used to construct a system.
	DimensionType d;//< space dimension

public:

	/** A default constructor of ConfigGen. */
	ConfigGen() { this->rho = 1.0; this->N = 0; this->d = 1; this->numThreads = 1; }

	ConfigGen(int dimensions, int N, double rho);

	ConfigGen(std::istream & ifile, std::ostream & ofile);

	/** A pure virtual function to construct a point configuration.
	 * @param[out] c	A point configuration. */
	virtual void GenerateC(Configuration & c) = 0;

	/** A pure virtual function to construct a sphere packing.
	 * @param[out] c	A sphere packing. */
	virtual void GenerateP(SpherePacking & c) = 0;
	
	//Change system parameters 
	/** A virtual member function to change number density (must be redefined).
	 * @param rho	new number density. */
	virtual void SetRho(double rho) { std::cerr << "SetRho() is undefined" << std::endl; }

	/** A virtual member function to change particle number (must be redefined).
	 * @param num	new particle numbers. */
	virtual void SetNumPrts(int num) {	this->N = num; } 
		//std::cerr << "SetNumPrts() is undefined" << std::endl; }

	/** A virtual member function to change the packing fraction (must be redefined).
	 * @param phi	new packing fraction. */
	virtual void SetPackingFraction(double phi) { std::cerr << "SetPackingFraction() is undefined" << std::endl; }	//must be redefined
	
	/** @return the number density.*/
	virtual double GetRho() { return this->rho; }
	
	/** @return the particle number.*/
	virtual size_t GetNumPrts() { return this->N; }
	
	/** Redefin the number of threads that are used to construct realizations.*/
	virtual void SetNumThreads(int numThreads) { 
		this->numThreads = numThreads; 
		omp_set_num_threads(this->numThreads);
		}

	/** A pure virutal member function to input extra parameters. 
	 * @param[in] option	the number of option.
	 * @param[in] input		An input stream.
	 * @param[out] output	An output stream.
	 */
	virtual void SetAdditionalOption(const std::string &option, std::istream & input, std::ostream &output) { output << "Uncrecognized command!\n"; /*by default, accept no additional option*/ }
	
	DimensionType GetDimension() {
		return this->d;
	}

	/** A default destructor. */
	virtual ~ConfigGen() {}
};

///* Generate an empty periodic fundamental cell */
//struct EmptyBox : public ConfigGen 
//{
//protected:
//	Configuration F_Cell;
//	double Cell_Volume;
//	virtual void SetBasis(const std::vector<GeometryVector> & basis);
//public:
//	/*** Constructors **********
//	 * Default Constructor 
//	 * Fully parametrized 
//	 * Name-parametrized
//	 * CMD interfaces
//	 * Copy
//	 ***************************/
//	EmptyBox();
//	EmptyBox(int dimension, const std::vector<GeometryVector> & BasisVectors);
//	EmptyBox(const std::string & CellName, double CellVolume);
//	EmptyBox(std::istream & ifile, std::ostream & ofile);
//	EmptyBox(const EmptyBox & copy);
//
//	void SetCellVolume(double new_cell_volume) { this->F_Cell.Resize(new_cell_volume); };
//
//	virtual void GenerateC(Configuration & c) { c = Configuration(this->F_Cell); };
//	virtual void GenerateP(SpherePacking & c) { c = SpherePacking(this->F_Cell, 0.0); };
//};


/** \breif A class to construct a simple Poisson point pattern. */
class PoissonPtGen : public ConfigGen {
protected:
	int mode = 0;	//<	0: particle number is fixed
					//< 1: particle number can vary
	RandomGenerator rng;	//< a random number generator
	
	inline size_t GetExpectedNumPrts(double mean) {
		if (this->mode == 0) {
			return (int)round(mean);
		}
		else if (this->mode == 1) {
			return rng.PoissonRNG(mean);
		}
		else{
			std::cerr<<"PoissonPtGen::invalid mode!\n";
			return 0;
		}
	}
public:
	/** A default constructor.
	 * @param rho			The average number density.
	 * @param mode(0/1)		0=fixed particle number/ 1=follows Poisson distribution.
	 * @param seed			A seed for a random numnber generator.	 */
	PoissonPtGen(double rho, int mode = 0, int seed = 0) : ConfigGen() {
		this->rho = rho;
		if (mode == 0 || mode == 1) {
			this->mode = mode;
		}
		else { std::cout << "PoissonPtGen::undefined mode!\n";	return; }
		this->rng.seed(seed);
	}

	/** A constructor via an input stream. */
	PoissonPtGen(std::istream & ifile, std::ostream & ofile);
	
	/** Reset the number density. */
	virtual void SetRho(double rho) { this->rho = rho; }

	/** Redefine the mode */
	void SetMode(int mode) { this->mode = mode; }
	
	/** Empty the simulation box and generate a Poisson point pattern of a given number density. 
	 *	@param[in,out] c	A periodic simulation box whose volume should be positive!	 */
	virtual void GenerateC(Configuration & c); 
	virtual void GenerateP(SpherePacking & c); 


	/** Undefined member functions !! */
	//virtual void SetNumPrts(int num);
	virtual void SetPackingFraction(double phi) { std::cout << "PoissonPtgen::packing fraction is unnecessary\n"; }
	
};

/** \brief A class to help to load and save a realization under the periodic boundary condition.
 *	This class currently supports 
		.txt	->    ReadConfiguration
	    .ConfigPack -> ConfigurationPack
		.SpherePack -> SpherePackingPack	*/
class LoadConfig : public ConfigGen {
protected:
	std::string LoadPrefix;
	std::string Format;

public:
	/* Load configurations
	 * numConfig	the number of repetitions, when prefix includes specifier (%d), searches from 0~numConfig.
	 * dataType		class type of the loaded data
					<0: no change	
					0: Configuration
					1: SpherePacking
					2: Dispersions (if possible)
	 * rho			number density of the load configuration.
					<0: no change for the original data.
	 * phi			packing fraction of the load configuration (only particle sizes are changed)
					<0
	*/
	//TODO: implementation of LoadConfig()
	LoadConfig(const std::string & prefix, const std::string & fileformat, int numConfig, int dataType,double rho =-1.0, double phi = -1.0);
	LoadConfig(std::istream & ifile, std::ostream & ofile);

//	virtual void GenerateC(Configuration & c);
//	virtual void GeenrateP(SpherePacking & c);

};




/*a class to use the GenerateConfig class in CLI
	there are three modes
		Load files
		Empty box
		Generate Configurations
*/
class GenerateConfigCLI {
private:
	ConfigGen * Model;
	std::string ModelName;
	std::function <Configuration(size_t)> getpointconfig;
	std::function <SpherePacking(size_t)> getpackingconfig;
public:
	GenerateConfigCLI() {
		Model = nullptr;
		ModelName = "";
		getpointconfig = nullptr;
		getpackingconfig = nullptr;
	}
	
	GenerateConfigCLI(std::istream & ifile, std::ostream & ofile) { SetModel(ifile, ofile); }

	//change the model that it will generate
	void SetModel(std::istream &ifile, std::ostream & ofile) {
		//Avaiable options
		std::string temp;
		//short count = 0;
		
		ofile << " Avaiable options (Load/EmptyBox/Generate)\n";
		ifile >> temp;
		std::transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
		if ((strcmp(temp.c_str(), "Generate") == 0) || (strcmp(temp.c_str(), "GenerateConfig") == 0)) {


		}
		else {
			if ((strcmp(temp.c_str(), "Load") == 0) || (strcmp(temp.c_str(), "LoadPrefix") == 0)) {

			}
			else if ((strcmp(temp.c_str(), "EmptyBox") == 0) || (strcmp(temp.c_str(), "EmptyUnitBox") == 0)) {

			}
			this->Model = nullptr;
		}
	
	}

	~GenerateConfigCLI() {
		if(Model != nullptr)
			delete Model;
	}

	void GetConfigC(Configuration &c, size_t i) { c = this->getpointconfig(i); }
	void GetConfigP(SpherePacking &p, size_t i) { p = this->getpackingconfig(i); }

	void SaveConfig(Configuration &c);
	void SaveConfig(SpherePacking &c);

};


#endif