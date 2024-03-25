/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	Feb. 2022 */

/** \file Lattice.h
 * \brief Header file for a derived class of ConfigGen, called Lattice.
 * They can be used as initial configurations.
	//What can we generate?
	-----  
	** lattice packing
	** Perturbed lattices
	** Imperfect lattices with vacancies
	** Paracrystals (in progress)
 */

#include "../GenerateConfigs.h"

/** Give basis vectors of well-known lattices
 * @param[in] latticeType	The name of Bravais lattice, e.g., 
								1D: Integer
								2D: Square, Triangle,
								3D: Cubic, BCC_cubic, FCC, FCC_cubic	 BCC_symmetric,
 * @param[out] basis	Basis vectors of corresponding Bravais lattice.
*/
void GetPredefinedBasis(const std::string & latticeType, std::vector<GeometryVector> & basis);

/** Give basis vectors and positions of atoms for some predefined lattices (including some Bravais lattices).
 *  @param[in] crystalName It includes some non-Bravais lattices, such as
 *				2D: Honeycomb, Kagome
 *				3D: Diamond.
 *	@param[out] basis	Basis vectors
 *  @param[out] atoms	Relative coordinates of particles (except for one at (0,0,0)). */
void GetPredefinedCrystals(const std::string & crystalName, std::vector<GeometryVector> & basis, std::vector<GeometryVector> & atoms);



/**	/brief A class to construct a periodic point pattern or sphere packing.
 *  It stores a fundamental cell of a sphere packing. */
class LatticeGen : public ConfigGen
{
	std::vector<size_t> NumInSide;	//< NumInSide[i] = the number of unit cells in the ith direction, for i=0, ..., d-1.

protected:
	double maxPhi;					//< maximal packing fraction. 
	SpherePacking F_Cell;			//<	A fundamental cell
//	virtual void SetBasis(const std::vector<GeometryVector> & basis);
public:
	/** A default constructor. */
	LatticeGen() : ConfigGen() { this->maxPhi = 0; this->F_Cell = SpherePacking(); };
	
	/** A constructor by a SpherePacking.
	 * If particle radii = 0, particles are considered "monodisperse" and the radius are set to be the maximum. */
	LatticeGen(const SpherePacking & unitcell) {
		NumInSide = std::vector<size_t>(unitcell.GetDimension(), 1);
		F_Cell = SpherePacking(unitcell);
		N = F_Cell.NumParticle();
		rho = N / F_Cell.PeriodicVolume();
		d = F_Cell.GetDimension();

		//Particle radii = 0.5 * nearest neighbor distance.
		{
			bool uninitializedPrtRadii = false;		//True = at least one particle has a nonpositive radius.
			for (int i = 0; i < this->N; i++)
				uninitializedPrtRadii |= (F_Cell.GetCharacteristics(i) <= 0.0);

			if (F_Cell.CheckOverlap()) {
				std::cout << "Particles are overlapping. Their sizes are adjusted.\n";
				uninitializedPrtRadii = false;
			}

			if (uninitializedPrtRadii) {
				std::cout	<< "Paritlce radii are not initialized."
							<<"Particle radii are set to be nearest neighbor distances.\n";
				/* at least one particle has a nonpositive radius 
				=> reset particle radii such that they are half of nearest neighbor distances.*/
				for (int i = 0; i < this->N; i++)
					F_Cell.UpdateCharacteristics(0.5*F_Cell.NearestParticleDistance(i), i);
			}
		}
		
		maxPhi = this->UpdateMaxPhi();	
	}

	/** A constructor by name of lattices. 
	 * @param name	The name of predefined periodic packings. */
	LatticeGen(const std::string & name);

	/** A constructor by an input stream. */
	LatticeGen(std::istream & ifile, std::ostream & ofile);

	/** A copy constructor. */
	LatticeGen(const LatticeGen & source) : LatticeGen(source.F_Cell) {};


	/** A member function to update and return maximum packing fraction.
	 * @return Maximal packing fraction.*/
	double UpdateMaxPhi();

	/** Insert a spherical particle in the fundamental cell.
	 * @param relativeCoordinates	Relative coordinates of a particle.
	 * @param radius	Particle radius. The default value is 0.0. */
	void AddPrt(const GeometryVector & relativeCoorinates, double radius = 0.0);

	/** Remove the ith particle in the fundamental cell. 
	 * @param i The index of a particle that is removed. */
	void RemovePrt(size_t i);

	/** Get basis vectors of a fundamental cell. 
	 * @param[out] basis[i]		The ith basis vector. */
	void GetBasisVectors(std::vector<GeometryVector> & basis) {
		basis.clear();
		for (int i = 0; i < this->d; i++)
			basis.push_back(this->F_Cell.GetBasisVector(i));
	}

	/** Rescale particle size.
	 * @param newPhi	A new packing fraction. */
	virtual void SetPackingFraction(double newPhi);

	/** Rescale the cell size. 
	 * @param newRho	A new number density. */
	virtual void SetRho(double newRho);

	/** Resize particle radii to attain a maximal packing fraction. */
	void Set2MaxPhi();

	/** Set the number of unit cells.
	 *@param numInEachSide The number of unit cells in each direction. */
	void SetNumCells(size_t numInEachSide) {	this->NumInSide = std::vector<size_t>(this->d, numInEachSide);	};

	/** Set the number of unit cells.
	 * @param numInSide[i]	The number of unite cells in the ith direction. */
	void SetNumCells(const std::vector<size_t> & numInSide) { this->NumInSide = std::vector<size_t>(numInSide); };
	
	
	/**	Generate a periodic point pattern.
 	 * @param[out] c	A periodic Configuration that contains N[0] x...x N[d-1] unit cells. */
	virtual void GenerateC(Configuration & c) { c = MultiplicateStructure(Configuration(this->F_Cell, "a"), this->NumInSide); }

	/**	Generate a periodic sphere packing.
	 * @param[out] c	A periodic SpherePacking that contains N[0] x...x N[d-1] unit cells. */
	virtual void GenerateP(SpherePacking & c);

};


/********* Transformation **********
 * Converts a given point process to another by a specified method
 * Require an initial configuration.
 * Default initial system = d-dimensional cubic lattice
 **********************************/

/** /brief A struct to introduce point defects to a given system. 
 *	It stores the mean of point defect concentration = (# of vacancies) / (total # of points ). 
 *  Unnecessary information; dimension, number density, the number of particles */
class GenPtsDefect :ConfigGen {
protected:
	double MeanPointDefects;		//< mean fraction of vacancies
	int dist_type;					//< 0 = fixed number of point defects
									//< 1 = fluctuating number of point defects (follow Poisson distribution).
	int DefectType;					//< 0 = vacancy
									//< 1 = interstitials.
	RandomGenerator RNG;			//< a random number generator	
	/** A member function to get the total number of point defects. 
	 * @param TotalPrts	The total number of particles in the initial configuration.
	 * @return The expected number of point defects.
	 */
	inline size_t GetNumDefects(size_t TotalPrts){
		double mean = TotalPrts * MeanPointDefects;
		if (dist_type == 0) {
			return (int)round(mean);
		}
		else if (dist_type == 1) {
			return RNG.PoissonRNG(mean);
		}
		else {
			std::cerr << "Undefined distributions\n";
			return 0;
		}
	}
public:
	/** A default constructor.
	 * @param mean				The mean of concentration of point defects (<1).
	 * @param distribution_type	Distribution of the number of point defects (0 => fixed, 1 => Poisson distribution).
	 * @param defect_type		Type of point defects (0 = vacancy, 1 = interstitials). 	 */
	GenPtsDefect(double mean, int distribution_type, int defect_type) : ConfigGen(), MeanPointDefects(mean), dist_type(distribution_type), DefectType(defect_type) {
		if (mean > 1.0) {	
			std::cout << "GenPtsDefect:: mean is too high!" << std::endl;	
			return;
		}
		if (distribution_type != 0 || distribution_type != 1) { 
			std::cout << "GenPtsDefect:: undefined distributions" << std::endl;
			return;
		}
		if (defect_type != 0 || defect_type != 1) {
			std::cout << "GenPtsDefect:: undefined defect types" << std::endl;
			return;
		}
	};
	
	/** A copy constructor. */
	GenPtsDefect(const GenPtsDefect & source) : GenPtsDefect(source.MeanPointDefects, source.dist_type, source.DefectType) {}

	/** A constructor via streams. */
	GenPtsDefect(std::istream & ifile, std::ostream & ofile);

	/** A member function to reset */
	virtual void SetRho(double numberdensity) { std::cout << "GenPtsDefect doesn't use number density\n"; }
	virtual void SetNumPrts(int N) { std::cout << "GenPtsDefect doesn't use the number of particles\n"; }

	/** A member function to introduce point defects in a Configuraiton object.
	 * @param[in,out] c	A point configuration in which point defects are introduced.
	 */
	virtual void GenerateC(Configuration &c) {
		if (c.NumParticle() > 0) {
			size_t num_defects = this->GetNumDefects(c.NumParticle());
			if (this->DefectType == 0) {
				//Point vacancies;
				int idx = 0;
				for (size_t i = 0; i < num_defects; i++) {
					idx = (int)floor(c.NumParticle()*RNG.RandomDouble());
					c.DeleteParticle(idx);
				}
			}
			else if (this->DefectType == 1) {
				//Interstitial
				GeometryVector rel_pos(static_cast<int>(c.GetDimension()));
				for (size_t i = 0; i < num_defects; i++) {
					for (int j = 0; j < c.GetDimension(); j++)
						rel_pos.x[j] = RNG.RandomDouble();

					c.Insert("a", rel_pos);
				}
			}
			else {
				std::cerr << "Undefined defect_type" << std::endl;	return;
			}
		}
		else {
			std::cout << "c should have at least 1 particle." << std::endl;
			return;
		}
	}

	/** A memebr function to introduce point defects in a SpherePacking object.
	 */
	virtual void GenerateP(SpherePacking & c) {
		std::cout << "GenPtsDefects::GenerateP is not implemented.\n";
		return;
	}
	/** Reset a seed. */
	void SetSeed(int seed) { this->RNG.seed(seed); }
};


/** \brief A class to stochastically and independently displace particle positions by a given distribution!
	This class doesn't require "the number of particles", "number density"
*/
class GenPerturb : public ConfigGen
{
protected:
	std::vector<double> parameters;			//< parameters for iCMD
	std::function<double (double) > iCMD;	//< the inverse of cummulative function of radial displacements.
											//< iCMD: [0,1]-> [0, infty].
	std::function<GeometryVector(RandomGenerator &)> GetRandomDispl;	//< A Lambda function to get a random displacement.
	RandomGenerator rng;					//< a number number generator
public:
	/** A constructor. */
	GenPerturb(size_t dimension, const std::function <double(double)> & inverseCMD, int seed = 0) : ConfigGen() {
		this->d = dimension;
		this->iCMD = inverseCMD;
		this->rng.seed(seed);

		this->GetRandomDispl 
			= [this](RandomGenerator &)->GeometryVector 
		{
			GeometryVector pos = RandomUnitVector(this->d, this->rng);
			pos = this->rng.RandomDouble(this->iCMD) * pos;
			return pos;
		};
	}

	/** A constructor via an input stream.
	 * Predefined distributions = uniform, gaussian, powerlaw.
	 * If Type of distribution = etc => it requries the inverse of cumulative distribution function.
	 * @param inverseCMD	The inverse function of cumulative distribution of radial displacement.  */
	GenPerturb(std::istream & ifile, std::ostream & ofile, const std::function <double(double)> & inverseCMD = nullptr) : ConfigGen() {
		ofile	<< " ----------------------------------\n "
				<< "             GenPerturb            \n"
				<< " ----------------------------------\n" << std::endl;
		ofile << "Dimension = ";
		ifile >> this->d;
		
		std::string name;
		ofile << "Type of distribution = ";
		ifile >> name;
		std::transform(name.begin(), name.end(), name.begin(), ::tolower);
		if (strcmp(name.c_str(), "uniform")==0) {
			double maxRadius;
			ofile << "Maximal radius (positive real) = ";
			ifile >> maxRadius;
			parameters.push_back(maxRadius);

			this->GetRandomDispl = [this](RandomGenerator & rng)->GeometryVector {
				GeometryVector result = RandomUnitVector(this->d, rng);
				double r;

				switch (this->d)
				{
				case 1:
					r = this->rng.RandomDouble();	break;
				case 2:
					r = sqrt(this->rng.RandomDouble());	break;
				case 3:
					r = cbrt(this->rng.RandomDouble()); break;
				default:
					r = pow(this->rng.RandomDouble(), 1.0 / this->d);
				}

				return this->parameters[0]* r *result;
			};
		}
		else if (strcmp(name.c_str(), "gaussian") == 0) {
			double standard_deviation;
			ofile << "Input the standard deviation (positive real)";
			ifile >> standard_deviation;
			parameters.push_back(standard_deviation);

			this->GetRandomDispl = [this](RandomGenerator & g)->GeometryVector {
				GeometryVector x(static_cast<int>(this->d));
				for (int i = 0; i < x.Dimension; i++) {
					x.x[i] = g.RandomDouble_normal(this->parameters[0]);
				}
				return x;
			};
		}
		else if (strcmp(name.c_str(), "powerlaw") == 0|| strcmp(name.c_str(), "power-law") == 0) {
			double delta, alpha;
			std::cout << "Input the length scale (positive real)";
			std::cin >> delta;
			std::cout << "Input the power exponent (positive real)";
			std::cin >> alpha;
			this->parameters.push_back(delta);	this->parameters.push_back(alpha);

			this->iCMD = [this](double y)->double {
				double delta = this->parameters[0], alpha = this->parameters[1];
				if (y  <  0.5 / (1.0 + alpha))
				{
					return  -delta / pow(2.0*(1.0 + alpha)*y, 1.0 / alpha);
				}
				else if (y  <  0.5*(1.0 + alpha / (1.0 + alpha)))
				{
					return  (2.0*y - 1)*(1 + alpha) / alpha * delta;
				}
				else
				{
					return  delta / pow(2.0*(1.0 - y)*(1.0 + alpha), 1.0 / alpha);
				}
			};

			/*auto f = [delta, alpha, dimension](double y)->double {
				if(y  <  0.5 / (1.0 + alpha)) 
				{	return  -delta / pow(2.0*(1.0 + alpha)*y, 1.0 / alpha); }
				else if(y  <  0.5*(1.0 + alpha / (1.0 + alpha))) 
				{	return  (2.0*y - 1)*(1 + alpha) / alpha * delta; }
				else 
				{	return  delta / pow(2.0*(1.0 - y)*(1.0 + alpha), 1.0 / alpha); }
			};*/
			this->GetRandomDispl = [this](RandomGenerator & g)->GeometryVector{
				GeometryVector x = RandomUnitVector(this->d, g);
				return g.RandomDouble(this->iCMD)*x;
			};
		}
		else if (strcmp(name.c_str(), "etc") == 0) {
			if (inverseCMD != nullptr) {
				this->iCMD = inverseCMD;
				this->GetRandomDispl = [this](RandomGenerator & g)->GeometryVector {
					return g.RandomDouble(this->iCMD)*RandomUnitVector(this->d, g);
				};
			}
			else {
				ofile << "Distribution is not defined!" << std::endl;
			}
		}		
		else{
			ofile << "Distribution is not defined!" << std::endl;
		}
	}

	/** A copy constructor */
	GenPerturb(const GenPerturb & source) : GenPerturb(source.d, source.iCMD){
		this->parameters = std::vector<double>(source.parameters);
		this->GetRandomDispl = source.GetRandomDispl;
	}
	

	/** Reset random seed. */
	void SetSeed(int seed) { this->rng.seed(seed); 	}

	/** Stochastically displace particles in c.
	 * @param[in,out] c	A point configuration that will be perturbed. It should have at least one particle, and positive volume.	 */
	virtual void GenerateC(Configuration & c) {
		if (c.GetDimension() == this->d) {
			if (c.NumParticle() > 0) {
				if (c.PeriodicVolume() > 0) {
					GeometryVector Displacement_Rel(static_cast<int> (this->d));

					for (size_t i = 0; i < c.NumParticle(); i++) {
						//A random and isotropic displacement in Cartesian coordinates.
						Displacement_Rel = this->GetRandomDispl(this->rng);//this->rng.RandomDouble(this->iCMD) * RandomUnitVector(this->d, this->rng);
						//Convert it into relative coordinates.
						Displacement_Rel = c.CartesianCoord2RelativeCoord(Displacement_Rel);
						//Displace particle i by this dispalcement.
						c.MoveParticle(i, c.GetRelativeCoordinates(i) + Displacement_Rel);
					}
				}
				else {
					std::cout << "c should have positive volume" << std::endl;
					return;
				}
			}
			else {
				std::cout << "c should have at least one particle." << std::endl;
				return;
			}
		}
		else {
			std::cout << "The dimension of c should be identical to the dimension of GenPerturb" << std::endl;
			return;
		}
	}
	
	
	/** Stochastically displace spherical particles in c.
	 * If particle i is overlapping with another at a new position, trial movements are made until particle i is not overlapping with others.
	 * @param[in,out] c	A SpherePacking that will be perutrbed. This member function ends up with a perturbed packing.*/
	virtual void GenerateP(SpherePacking & c) {
		if (c.GetDimension() == this->d) {
			if (c.NumParticle() > 0) {
				if (c.PeriodicVolume() > 0) {
					GeometryVector	TrialPos_Rel(static_cast<int> (this->d)),	//Trial position in relative coordinates.
									OldPos_Rel(static_cast<int> (this->d));		//Initial position in relative coordinates.			

					double succ_rate = 0.0;
					for (size_t i = 0; i < c.NumParticle(); i++) {
						OldPos_Rel = c.GetRelativeCoordinates(i);	//Initial particle position.
						double r = c.GetCharacteristics(i);			//Particle radius.
						
						if (c.NumParticle() > 1) {	//There are multiple particles!!
							//Temporarily set aside particle i to either particle 0 or particle i-1.
							int nextPrtIdx = (i == 0) ? 1 : i - 1;
							c.MoveParticle(i, c.GetRelativeCoordinates(nextPrtIdx));
							c.UpdateCharacteristics(0.0, i);
							//Try to move particle i until it is not overlapping with others.
							do {
								succ_rate++;
								//A random and isotropic displacement in Cartesian coordinates.
								TrialPos_Rel = this->GetRandomDispl(this->rng);  //this->rng.RandomDouble(this->iCMD) * RandomUnitVector(this->d, this->rng);
								//Convert it into relative coordinates, and add initial position.
								TrialPos_Rel = c.CartesianCoord2RelativeCoord(TrialPos_Rel) + OldPos_Rel;
							} while (c.CheckOverlap(TrialPos_Rel, r));
						}
						else {	//There is only one particle.
							//A random and isotropic displacement in Cartesian coordinates.
							TrialPos_Rel = this->GetRandomDispl(this->rng); //this->rng.RandomDouble(this->iCMD) * RandomUnitVector(this->d, this->rng);
							//Convert it into relative coordinates, and add initial position.
							TrialPos_Rel = c.CartesianCoord2RelativeCoord(TrialPos_Rel) + OldPos_Rel;
						}
						//Move particle i into a new position and recover its radius.
						c.MoveParticle(i, TrialPos_Rel);
						c.UpdateCharacteristics(r, i);
					}
					std::cout << "Success rate of displacing particles = " << c.NumParticle() / succ_rate << std::endl;
				}
				else {
					std::cout << "c should have positive volume" << std::endl;
					return;
				}
			}
			else {
				std::cout << "c should have at least one particle." << std::endl;
				return;
			}
		}
		else {
			std::cout << "The dimension of c should be identical to the dimension of GenPerturb" << std::endl;
			return;
		}
	}


	/** Undefined member functions !*/
	virtual void SetRho(double rho) { std::cout << "GenPerturb::SetRho is not defined\n"; };
	virtual void SetNumPrts(int num) { std::cout << "GenPerturb::SetNumPrts is not defined\n"; };
	virtual void SetPackingFraction(double phi) { std::cout << "GenPerturb::SetPackingFraction is not defined\n"; };
};