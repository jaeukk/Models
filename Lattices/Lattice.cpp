#include "Lattice.h"

/** Give basis vectors of some Bravais lattices. */
void GetPredefinedBasis(const std::string & lattice, std::vector<GeometryVector> & basis) {
	basis.resize(0);

	if (strcmp(lattice.c_str(), "integer") == 0 || strcmp(lattice.c_str(), "Integer") == 0) {
		basis.emplace_back(GeometryVector(static_cast<double>(1.0)));
	}
	else if (strcmp(lattice.c_str(), "square") == 0 || strcmp(lattice.c_str(), "Square") == 0) {
		basis.emplace_back(GeometryVector(1.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 1.0));
	}
	else if (strcmp(lattice.c_str(), "triangle") == 0 || strcmp(lattice.c_str(), "triangular") == 0 || strcmp(lattice.c_str(), "Triangular") == 0) {
		basis.emplace_back(GeometryVector(1.0, 0.0));
		basis.emplace_back(GeometryVector(cos(pi / 3.0), sin(pi / 3.0)));
	}
	else if (strcmp(lattice.c_str(), "cubic") == 0||strcmp(lattice.c_str(), "Cubic") == 0) {
		basis.emplace_back(GeometryVector(1.0, 0.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 1.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 0.0, 1.0));
	}
	else if (strcmp(lattice.c_str(), "fcc") == 0 || strcmp(lattice.c_str(), "FCC") == 0) {//symmetric form
		basis.emplace_back(GeometryVector(0.0, 0.5, 0.5));
		basis.emplace_back(GeometryVector(0.5, 0.0, 0.5));
		basis.emplace_back(GeometryVector(0.5, 0.5, 0.0));
	}
	else if (strcmp(lattice.c_str(), "bcc") == 0 || strcmp(lattice.c_str(), "BCC") == 0) {//symmetric form
		basis.emplace_back(GeometryVector(-0.5, 0.5, 0.5));
		basis.emplace_back(GeometryVector(0.5, -0.5, 0.5));
		basis.emplace_back(GeometryVector(0.5, 0.5, -0.5));
	}
	else if (strcmp(lattice.c_str(), "fcc_cubic") == 0 || strcmp(lattice.c_str(), "FCC_cubic") == 0) {//non-symmetric form
		basis.emplace_back(GeometryVector(1.0, 0.0, 0.0));
		basis.emplace_back(GeometryVector(0.5, sqrt(3.0)*0.5, 0.0));
		basis.emplace_back(GeometryVector(0.5, 0.5 / sqrt(3.0), sqrt(2.0 / 3.0)));
	}
	else if (strcmp(lattice.c_str(), "bcc_cubic") == 0 || strcmp(lattice.c_str(), "BCC_cubic") == 0) {//Non-symmetric form
		basis.emplace_back(GeometryVector(1.0, 0.0, 0.0));
		basis.emplace_back(GeometryVector(0.0, 1.0, 0.0));
		basis.emplace_back(GeometryVector(0.5, 0.5, 0.5));
	}

	//else {
	//	std::cout << "Incorrect cell type" << std::endl;
	//}
}

/** Give basis vectors and positions of atomes for predefined periodic point patterns.*/
void GetPredefinedCrystals(const std::string & crystalName, std::vector<GeometryVector> & basis, std::vector<GeometryVector> & atoms) {
	atoms.resize(0);

	GetPredefinedBasis(crystalName, basis);
	//search for nonBravais crystals
	if (basis.size() == 0) {
		if (strcmp(crystalName.c_str(), "honeycomb") == 0 || strcmp(crystalName.c_str(), "HoneyComb") == 0) {
			GetPredefinedBasis("triangular", basis);
			atoms.emplace_back(GeometryVector(1.0/3.0, 1.0/3.0));
		}
		else if (strcmp(crystalName.c_str(), "kagome") == 0 || strcmp(crystalName.c_str(), "Kagome") == 0) {
			GetPredefinedBasis("triangular", basis);
			atoms.emplace_back(GeometryVector(0.5, 0.0));
			atoms.emplace_back(GeometryVector(0.0, 0.5));
		}
		else if (strcmp(crystalName.c_str(), "diamond") == 0 || strcmp(crystalName.c_str(), "Diamond") == 0) {
			GetPredefinedBasis("fcc", basis);
			atoms.emplace_back(GeometryVector(0.25, 0.25, 0.25));
		}
		else if (strcmp(crystalName.c_str(), "hcp") == 0 || strcmp(crystalName.c_str(), "HCP") == 0) {
			basis.clear();
			basis.emplace_back(0.5, -0.5*sqrt(3.0), 0.0);
			basis.emplace_back(0.5, +0.5*sqrt(3.0), 0.0);
			basis.emplace_back(0.0, 0.0, 2.0*sqrt(2.0/3.0));

			atoms.emplace_back(1.0 / 3.0, 1.0 / 3.0, 0.50);
		}
		else {
			std::cout << "Undefined crystal structure "<<std::endl;
		}
	}
}



/***********************************
	LatticeGen: Constructors
***********************************/

/** Construct a LatticeGen object by lattice name. 
 *	@param name (case-insenstive)	1D:	Integer
									2D: Square, Triangular (triangle), Honeycomb, Kagome
									3D: Cubic, BCC_cubic, FCC, FCC_cubic, BCC_symmetric, Diamond */
LatticeGen::LatticeGen(const std::string & Name) {
	//Change all alphabets to lower cases.
	std::string name(Name);
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);	
	//Obtain basis vectors and positions of atoms
	std::vector<GeometryVector> basis, atoms;
	GetPredefinedCrystals(name, basis, atoms);
	this->F_Cell = SpherePacking((DimensionType)basis.size(), &basis[0], 1.0);

	//a particle at the origin
	this->F_Cell.Insert(0.0, GeometryVector(static_cast<int>(basis.size())));
	//Insert other particles
	for (auto rel_pos = atoms.begin(); rel_pos != atoms.end(); rel_pos++)
		this->F_Cell.Insert(0.0, *rel_pos);

	this->d	= this->F_Cell.GetDimension();
	this->N = this->F_Cell.NumParticle();
	this->rho = this->N / this->F_Cell.PeriodicVolume();
	this->SetNumCells(1);

	//Max packing fraction
	this->Set2MaxPhi();
}

/** Construct a LatticeGen object by an input stream. 
 *	@param[in] ifile	An input stream.
 *	@param[out] ofile	An output stream.
 */
LatticeGen::LatticeGen(std::istream & ifile, std::ostream &ofile) {
	char temp[300] = {};
	std::string cell_type;	double value;	size_t NumInEachSide;
	ofile << "-----------------------------------\n"
		  << "	Choose a lattice to generate     \n"
		  << "-----------------------------------" << std::endl;
	ofile << "Cell Type = ";	Echo<std::string>(ifile, ofile, cell_type);
	{//Define a fundamental cell!
	 //Change all alphabets to the lower cases.
		std::string name(cell_type);
		std::transform(name.begin(), name.end(), name.begin(), ::tolower);
		//Obtain basis vectors and positions of atoms
		std::vector<GeometryVector> basis, atoms;
		GetPredefinedCrystals(name, basis, atoms);
		this->F_Cell = SpherePacking(basis.size(), &basis[0], 1.0);

		//a particle at the origin
		this->F_Cell.Insert(0.0, GeometryVector(static_cast<int>(basis.size())));
		//Insert other particles
		for (auto rel_pos = atoms.begin(); rel_pos != atoms.end(); rel_pos++) {
			this->F_Cell.Insert(0.0, *rel_pos);			
		}


		this->d = this->F_Cell.GetDimension();
		
		ofile << "Add additional particles (y/n)?";
		ifile >> temp;
		if (strcmp(temp, "y") == 0) {
			GeometryVector pos(static_cast<int>(this->d));
			size_t NumAdditionalPrts = 0, 
				NumInit = this->F_Cell.NumParticle();
			ofile << "How many would you add?";
			ifile >> NumAdditionalPrts;
			for (size_t i = 0; i < NumAdditionalPrts; i++) {
				ofile << "Please input relative coordinates of particle " << i+NumInit <<":\t";
				for (int j = 0; j < this->d; j++)
					ifile >> pos.x[j];

				this->AddPrt(pos);
			}
		}
		this->N = this->F_Cell.NumParticle();
		
		//Change particle radii to be maximized.
		this->Set2MaxPhi();
	}
	
	//Resize the fundamental cell and particles.
	{
		ofile << "Rescale the fundamental cell with respect to \n"
			  << "\t D (nearest neighbor distance) \n"
			  << "\t rho (number density) \n"
			  << "\t a (lattice constant).";
		ifile >> temp;
		ofile << temp << " = ";
		ifile >> value; 
	
		double factor;
		if (strcmp(temp, "D") == 0) {
			double D = this->F_Cell.NearestParticleDistance(0);
			for (size_t i = 1; i < this->N; i++)
				D = std::min(D, this->F_Cell.NearestParticleDistance(i));

			factor = value / D;
		}
		else if (strcmp(temp, "rho") == 0) {
			double currRho = this->F_Cell.NumParticle() / this->F_Cell.PeriodicVolume();
			factor = 1.0/pow(value / currRho, 1.0 / this->d);
		}
		else if (strcmp(temp, "a") == 0) {
			double currA = sqrt(this->F_Cell.GetBasisVector(0).Modulus2());
			for (int i = 1; i < this->d; i++)
				currA = std::max(currA, this->F_Cell.GetBasisVector(i).Modulus2());

			factor = value / currA;
		}
		else {
			ofile << "Please choose a correct standard." << std::endl;
			exit(1);
		}

		this->F_Cell.Rescale(factor);
	}




	ofile << "Number of unit cells in each side = ";	Echo<size_t>(ifile, ofile, NumInEachSide);
	assert(NumInEachSide > 0); 
	this->SetNumCells(NumInEachSide);	
}

/***********************************
	LatticeGen: Member functions
***********************************/

/** Generate a periodic sphere packing.
 * @param[out] c	A periodic SpherePacking */
void LatticeGen::GenerateP(SpherePacking & c) {
	size_t NumCells = 1;
	for (int i = 0; i < this->d; i++)
		NumCells *= this->NumInSide[i];
	
	if (NumCells == 1) {
		//There is only one unit cell
		c = SpherePacking(this->F_Cell);
	}
	else {
		//There are multiple unit cells.
		Configuration cc = MultiplicateStructure(Configuration(this->F_Cell, "a"), this->NumInSide);
		c = SpherePacking(cc, 0);

		//Assign particle radii
		int idx = 0;
		for (size_t i = 0; i < this->N; i++) {
			for ( ; idx % NumCells > 0 || idx == 0; idx++) {
				c.UpdateCharacteristics(this->F_Cell.GetCharacteristics(i), idx);
			}
		}
	}
}

/** Update the maximal packing fraction, but doesn't change particle radii.*/
double LatticeGen::UpdateMaxPhi() {
	double cellLeng = sqrt(this->F_Cell.GetBasisVector(0).Modulus2());
	for (int i = 1; i < this->d; i++)
		cellLeng = std::min(cellLeng, sqrt(this->F_Cell.GetBasisVector(i).Modulus2()));

	double PrtVolume = 0.0;	//Particle volume
	//		NND;		//nearest neighbor distance
	for (size_t i = 0; i < this->N; i++) {
		PrtVolume += HyperSphere_Volume(this->d, 0.5*this->F_Cell.NearestParticleDistance(i));
	}

	this->maxPhi = PrtVolume/this->F_Cell.PeriodicVolume();
	return this->maxPhi;
}

/** Check particle radii to 0.5 of their nearest neighbor distances. */
void LatticeGen::Set2MaxPhi() {
	for (size_t i = 0; i < this->N; i++) {
		this->F_Cell.UpdateCharacteristics(0.5*this->F_Cell.NearestParticleDistance(i), i);
	}
	if (this->F_Cell.CheckOverlap())
		std::cout << "Something wrong\n";
	else {
		this->maxPhi = this->F_Cell.PackingFraction();
		std::cout << "Particle radii are adjusted to 0.5 of nearest neighbor distances.\n"
			<< "New maximal packing fraction is " << this->maxPhi << "\n";
	}
	

}

/** Add a spherical particle in the fundamental cell.
 * If the added particle is overlapping with other particles, particle radii are adjusted.
 * @param relativePosition	Relative coordinates at which a spherical particle is added.
 * @param radius			Radius of an added particle. */
void LatticeGen::AddPrt(const GeometryVector & relativePosition, double radius) {
	this->N++;
	this->F_Cell.Insert(radius, relativePosition);
	this->rho += 1.0 / this->F_Cell.PeriodicVolume();
	this->UpdateMaxPhi();

	if (this->F_Cell.CheckOverlap()) {
		std::cout << "Particles are overlapping.\n";
		this->Set2MaxPhi();
	}
}

/** Remove the i th particle in the fundamental cell. */
void LatticeGen::RemovePrt(size_t i) {
	this->N--;
	this->F_Cell.DeleteParticle(i);
	this->rho -= 1.0 / this->F_Cell.PeriodicVolume();
	this->UpdateMaxPhi();
}

/** Resize particles. */
void LatticeGen::SetPackingFraction(double targetPackingFraction) {
	if (targetPackingFraction <= this->maxPhi) {
		double factor = 1.0 / pow(targetPackingFraction / this->F_Cell.PackingFraction(), 1.0 / this->d);
		this->F_Cell.Rescale_prtOnly(factor);
		if (this->F_Cell.CheckOverlap()) {
			this->Set2MaxPhi();
			factor = 1.0 / pow(targetPackingFraction / this->F_Cell.PackingFraction(), 1.0 / this->d);
			this->F_Cell.Rescale_prtOnly(factor);
		}
	}
	else {
		std::cout << "target packing fraction is too high!\n";
	}
}

/** Resize fundamental cell. */
void LatticeGen::SetRho(double numDensity) {
	double factor = 1.0/pow(numDensity / this->rho, 1.0 / this->d);
	this->rho = numDensity;
	this->F_Cell.Rescale(factor);
}

/** Construct a GenPtsDefect class via an input stream. */
GenPtsDefect::GenPtsDefect(std::istream & ifile, std::ostream & ofile) : ConfigGen() {
	int temp1;
	double temp2;
	ofile << "----------------------------------\n"
		  << "		      Point defects		    \n"
		  << "----------------------------------\n" << std::endl;
	ofile << "The expected concentration of point defects (# of point defects / total # of particles ) = ";
	ifile >> temp2;
	if (temp2 > 0 && temp2 < 1.0) {
		this->MeanPointDefects = temp2;
	}
	else {
		ofile << "This value should be between 0 and 1" << std::endl;
		return;
	}

	ofile << "Type of distribution of the number of defects (0 = fixed/ 1 = Poisson) = ";
	ifile >> temp1;
	if (temp1 == 0 || temp1 == 1)
		this->dist_type = temp1;
	else {
		std::cout << "Please choose either of 0 and 1\n";
	}

	ofile << "Type of point defects (0 =Vacancies/ 1 = interstitial)";
	ifile >> temp1;
	if (temp1 == 0 || temp1 == 1)
		this->DefectType = temp1;
	else {
		std::cout << "Please choose either of 0 and 1\n";
	}
}





#ifdef STRUCTUREFACTOR_INCLUDED
#include "StructureFactor.h"
#endif

size_t numK_points(PeriodicCellList<Empty> & FundamentalCell, std::vector<GeometryVector> & ks) {
#ifdef STRUCTUREFACTOR_INCLUDED

	if (FundamentalCell.PeriodicVolume() > 0) {
		double CircularKMax, LinearKMax;
		std::cout << "Circular Kmax = ";
		std::cin >> CircularKMax;
		std::cout << "Linear Kmax = ";
		std::cin >> LinearKMax;

		ks = GetKs(FundamentalCell, CircularKMax, LinearKMax);
		return ks.size();

	}
	else {
		std::cerr << "FundamentalCell is not valid." << std::endl;
		return 0;
	}
#else

	std::cerr << "Please include \"StructureFactors.h\"\n";
	return 0;
#endif // STRUCTUREFACTOR_INCLUDED

}

