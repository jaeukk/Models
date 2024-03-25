/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file GenerateConfigs.cpp
 *	\brief	Implementations of functions declared in GenerateConfigs and derived classes. */


#include "GenerateConfigs.h"

ConfigGen::ConfigGen(int dimensions, int N, double rho)
{
	this->d = dimensions;
	this->N = N;
	this->rho = rho;
}

ConfigGen::ConfigGen(std::istream &ifile, std::ostream &ofile)
{
	ofile << "Dimensions = "; 		Echo (ifile, ofile, this->d);
	ofile << "Particle Number = "; 	Echo (ifile, ofile, this->N);
	ofile << "Number Density = "; 	Echo (ifile, ofile, this->rho);
}

///***********************************
//	Implementation of an empty box
// ***********************************/
//
//void EmptyBox::SetBasis(const std::vector<GeometryVector> & basis) {
//	std::vector<GeometryVector> tBasis(basis);
//	this->d = basis.size();
//	this->F_Cell = Configuration(tBasis.size(), &tBasis[0], 1.0);
//	
//	//Check whether basis vectors are linearly independent or not.
//	if (this->F_Cell.PeriodicVolume() > 0) {
//		this->Cell_Volume = this->F_Cell.PeriodicVolume();
//	}
//	else {
//		std::cerr << "Invalid basis vectors!\n";
//		this->F_Cell = Configuration();
//	}
//}
//
//EmptyBox::EmptyBox() {
//	/* inherited parameters */
//	this->d = 0;	this->rho = 0;	this->N = 0;	this->numThreads = 1;
//	/* non-inherited ones */
//	this->Cell_Volume = 0;	this->F_Cell = Configuration();
//}
//
//EmptyBox::EmptyBox(int dimension, const std::vector<GeometryVector> & BasisVectors) : EmptyBox(){
//	assert(dimension == BasisVectors.size());
//	this->SetBasis(BasisVectors);
////	else
////		std::cerr << "inconsistent: the prescribed dimension and size of basis vectors\n";
//}
//
//EmptyBox::EmptyBox(const std::string & cell_name, double cell_volume) : EmptyBox() {
//	std::vector<GeometryVector> basis;
//	GetPredefinedBasis(cell_name, basis);
//	if (basis.size() != 0) {
//		this->SetBasis(basis);
//		this->SetCellVolume(cell_volume);
//	}
//	else {
//		std::cout << "Wrong Cell Type " << cell_name << std::endl;
//	}
//}
//
//EmptyBox::EmptyBox(std::istream & ifile, std::ostream & ofile) : EmptyBox() {
//	std::string cell_name; double cell_volume; std::vector<GeometryVector> basis;
//	ofile << "Dimension = ";	Echo<int>(ifile, ofile, this->d);
//	ofile << "Cell Type = ";	Echo<std::string>(ifile, ofile, cell_name);
//	ofile << "Cell Volume = ";	Echo<double>(ifile, ofile, cell_volume);
//	
//	GetPredefinedBasis(cell_name, basis);
//	if (this->d == basis.size()) 
//	{	this->SetBasis(basis);	this->SetCellVolume(cell_volume);}
//	else
//		std::cerr << "wrong Cell Type " << cell_name <<"\n";
//}
//
//EmptyBox::EmptyBox(const EmptyBox & copy) :EmptyBox() {
//	this->d = copy.d;
//	this->F_Cell = Configuration(copy.F_Cell);
//	this->Cell_Volume = copy.Cell_Volume;
//}

/* -----------------------------
	Poisson point process
 -------------------------------*/


PoissonPtGen::PoissonPtGen(std::istream & ifile, std::ostream & ofile) : ConfigGen(ifile, ofile)
 {
	// ofile	<< "------------------------------------------\n"
	// 		<< "	Spatially uncorrelated point patterns \n"
	// 		<< "------------------------------------------\n" << std::endl;
	ofile << "Modes (0 = cannonical/ 1 = grand cannonical)= ";
	Echo(ifile, ofile, this->mode);
	if (this->mode == 1 || this->mode == 0)
		ofile << this->mode << "\n";
	else {
		std::cout << "Wrong mode!\n";
		return;
	}
}

void PoissonPtGen::GenerateC(Configuration &c)
{
	c = GetUnitCubicBox(this->d, 0.1);
	double volume = this-> N / this->rho;
	c.Resize(volume);

	size_t N_ = this->GetExpectedNumPrts((double)this->N);
	for (size_t i = 0; i < N_; i++) {
		//Randomly place a single point inside the fundamental cell.
		c.Insert("a", this->rng);	
	}
}

void PoissonPtGen::GenerateP(SpherePacking &c)
{
	Configuration c_;
	this->GenerateC(c_);
	c = SpherePacking(c_,0.0);
//	std::cout << "PoissonPtGen::GenerateP is not defined.\n"; 
}
