/** 
 *  Author	:	Jaeuk Kim
 *  Email	: phy000.kim@gmail.com
 *  Date	: Febrauray 2022
 */


#ifndef HSF_Gen_H
#define HSF_Gen_H

/* in $(core) */
#include <PeriodicCellList.h>
#include <GeometryVector.h>

#include "HardSphereFluids.h"
#include "../GenerateConfigs.h"

class HSF_Gen : public ConfigGen {
private:
	int samplingInterval;
public:
	HSF setting;

	HSF_Gen() : ConfigGen::ConfigGen() {};

	HSF_Gen(int dimension, int N, double packingFraction, double numberDensity);

	HSF_Gen(std::istream & ifile, std::ostream & ofile);

	virtual void SetNumPrts(int num) {
		std::cout << "SetNumPrts(int) is not supported in HSF\n";
	}
	virtual void SetRho(double rho) {
		this->setting.RescaleConfig(rho);
	}
	virtual void SetPackingFraction(double phi) {
		double factor = pow(phi / this->setting.PackingFraction, 1.0 / this->d);
		SpherePacking temp(*this->setting.pSpheres);
		temp.Rescale_prtOnly(factor);
		this->setting.ReplaceConfig(temp);
		this->setting.Prepare();
	}

	void SetInterval(int sampleinterval) { this->samplingInterval = sampleinterval; }

	virtual void GenerateC(Configuration &c);
	virtual void GenerateP(SpherePacking &c);
};

HSF_Gen::HSF_Gen(int dimension, int N, double packingFraction, double numberDensity) {
	this->d = dimension;
	this->N = N;
	this->rho = numberDensity;

	setting = HSF(dimension, N, packingFraction);	setting.SetNumThreads(1);
	setting.RescaleConfig(numberDensity);
	//setting.Prepare();
	this->samplingInterval = 1000;
	std::cout << "Default sampling interval = " << this->samplingInterval << " MCcycles\n";
}

HSF_Gen::HSF_Gen(std::istream & ifile, std::ostream & ofile) {
	double packingFraction;
	ofile << "dimension = ";				ifile >> this->d;
	ofile << "the number of particles = ";	ifile >> this->N;
	ofile << "the number density = ";		ifile >> this->rho;
	ofile << "packing fraction = ";			ifile >> packingFraction;
	ofile << "sampling interval (MC  cycles) = "; ifile >> this->samplingInterval;

	setting = HSF(this->d, this->N, packingFraction);	setting.SetNumThreads(1);
	setting.RescaleConfig(this->rho);
	//setting.Prepare();
	//this->samplingInterval = 100;
	std::cout << "Default sampling interval = " << this->samplingInterval << " MCcycles\n";
}

void HSF_Gen::GenerateP(SpherePacking &c) {
	if (!this->setting.IsEquilibrated()) {
		this->setting.Prepare();
		this->setting.Evolve(200);
	}
	this->setting.Evolve(this->samplingInterval);
	c = SpherePacking(*this->setting.pSpheres);
}

void HSF_Gen::GenerateC(Configuration &c) {
	SpherePacking temp;
	this->GenerateP(temp);
	c = Configuration(temp, "a");
}

#endif