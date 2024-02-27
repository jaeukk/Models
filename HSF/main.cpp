/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: February 2022 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

/** \file main.cpp 
 * \brief Simple CLI for generating equilibrium hard sphere fluids in 1,2,3 dimensions. */

#if defined(DELLA)
std::string masterDir("/tigress/jaeukk");
#else
std::string masterDir(".");
#endif

#include "../GeometryVector.h"
#include "../etc.h"
#include "../RandomGenerator.h"
#include "../PeriodicCellList.h"
#include "../StructureFactor.h"

//#include "../RandomSequentialAddition.h"
#include "../GenerateConfigs.h"
#include "../HardSphereFluids.h"
#include "HSF_Gen.h"

int main(int argc, char ** argv){
	char tempstring[1000];
//	char name[300] = {};
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	RandomGenerator rngGod(0), rng(0);

	/* CLI for equilibrum hard spheres. */
	ofile << "binomial? (y/n)\n";
	ifile >> tempstring;
	bool isBinary = (tempstring[0] == 'y');

	if(isBinary){
        ofile << "A simple code for generating binary mixture of hard sphere fluids with a given packing frations in canonical ensemble\n";
		std::string dir = masterDir + std::string("/results/HS_Fluids/");
		int dimension = 0, N1 = 0, N2 = 0, numConfig = 0, samInter = 0, numThreads=1, startIdx=0;
		double phi = 0.0, rho = 1.0;
		double sigma1_to_sigma2 = 1.0; /* diameter ratio */
		double sigma1 = 0.0, sigma2 = 0.0; // particle diameter
		/* Input parameters */
		{
			ofile << "dimension = ";                                ifile >> dimension;     ofile << dimension <<"\n";

			ofile << "the number of particles 1 = ";  ifile >> N1;             ofile << N1 << "\n";
			ofile << "the number of particles 2 = ";  ifile >> N2;             ofile << N2 << "\n";

			ofile << "sigma1 / sigma2 = ";	ifile >> sigma1_to_sigma2;	ofile << sigma1_to_sigma2 <<"\n";

			ofile << "the number density = ";               ifile >> rho;   ofile << rho << "\n";

			ofile << "packing fraction = ";                 ifile >> phi;   ofile << phi << "\n";

			ofile << "sampling interval (MC  cycles) = "; 	ifile >> samInter; ofile << samInter << "\n";

			ofile << "number of configurations = "; ifile >> numConfig; ofile << numConfig << "\n";

			ofile << "numThreads = ";       ifile >> numThreads; ofile << numThreads << "\n";

			ofile << "start idx = ";        ifile >> startIdx; ofile << startIdx <<"\n";
		}
		double phi_temp = phi * (N1+N2)/(N1+N2/pow(sigma1_to_sigma2, dimension));
		int iterPerSeed = numConfig / numThreads;
		
		int start = (int)floor(startIdx / (double) numThreads + 0.1);
#pragma omp parallel for num_threads(numThreads)
		for (int i = 0; i < numThreads; i++) {
			HSF_Gen Seed(dimension, N1+N2, phi_temp, rho);
			/* resize smaller particles */
			double R2 = Seed.setting.pSpheres->GetCharacteristics(0)/sigma1_to_sigma2;
			for (int j=N1; j<N1+N2; j++)
				Seed.setting.pSpheres->UpdateCharacteristics(R2, j);
			Seed.SetInterval(samInter);
			Seed.setting.SetRandomSeed(i+std::time(NULL));
			Seed.setting.Switch2Polydisperse();
			SpherePacking c;

			char filename[300];
			sprintf(filename, "%s%dD_HSF_%0.4f_N-%d-%d_sigmas-%0.2f_%04d", dir.c_str(), dimension, phi, N1, N2, sigma1_to_sigma2, i);
			if(start > 0){
					ReadConfiguration(c, filename);
					Seed.setting.ReplaceConfig(c);
			}

			for (int j = start; j < iterPerSeed; j++) {
				Seed.GenerateP(c);
				sprintf(filename, "%s%dD_HSF_%0.4f_N-%d-%d_sigmas-%0.2f_%04d", dir.c_str(), dimension, phi, N1, N2, sigma1_to_sigma2, i + j * numThreads);
				if(c.CheckOverlap())
				#pragma omp critical
				{
						ofile << i + j * numThreads <<"th config has overlapping particles).\n";
				}
				else{
					/* entire */
					WriteConfiguration(c, filename);
					/* particle 1 */
					sprintf(tempstring, "%s-prt1", filename);
					Configuration c_temp(c, "a");
					for (int k=0; k<N2; k++)
						c_temp.DeleteParticle(N1);
					WriteConfiguration(c_temp, tempstring); 
					/* particle 2 */
					sprintf(tempstring, "%s-prt2", filename);
					c_temp = Configuration(c, "a");
					for (int k=0; k<N1; k++)
						c_temp.DeleteParticle(0);
					WriteConfiguration(c_temp, tempstring);
#pragma omp critical
					{
						ofile << i + j * numThreads <<std::endl;
					}
				}
			}
		}

		// auto getC = [&dir, &dimension, &phi, &N](size_t i)->Configuration {
		// 		char filename[300] = "";
		// 		sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04zu", dir.c_str(), dimension, phi, N, i);
		// 		SpherePacking c;
		// 		ReadConfiguration(c, filename);
		// 		Configuration result(c,"a");
		// 		return result;
		// };
		// omp_set_num_threads(numThreads);
		// double K = 4.0*pi / pow(N, 1.0 / (double)dimension);
		// std::vector<GeometryVector> result;
		// IsotropicStructureFactor(getC, numConfig, 3*K, 3*K, result, K);
		// WriteFunction(result, (std::string("temp_")+std::to_string(phi)).c_str());
		// ofile << "S(0): simulation \n";
		// ofile << result[0];
		// ofile << "S(0): theory \n";
		// if(dimension==3)
		// 		ofile << pow(1.0 - phi, 3) / (1.0 + phi + 0.34406 * phi*phi - 0.12802 * phi*phi*phi) << std::endl;	
	}
	else
	{	
        ofile << "Simple code for generating hard sphere fluids with a given packing frations in canonical ensemble\n";
		std::string dir = masterDir + std::string("/results/HS_Fluids/");
		int dimension = 0, N = 0, numConfig = 0, samInter = 0, numThreads=1, startIdx=0;
		double phi = 0.0, rho = 1.0;

		/* Input parameters */
		{
			ofile << "dimension = ";                                ifile >> dimension;     ofile << dimension <<"\n";

			ofile << "the number of particles = ";  ifile >> N;             ofile << N << "\n";

			ofile << "the number density = ";               ifile >> rho;   ofile << rho << "\n";

			ofile << "packing fraction = ";                 ifile >> phi;   ofile << phi << "\n";

			ofile << "sampling interval (MC  cycles) = "; ifile >> samInter; ofile << samInter << "\n";

			ofile << "number of configurations = "; ifile >> numConfig; ofile << numConfig << "\n";

			ofile << "numThreads = ";       ifile >> numThreads; ofile << numThreads << "\n";

			ofile << "start idx = ";        ifile >> startIdx; ofile << startIdx <<"\n";
		}

		int iterPerSeed = numConfig / numThreads;
		
		int start = (int)floor(startIdx / (double) numThreads + 0.1);
#pragma omp parallel for num_threads(numThreads)
		for (int i = 0; i < numThreads; i++) {
				HSF_Gen Seed(dimension, N, phi, rho);
				Seed.SetInterval(samInter);
				Seed.setting.SetRandomSeed(i+std::time(NULL));
				SpherePacking c;

				char filename[300];
				sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i);
				if(start > 0){
						ReadConfiguration(c, filename);
						Seed.setting.ReplaceConfig(c);
				}

				for (int j = start; j < iterPerSeed; j++) {
						Seed.GenerateP(c);
						sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04d", dir.c_str(), dimension, phi, N, i + j * numThreads);
						if(c.CheckOverlap())
						#pragma omp critical
						{
								ofile << i + j * numThreads <<"th config has overlapping particles).\n";
						}
						else{
							WriteConfiguration(c, filename);
		#pragma omp critical
						{
								ofile << i + j * numThreads <<std::endl;
						}
						}
				}
		}

		auto getC = [&dir, &dimension, &phi, &N](size_t i)->Configuration {
				char filename[300] = "";
				sprintf(filename, "%s%dD_HSF_%0.4f_%d_init_config_%04zu", dir.c_str(), dimension, phi, N, i);
				SpherePacking c;
				ReadConfiguration(c, filename);
				Configuration result(c,"a");
				return result;
		};
		omp_set_num_threads(numThreads);
		double K = 4.0*pi / pow(N, 1.0 / (double)dimension);
		std::vector<GeometryVector> result;
		IsotropicStructureFactor(getC, numConfig, 3*K, 3*K, result, K);
		WriteFunction(result, (std::string("temp_")+std::to_string(phi)).c_str());
		ofile << "S(0): simulation \n";
		ofile << result[0];
		ofile << "S(0): theory \n";
		if(dimension==3)
				ofile << pow(1.0 - phi, 3) / (1.0 + phi + 0.34406 * phi*phi - 0.12802 * phi*phi*phi) << std::endl;	


	}
	return 0;
}