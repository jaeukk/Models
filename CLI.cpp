/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2022 */

/** \file GenerateConfigs.cpp
 *	\brief	Implementations of functions declared in GenerateConfigs and derived classes. */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>


#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>

#include "GenerateConfigs.h"
//#ifdef RSA_gen
#include "RSA/RSA_Gen.h"
//#endif


size_t Verbosity = 2;
int main(int argc, char ** argv){
	/* An extern variable defined in ${core}/etc.h */
	std::vector<std::string> models {"Poisson"};
	models.push_back("RSA");

	char tempstring[1000];

	std::istream & ifile=std::cin;
	std::ostream & ofile=std::cout;
	ConfigGen * model;

	if (argc > 1){
		if (strcmp(argv[1], "-h")==0){
			/* output available models */
			ofile << "Available models\n";
			for(auto iter = models.begin(); iter != models.end(); ++iter){
				ofile << *iter << "\t";
			}
			ofile << "\n";
			return 0;
		}
		else if(strcmp(argv[1], "models")==0){
			ofile << "Model name = ";
			ifile >> tempstring;

			if (std::find(models.begin(), models.end(), std::string(tempstring)) != models.end() ){
				/* Input  */
				if (strcmp(tempstring, "RSA") == 0){
					model = new RSA_gen(ifile, ofile);
				}
				// else if (strcmp(tempstring, "HSF") == 0){

				// }
			}
			else{
				ofile << tempstring << " is not in the available model list.\n";
			}

			/* Output ... */
			
			size_t Nc;
			std::string output_name, output_type;

			ofile << "Number of configurations = ";
			Echo(ifile, ofile, Nc);

			ofile << "output prefix = ";
			Echo(ifile, ofile, output_name);

			ofile << "output type (ConfigPack/SpherePacking[.txt]) = ";
			Echo(ifile, ofile, output_type);

			if (strcmp(output_type.c_str(), "ConfigPack") == 0){
				ConfigurationPack save(output_name);

				{
					/* check the particle radius */
					SpherePacking c;
					model->GenerateP(c);
					ofile << "Particle radius = " << c.GetMaxRadius() <<"\n";
					save.AddConfig(Configuration(c, "a"));
				}

				ofile << "output file = "<< output_name  << ".ConfigPakc\n";
				ofile << "Generating Configurations ... \n";
				progress_display pd(Nc);
				pd ++ ;
				for(int i = 1; i < Nc; i++){
					Configuration c;
					model->GenerateC(c);
					save.AddConfig(c);
					pd++;	
				}

			}
			else if (strcmp(output_type.c_str(), "SpherePacking") == 0){
				output_name += std::string("__%05d");
				ofile << "output file = "<< output_name <<"\n";
				progress_display pd(Nc);
				
				for(int i = 0; i < Nc; i++){
					sprintf(tempstring, output_name.c_str(), i);
					SpherePacking c;
					model->GenerateP(c);
					WriteConfiguration(c, tempstring);
					pd++;	
				}
			}
			ofile << "\n Done!\n";



			return 0;
		}
		else{

		}
	}
	

}