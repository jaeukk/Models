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
size_t Verbosity = 2;

#include "GenerateConfigs.h"
std::vector<std::string> models({std::string("Poisson")});

//#ifdef RSA_Gen
#include "RSA/RSA_Gen.h"
models.emplace_back("RSA")
//#endif

int main(int argc, char ** argv){
	char tempstring[1000];
	/* An extern variable defined in ${core}/etc.h */
	Verbosity = 3;

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
					model = new RSA_Gen(ifile, ofile);
				}
				// else if (strcmp(tempstring, "HSF") == 0){

				// }
			}
			else{
				ofile << tempstring << " is not in the available model list.\n";
			}

			/* Output ... */


			return 0;
		}
		else{

		}
	}
	

}