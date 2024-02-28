# Summary 
CLI of generating some selective disordered models.

# Requirements
	* gsl library
	* boost library
		* The `boost/random/*` must be available
		
	* `core` from git@github.com:jaeukk/cores.git    
	* Make

# Compilation
	In the makefile, change the value of $(core) to /path/to/core/from/git .
	
	Then, the code can be compiled with the provided makefile using the standard `make` command.
	This creates an executable ~/EXC/model.out , which is a CLI for handling ConfigPack files.

# Available models:
	* Uncorrelated point patterns, i.e., Poisson (or Binomial) Point Patterns
	* (Saturated) Random Sequential Addition (RSA) Packings of Spheres.
		G. Zhang and S. Torquato, Precise Algorithm to Generate Random Sequential Addition of Hard Hyperspheres at Saturation, Physical Review E, 88, 053312 (2013).
	
	
	
	// TODO: equilibrium hard spheres, lattices, ...

# Usages
	./models.out -h 	# show a list of available models

	* RSA:
	./models.out models	<<< "RSA ${dimensions} ${particle_number} ${number_density} ${packing_fraction} ${seed} ${num_threads} ${configuration_number} ${output_prefix} ${output_type}"
		# the code generates the saturated RSA, if ${packing_fraction} is greater than or equal to the known saturated packing fraction (e.g., 1.0)
		output_type = (ConfigPack/SpherePacking). 
			ConfigPack: save `point configurations' in a ConfigPack file.
			SpherePacking: save `sphere packings' in separated text files.
		