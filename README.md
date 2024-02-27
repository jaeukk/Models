# Summary 
CLI of generating some selective disordered models.

# Install


# Available models:
	* Uncorrelated point patterns, i.e., Poisson (or Binomial) Point Patterns
	* (Saturated) Random Sequential Addition (RSA) Packings of Spheres
	// TODO: equilibrium hard spheres, lattices, ...

# Usages
	./models.out -h 	# show a list of available models

	* RSA:
	./models.out models	<<< "RSA ${dimensions} ${particle_number} ${number_density} ${packing_fraction} ${seed} "
		# the code generates the saturated RSA, if ${packing_fraction} is greater than or equal to the known saturated packing fraction (e.g., 1.0)
		