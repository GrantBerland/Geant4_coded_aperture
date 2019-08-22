## Geant4 Shielding Simulation 

For determining background radiation-inducing noise counts on an spacecraft-based X-ray detector.

This branch uses native Geant geometry.

Mission parameters:
* 500 km, near polar orbit
* loss cone angle of approximately 64 degrees during radiation belt pass-thrus
* 10^5 e-/cm^2/s flux event

* In numberOfParticles.txt the 4 lines describe:
	* Multiplicative factor to uniformly modify number of particles generated
	* 0: auto-calculated, non-0: number of trapped particles to generate
	* 0: auto-calculated, non-0: number of loss cone particles to generate
	* 0: auto-calculated, non-0: number of signal photons to generate
	* 0: debug flag off, 1: debug flag on, prints autocalculated particle number and exits
