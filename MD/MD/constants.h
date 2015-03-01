
#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cmath>
#include <cstring>

// math constants
const double PI = acos(-1.0);
// math on cuda
const float fPI = 3.1415926535897932384626;

// physical constants
const double c0 = 299792458; // speed of light
const double eps0 = 8.854187817e-12; // Vacuum permittivity
const double e = 1.602176487e-19; // electron charge in C
const double hbar = 1.054571628e-34; // hbar
const double Kb = 1.380650424e-23; // Boltzman konstant
const double u2kg = 1.66053878283e-27; // unit to kg
// physicsal constants as float - for cuda
const float fc0 = 299792458; // speed of light
const float feps0 = 8.854187817e-12; // Vacuum permittivity
const float fe = 1.602176487e-19; // electron charge in C
const float fhbar = 1.054571628e-34; // hbar
const float fKb = 1.380650424e-23; // Boltzman konstant
const float fu2kg = 1.66053878283e-27; // unit to kg

// trap constants
//const double OmegaRF = 2*PI*4.01e6; // rf freq. in Hz
//const double r0 = 0.00235;	// electrode inscribed radius in m
//const double z0 = 0.0025;  	// center electrode length in m
//const double eta = 0.342; 	// Axial geometrical constant



// Grid spacing in String Crystal generator - outdated
//const double GridSpacing = 20e-6;
//const double Tinitial = 0.010*0; // initial temperature, Skal den være 0 ?

// friction coefficient (laser-cooling)
//const double beta=1.430e-20;
const double beta=2e-22; //The real factor

// random velocity kick
//const double vkick = 0*2.68e-6;

// Default temperature - when not specified by SetSteadyStateTemperature or not using friction force!
//const double SteadyStateTzSec = 0.01; // in Kelvin

// timesteps pr rf-cycle
//const double StepsPrPeriode = 105;//105;  // 105 er fint. evt 25 45?

// time step
//const double dt = 1/(OmegaRF/2/PI)/StepsPrPeriode;
//const double dt = 1/(OmegaRF/2/PI)/StepsPrPeriode; //psedoTime

// Old bin size : 200 - 200 - 520 

// Number of bins in 3d histogram

// For sphere
//const int HistNx = 250;
//const int HistNy = 250; // x and y should be the same...
//const int HistNz = 250;
/*
// For Oval (Normal crystal
// Number of bins in 3d histogram
const int HistNx = 250;
const int HistNy = 250; // x and y should be the same...
const int HistNz = 250;
*/

// Pixel size - that is the length of each bin in the 3d histogram
//const double PixelToDistance = 0.89e-6;

// After this timestep we should start building the histogram
//const int StartRecordingHistogram = 360000;

#endif /* CONSTANTS_H_ */