//integrator.h
#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "FastIon.h"
#include "FastEnsemble.h"
#include <cstring>
#include <iostream>
#include "Trap.h"
#include <algorithm>

void TemperatureRescaleLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, int StartHistograms, double Vrf, double Vend, double StepsPrPeriode, Trap & trap);
void TwoIonTypesTemperatureRescaleLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, int StartHistograms, double Vrf, double Vend, double StepsPrPeriode, Trap & trap);

void NewTemperatureRescaleLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, int StartHistograms, double Vrf, double Vend, double StepsPrPeriode, Trap & trap);

//void MADSDynamicTemperatureLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend); // The best one.
//void TauPeriodeCudaLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, double Vrf, double Vend);
#endif
