// forces.h

#ifndef FORCES_H_
#define FORCES_H_


#include "FastIon.h"
#include "FastEnsemble.h"
#include <cstring>
#include "Trap.h"

// Direct calculation forces

double Ftot(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend,double dt, Trap trap);
double Fcoulumb(FastEnsemble & ensemble, int N, int dim);
double Ftrap(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend, double dt, Trap trap);
double Ffriction(FastEnsemble & ensemble, int N, int dim);
double Fpseudo(FastEnsemble & ensemble, int N, int dim, double Vrf, double Vend, Trap trap);




#endif
