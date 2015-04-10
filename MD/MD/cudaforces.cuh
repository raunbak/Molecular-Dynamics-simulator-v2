#ifndef __cudaforces_H_
#define __cudaforces_H_

//void CoulombWrapper(float * ForceX, float * ForceY, float * ForceZ, const float * PosX,const float * PosY,const float * PosZ, const int Nions);


void CudaCoulombAlloc(float ** ForceX_d, float ** ForceY_d, float ** ForceZ_d, float ** PosX_d, float ** PosY_d, float ** PosZ_d, float ** Charge_d, const int Nions);

void FastCoulombWrapper(float * ForceX, float * ForceY, float * ForceZ,const float * PosX,const float * PosY,const float * PosZ, float * ForceX_d, float * ForceY_d, float * ForceZ_d, float * PosX_d, float * PosY_d, float * PosZ_d,const float * Charge,float * Charge_d, const int Nions);

void CudaCoulombFree(float * ForceX_d, float * ForceY_d, float * ForceZ_d, float * PosX_d, float * PosY_d, float * PosZ_d, float * Charge_d);
#endif
