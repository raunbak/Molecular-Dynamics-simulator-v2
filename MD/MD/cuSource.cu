#include "stdafx.h" // Pre-compiled header.

#include <stdio.h>
#include <cuda.h>
#include "cudaforces.cuh"
//#include "constants.h"


// Kernel that executes on the CUDA device
__global__ void CoulombForce(float *ForceX, float *ForceY, float *ForceZ, const float *PosX,const float *PosY,const float *PosZ,const float *Charge, int N)
{

	// Constants
	 float fPI = acos(-1.0);
	 float feps0 = 8.854187817e-12; // Vacuum permittivity
	 float fe = 1.602176487e-19; // electron charge in C



	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	ForceX[idx] = 0.0;
	ForceY[idx] = 0.0;
	ForceZ[idx] = 0.0;
	// do the first "half"
	for(int n = 0; n  < idx; n++)
	{
		float dcubed = pow(pow(PosX[n]-PosX[idx], 2)+pow(PosY[n]-PosY[idx], 2)+pow(PosZ[n]-PosZ[idx], 2),((float) 1.5));
		ForceX[idx] += (PosX[idx]-PosX[n])/dcubed  *(fe*Charge[idx])*(fe*Charge[n])/(4*fPI*feps0);
		ForceY[idx] += (PosY[idx]-PosY[n])/dcubed  *(fe*Charge[idx])*(fe*Charge[n])/(4*fPI*feps0);
		ForceZ[idx] += (PosZ[idx]-PosZ[n])/dcubed  *(fe*Charge[idx])*(fe*Charge[n])/(4*fPI*feps0);

	}

	// do the rest
	for(int n = idx+1; n < N; n++)
	{
		float dcubed = pow(pow(PosX[n]-PosX[idx], 2)+pow(PosY[n]-PosY[idx], 2)+pow(PosZ[n]-PosZ[idx], 2),((float) 1.5));
		ForceX[idx] += (PosX[idx]-PosX[n])/dcubed  *(fe*Charge[idx])*(fe*Charge[n])/(4*fPI*feps0);
		ForceY[idx] += (PosY[idx]-PosY[n])/dcubed  *(fe*Charge[idx])*(fe*Charge[n])/(4*fPI*feps0);
		ForceZ[idx] += (PosZ[idx]-PosZ[n])/dcubed  *(fe*Charge[idx])*(fe*Charge[n])/(4*fPI*feps0);

	}

	// fixing the units
	//ForceX[idx] = ForceX[idx]*fe*fe/(4*fPI*feps0);
	//ForceY[idx] = ForceY[idx]*fe*fe/(4*fPI*feps0);
	//ForceZ[idx] = ForceZ[idx]*fe*fe/(4*fPI*feps0);
}

/*
void CoulombWrapper(float * ForceX, float * ForceY, float * ForceZ, const float * PosX,const float * PosY,const float * PosZ,const int Nions)
{
	//printf("Hello World\nThis is the CU file\n");


	float *ForceX_d, *ForceY_d, *ForceZ_d, *PosX_d, *PosY_d, *PosZ_d; // pointers to device
	size_t size = Nions * sizeof(float);
	cudaMalloc((void **) &ForceX_d, size);   // Allocate array on device
	cudaMalloc((void **) &ForceY_d, size);   // Allocate array on device
	cudaMalloc((void **) &ForceZ_d, size);   // Allocate array on device
	cudaMalloc((void **) &PosX_d, size);   // Allocate array on device
	cudaMalloc((void **) &PosY_d, size);   // Allocate array on device
	cudaMalloc((void **) &PosZ_d, size);   // Allocate array on device

	// copying data to device
	cudaMemcpy(PosX_d, PosX, size, cudaMemcpyHostToDevice);
	cudaMemcpy(PosY_d, PosY, size, cudaMemcpyHostToDevice);
	cudaMemcpy(PosZ_d, PosZ, size, cudaMemcpyHostToDevice);

	// Do calculation on device:
	int block_size = 4;
	int n_blocks = Nions/block_size + (Nions%block_size == 0 ? 0:1);
	CoulombForce <<< n_blocks, block_size >>> (ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d, Nions);
	// Retrieve result from device and store it in host array
	cudaMemcpy(ForceX, ForceX_d, sizeof(float)*Nions, cudaMemcpyDeviceToHost);
	cudaMemcpy(ForceY, ForceY_d, sizeof(float)*Nions, cudaMemcpyDeviceToHost);
	cudaMemcpy(ForceZ, ForceZ_d, sizeof(float)*Nions, cudaMemcpyDeviceToHost);
	// Print results
	//for (int i=0; i<Nions; i++)
	//	  printf("%d %f\t %f\t %f\t\n", i, ForceX[i], ForceY[i], ForceZ[i]);
	// Cleanup
	cudaFree(PosX_d);
	cudaFree(PosY_d);
	cudaFree(PosZ_d);
	cudaFree(ForceX_d);
	cudaFree(ForceY_d);
	cudaFree(ForceZ_d);
}
*/
void FastCoulombWrapper(float * ForceX, float * ForceY, float * ForceZ,const float * PosX,const float * PosY,const float * PosZ, float * ForceX_d, float * ForceY_d, float * ForceZ_d, float * PosX_d, float * PosY_d, float * PosZ_d,const float * Charge,float * Charge_d, const int Nions)
{
	size_t size = Nions * sizeof(float);
	// copying data to device
	cudaMemcpy(PosX_d, PosX, size, cudaMemcpyHostToDevice);
	cudaMemcpy(PosY_d, PosY, size, cudaMemcpyHostToDevice);
	cudaMemcpy(PosZ_d, PosZ, size, cudaMemcpyHostToDevice);

	cudaMemcpy(Charge_d, Charge, size, cudaMemcpyHostToDevice);

	// Do calculation on device:
	int block_size = 32;
	int n_blocks = Nions/block_size + (Nions%block_size == 0 ? 0:1);
	CoulombForce <<< n_blocks, block_size >>> (ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d,Charge_d, Nions);
	// Retrieve result from device and store it in host array
	cudaMemcpy(ForceX, ForceX_d, sizeof(float)*Nions, cudaMemcpyDeviceToHost);
	cudaMemcpy(ForceY, ForceY_d, sizeof(float)*Nions, cudaMemcpyDeviceToHost);
	cudaMemcpy(ForceZ, ForceZ_d, sizeof(float)*Nions, cudaMemcpyDeviceToHost);

	//for (int i=0; i<Nions; i++)
		//	  printf("%d %f\t %f\t %f\t\n", i, ForceX[i], ForceY[i], ForceZ[i]);

}

void CudaCoulombAlloc(float ** ForceX_d, float ** ForceY_d, float ** ForceZ_d, float ** PosX_d, float ** PosY_d, float ** PosZ_d, float ** Charge_d, const int Nions)
{


	size_t size = Nions * sizeof(float);
	cudaMalloc((void **) ForceX_d, size);   // Allocate array on device
	cudaMalloc((void **) ForceY_d, size);   // Allocate array on device
	cudaMalloc((void **) ForceZ_d, size);   // Allocate array on device
	cudaMalloc((void **) PosX_d, size);   // Allocate array on device
	cudaMalloc((void **) PosY_d, size);   // Allocate array on device
	cudaMalloc((void **) PosZ_d, size);   // Allocate array on device

	cudaMalloc((void **) Charge_d, size);   // Allocate array on device
}

void CudaCoulombFree(float * ForceX_d, float * ForceY_d, float * ForceZ_d, float * PosX_d, float * PosY_d, float * PosZ_d, float * Charge_d)
{
	cudaFree(PosX_d);
	cudaFree(PosY_d);
	cudaFree(PosZ_d);
	cudaFree(ForceX_d);
	cudaFree(ForceY_d);
	cudaFree(ForceZ_d);

	cudaFree(Charge_d);

	// Testing this magic clean up line...
	cudaDeviceReset();
}



