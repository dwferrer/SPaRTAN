#include <stdio.h>
#include <gsl/gsl_rng.h>

//shamelessly stolen from cuda gems nbody code
#define EPS2 0.00000001
#define NThreads 1024
__device__ float4
bodyBodyInteraction(float4 bi, float4 bj, float4 ai)
{
        float3 r;
        // r_ij [3 FLOPS]
        r.x = bj.x - bi.x;
        r.y = bj.y - bi.y;
        r.z = bj.z - bi.z;
        // distSqr = dot(r_ij, r_ij) + EPS^2 [6 FLOPS]
        float distSqr = r.x * r.x + r.y * r.y + r.z * r.z + ai.w;
        // invDistCube =1/distSqr^(3/2) [4 FLOPS (2 mul, 1 sqrt, 1 inv)]
        float distSixth = distSqr * distSqr * distSqr;
        float invDistCube = 1.0f/sqrtf(distSixth);
        // s = m_j * invDistCube [1 FLOP]
        float s = bj.w * invDistCube;
        // a_i = a_i + s * r_ij [6 FLOPS]
        ai.x += r.x * s;
        ai.y += r.y * s;
        ai.z += r.z * s;
        return ai;
}


template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};


__device__ float4
tile_calculation(float4 myPosition, float4 accel){
        long long  int i;
        float4 *shPosition = SharedMemory<float4>();
        #pragma unroll 32
        for (i = 0; i < blockDim.x; i++) {
        accel = bodyBodyInteraction(myPosition, shPosition[i], accel);
        }
        return accel;
        }

__global__ void
calculate_forces(void *devXsource, void * devXsink, void *devA, int Nsource, int Nsink, int numdevs)
{
        float4 *shPosition = SharedMemory<float4>();
        float4 *globalXsource = (float4 *)devXsource;
        float4 *globalXsink = (float4 *)devXsink;
        float4 *globalA = (float4 *)devA;
        float4 myPosition;
        int i, tile;
        float4 acc;
        int gtid = blockIdx.x * blockDim.x + threadIdx.x;
        if (gtid >= Nsink) return;
        myPosition = globalXsink[gtid];
        acc.x = globalA[gtid].x; acc.y = globalA[gtid].y; acc.z = globalA[gtid].z; acc.w = globalA[gtid].w;
        for (i = 0, tile = 0; i < Nsource; i += NThreads, tile++) {
                int idx = tile * blockDim.x + threadIdx.x;
                shPosition[threadIdx.x] = globalXsource[idx];
                __syncthreads();
                acc = tile_calculation(myPosition, acc);
                __syncthreads();
        }
        // Save the result in global memory for the integration step.
        float4 acc4 = {acc.x, acc.y, acc.z, acc.w};
        globalA[gtid] = acc4;
}


#include <cassert>
#include <stdio.h>

void gpugravity(float * pos, float *accel, long long int N){
        float4 *positions = (float4 *) pos;
        float4 *acc = (float4 *) accel;
        
        int numdevs = 0;
        cudaGetDeviceCount(&numdevs);
        cudaStream_t * streams = new cudaStream_t[numdevs];
        cudaEvent_t * events = new cudaEvent_t[numdevs];
        
        int * devicesinks = new int[numdevs];
        size_t *offset = new size_t[numdevs];
        size_t total_offset = 0;
        int remainingsinks = N;
        int allotment = N/numdevs;
        int d_sourcesize = N*sizeof(float4);
        
        float4 ** d_pos = new float4 *[numdevs];
        float4 ** d_acc = new float4 *[numdevs];
        
        for(int i = 0; i < numdevs; i++){
        	
        	//create the streams and events
        	cudaStreamCreate(&streams[i]);
        	cudaEventCreate(&events[i]);
        	
        	//figure out how many sinks to give each device
        	if (remainingsinks > allotment) devicesinks[i] = allotment;
        	else devicesinks[i] = allotment;
        	remainingsinks -= devicesinks[i];
        	
        	//calculate the offset for each device
        	
        	offset[i] = total_offset;
        	total_offset += devicesinks[i];
        	
		printf("Device %d has %d sinks and an offset of %d\n There are %d particles remaining\n\n",i,devicesinks[i],offset[i],remainingsinks);
        	
		cudaSetDevice(i);
        	
        	int d_sinksize = devicesinks[i] *sizeof(float4);
        	cudaMalloc((void **) &d_pos[i],d_sourcesize);
        	cudaMalloc((void **) &d_acc[i],d_sinksize);      	
        }
        
        
        
        for(int i = 0; i < numdevs; i++){
        	cudaSetDevice(i);
        	int d_sinksize = devicesinks[i] *sizeof(float4);
        	cudaMemcpyAsync(d_pos[i], positions, d_sourcesize,cudaMemcpyHostToDevice,streams[i]);
        	cudaMemcpyAsync(d_acc[i], &(acc[offset[i]]), d_sinksize,cudaMemcpyHostToDevice,streams[i]);    	
        }
        
        for(int i = 0; i < numdevs; i++){
        	cudaSetDevice(i);
        	calculate_forces<<<(devicesinks[i]+NThreads-1)/NThreads,NThreads,NThreads*sizeof(float4),streams[i]>>>(d_pos[i],&((d_pos[i])[offset[i]]), d_acc[i] ,N, devicesinks[i],numdevs);
        }
        
        for(int i = 0; i < numdevs; i++){
        	cudaSetDevice(i);
        	int d_sinksize = devicesinks[i] *sizeof(float4);
        	cudaMemcpyAsync(&(acc[offset[i]]),d_acc[i],d_sinksize,cudaMemcpyDeviceToHost,streams[i]);
        	cudaEventRecord(events[i],streams[i]);        	
        }
        
        //wait for all devices to complete
        for(int i = 0; i < numdevs; i++) cudaEventSynchronize(events[i]);    
	for(int i = 0; i < numdevs; i++){
		cudaSetDevice(i); cudaFree(d_pos[i]); cudaFree(d_acc[i]);
	}
	delete[] offset; delete[] devicesinks;
	
}
