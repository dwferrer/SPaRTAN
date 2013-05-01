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
        return -accel;
        }

__global__ void
calculate_forces(void *devX, void *devA,long long int N)
{
        float4 *shPosition = SharedMemory<float4>();
        float4 *globalX = (float4 *)devX;
        float4 *globalA = (float4 *)devA;
        float4 myPosition;
        int i, tile;
        float4 acc;
        int gtid = blockIdx.x * blockDim.x + threadIdx.x;
        myPosition = globalX[gtid];
        acc.x = globalA[gtid].x; acc.y = globalA[gtid].y; acc.z = globalA[gtid].z; acc.w = globalA[gtid].w;
        for (i = 0, tile = 0; i < N; i += NThreads, tile++) {
                int idx = tile * blockDim.x + threadIdx.x;
                shPosition[threadIdx.x] = globalX[idx];
                __syncthreads();
                acc = tile_calculation(myPosition, acc);
                __syncthreads();
        }
        // Save the result in global memory for the integration step.
        float4 acc4 = {acc.x, acc.y, acc.z, 0.0f};
        globalA[gtid] = acc4;
}


void gpugravity(float * pos, float *accel, long long int N){
        float4 *positions = (float4 *) pos;
        float4 *acc = (float4 *) accel;
        int size = N*sizeof(float4);

        float4 * d_pos, *d_acc;
        int d_size = N*sizeof(float4);
        cudaMalloc((void **) &d_pos,d_size);
        cudaMalloc((void **) &d_acc,d_size);


        cudaMemcpy(d_pos,positions,size,cudaMemcpyHostToDevice);
        cudaMemcpy(d_acc,positions,size,cudaMemcpyHostToDevice);

        calculate_forces<<<(N+NThreads-1)/NThreads,NThreads,NThreads*sizeof(float4)>>>(d_pos,d_acc,N);
        cudaMemcpy(acc,d_acc,size,cudaMemcpyDeviceToHost);
}
