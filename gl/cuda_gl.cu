#include "cuda_gl.cuh"
#include <iostream>
#include <algorithm>

namespace custom_gl
{
    __global__ void map_to_gpu(float * output, float* input, int size)
    {
        int numThreads = blockDim.x * gridDim.x;
        int global_id = threadIdx.x + blockIdx.x * blockDim.x;

        for (int id = global_id; id < size; id+=numThreads)
        {
            output[4*(id)+3] = input[id];
        }
        return;
    }

    cuda_gl::cuda_gl(unsigned int * VBO, int size)
    {
        this->size = size;
        block = 256;
        grid = (size + block.x-1)/block.x;
        cudaMallocManaged(&data_device, sizeof(float)*size);
        cudaGraphicsGLRegisterBuffer(&vbo_res, *VBO,
                                 cudaGraphicsMapFlagsNone);
        cudaGraphicsMapResources(1, &vbo_res, 0);
        cudaGraphicsResourceGetMappedPointer((void**)&gl_data,
                   nullptr, vbo_res);
    }

    cuda_gl::~cuda_gl()
    {
        cudaGraphicsUnmapResources(1, &vbo_res, 0);
        cudaFree(data_device);
    }

    void cuda_gl::update_gpu_data(float * data)
    {
        //auto max_val = *(std::max_element(data, data + size));
        //if (max_val < 1e-7)
        //    max_val = 1.0;
        //std::cout << max_val;
        //std::copy(data, data + size, data_device);
        #pragma omp parallel for
        for (int i = 0; i < size; i++)
            data_device[i] = data[i];//*100;//max_val;
        map_to_gpu<<<grid, block>>>(gl_data, data_device, size);
        cudaDeviceSynchronize();
    }
}
