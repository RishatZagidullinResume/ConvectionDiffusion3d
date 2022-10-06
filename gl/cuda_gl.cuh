#pragma once
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

namespace custom_gl
{
    class cuda_gl
    {
    private:
        cudaGraphicsResource *vbo_res;
        float * gl_data;
        float * data_device;
        int size;
        dim3 block;
        dim3 grid;
    public:
        void update_gpu_data(float * data);
        cuda_gl(unsigned int * VBO, int size);
        ~cuda_gl();
    };
}
