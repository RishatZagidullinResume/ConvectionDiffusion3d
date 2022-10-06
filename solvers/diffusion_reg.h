#pragma once
#include <iostream>
#include "../geometry/vector3d.h"

namespace solvers
{
    template <typename T>
    class Diffusion3d_reg
    {
    private:
        //N->z, M->y, K->x
        int N, M, K;
        T dd, dt, D;
        Vector3d<T> * vels;
        void init_data();
        void solve_matrix (int n, T *a, T *c, T *b, T *f, T *x);
    public:
        void iteration(T * data, T * data_new, int coef, bool bound);
        Diffusion3d_reg(T D, int N, int M, int K, T dd, T dt):
                    D(D), N(N), M(M), K(K), dd(dd), dt(dt)
                    {
                        init_data();
                    };
        ~Diffusion3d_reg();
    };
}
