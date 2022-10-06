#pragma once
#include <iostream>
#include <vector>
#include "../geometry/vector3d.h"

namespace solvers
{
    template <typename T>
    class Advection3d_reg
    {
    private:
        //N->z, M->y, K->x
        int N, M, K;
        T dd, dt;
        Vector3d<T> * vels;
        std::vector<T> * weights;
        std::vector<int> * inds;
        std::pair<T, int> find_inds(T src_pos);
        void init_data();
        void find_src();
    public:
        void iteration(T * data, T * data_new, int coef, bool cond);
        Advection3d_reg(int N, int M, int K, T dd, T dt):
                    N(N), M(M), K(K), dd(dd), dt(dt)
                    {
                        init_data();
                        find_src();
                    };
        ~Advection3d_reg();
    };
}
