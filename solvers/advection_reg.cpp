#include "advection_reg.h"

namespace solvers
{
    template<typename T>
    std::pair<T, int> Advection3d_reg<T>::find_inds(T src_pos)
    {
        int ind = -1;
        T distance = 0.0;
        for (int i = 0; i < N*M*K; i++)
        {
            if (src_pos <= 0.0 || src_pos >= (N*M*K-1)*dd)
            {
                ind = 0;
                distance = -1.0;
                break;
            }
            else if (src_pos < i*dd)
            {
                ind = i;//between i and i-1
                distance = i*dd - src_pos;
                break;
            }
            else if (src_pos == i*dd)
            {
                ind = i;
                distance = 0.0;
                break;
            }
        }
        if (ind == -1) std::cout << "advection::find_inds error\n";
        return std::pair<T, int>{distance, ind};
    }

    template<typename T>
    void Advection3d_reg<T>::init_data()
    {
        this->vels = new Vector3d<T> [N*M*K];
        int c = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    //data_new[c] = exp(-5.0*((i*dd-2.0)*(i*dd-2.0)+
                    //                    (j*dd-1.0)*(j*dd-1.0)+
                    //                    (k*dd-1.0)*(k*dd-1.0))
                    //          );
                    vels[c].x = 0.0;
                    vels[c].y = 0.0;
                    T r_sq = pow(j*dd-1., 2)+pow(k*dd-1., 2);
                    if (r_sq<=1.0)
                        vels[c].z = 1.0*(1.-r_sq*r_sq);
                    else
                        vels[c].z = 0.0;
                    c++;
                }
            }
        }
    }

    template<typename T>
    void Advection3d_reg<T>::find_src()
    {
        this->weights = new std::vector<T> [N*M*K];
        this->inds = new std::vector<int> [N*M*K];
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    int idx = k+j*K+i*K*M;
                    Vector3d<T> cur_pos = {k*dd, j*dd, i*dd};
                    Vector3d<T> src_pos = {cur_pos.x-vels[idx].x*dt, 
                                           cur_pos.y-vels[idx].y*dt, 
                                           cur_pos.z-vels[idx].z*dt};
                    auto pair_h = find_inds(src_pos.z);
                    auto pair_w = find_inds(src_pos.y);
                    auto pair_d = find_inds(src_pos.x);
                    //boundary source
                    if (pair_h.first < 0.0)
                    {
                        this->weights[idx].push_back(2.0);
                        this->inds[idx].push_back(
                                        pair_d.second
                                        +pair_w.second*K
                                        +pair_h.second*K*M);
                    }
                    else if (pair_w.first < 0.0 ||
                        pair_d.first < 0.0)
                    {
                        this->weights[idx].push_back(0.0);
                        this->inds[idx].push_back(0);
                    }
                    else if (pair_d.first == 0 &&
                             pair_w.first == 0 &&
                             pair_h.first == 0)
                    {
                        this->weights[idx].push_back(1.0);
                        this->inds[idx].push_back(
                                        pair_d.second
                                        +pair_w.second*K
                                        +pair_h.second*K*M);
                    }
                    else
                    {
                        T total_distance = 0.0;
                        //==================i===j===k==============
                        int cur_ind = pair_d.second+
                                      pair_w.second*K+
                                      pair_h.second*K*M;
                        T cur_distance = pow(
                                pair_d.first*pair_d.first+
                                pair_w.first*pair_w.first+
                                pair_h.first*pair_h.first, 0.5);
                        this->inds[idx].push_back(cur_ind);
                        this->weights[idx].push_back(1./cur_distance);
                        total_distance += 1./cur_distance;
                        //==================i===j===k-1============
                        if (pair_d.first !=0)
                        {
                            cur_ind = (pair_d.second-1)+
                                      pair_w.second*K+
                                      pair_h.second*K*M;
                            cur_distance = pow(
                                 (dd-pair_d.first)*(dd-pair_d.first)+
                                 pair_w.first*pair_w.first+
                                 pair_h.first*pair_h.first, 0.5);
                            this->inds[idx].push_back(cur_ind);
                            this->weights[idx].push_back(
                                               1./cur_distance);
                            total_distance += 1./cur_distance;
                        }
                        //==================i===j-1=k==============
                        if (pair_w.first !=0)
                        {
                            cur_ind = pair_d.second+
                                      (pair_w.second-1)*K+
                                      pair_h.second*K*M;
                            cur_distance = pow(
                                  pair_d.first*pair_d.first+
                                  (dd-pair_w.first)*(dd-pair_w.first)+
                                  pair_h.first*pair_h.first, 0.5);
                            this->inds[idx].push_back(cur_ind);
                            this->weights[idx].push_back(
                                               1./cur_distance);
                            total_distance += 1./cur_distance;
                        }
                        //==================i-1=j===k==============
                        if (pair_h.first !=0)
                        {
                            cur_ind = pair_d.second+
                                      pair_w.second*K+
                                      (pair_h.second-1)*K*M;
                            cur_distance = pow(
                             pair_d.first*pair_d.first+
                             pair_w.first*pair_w.first+
                             (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                            this->inds[idx].push_back(cur_ind);
                            this->weights[idx].push_back(
                                               1./cur_distance);
                            total_distance += 1./cur_distance;
                        }
                        //==================i===j-1=k-1============
                        if (pair_d.first !=0 && pair_w.first !=0)
                        {
                            cur_ind = (pair_d.second-1)+
                                      (pair_w.second-1)*K+
                                      pair_h.second*K*M;
                            cur_distance = pow(
                                 (dd-pair_d.first)*(dd-pair_d.first)+
                                 (dd-pair_w.first)*(dd-pair_w.first)+
                                 pair_h.first*pair_h.first, 0.5);
                            this->inds[idx].push_back(cur_ind);
                            this->weights[idx].push_back(
                                               1./cur_distance);
                            total_distance += 1./cur_distance;
                        }
                        //==================i-1=j-1=k==============
                        if (pair_w.first !=0 && pair_h.first !=0)
                        {
                            cur_ind = pair_d.second+
                                      (pair_w.second-1)*K
                                      +(pair_h.second-1)*K*M;
                            cur_distance = pow(
                             pair_d.first*pair_d.first+
                             (dd-pair_w.first)*(dd-pair_w.first)+
                             (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                            this->inds[idx].push_back(cur_ind);
                            this->weights[idx].push_back(
                                               1./cur_distance);
                            total_distance += 1./cur_distance;
                        }
                        //==================i-1=j===k-1============
                        if (pair_d.first !=0 && pair_h.first !=0)
                        {
                            cur_ind = (pair_d.second-1)+
                                      pair_w.second*K+
                                      (pair_h.second-1)*K*M;
                            cur_distance = pow(
                             (dd-pair_d.first)*(dd-pair_d.first)+
                             pair_w.first*pair_w.first+
                             (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                            this->inds[idx].push_back(cur_ind);
                            this->weights[idx].push_back(
                                               1./cur_distance);
                            total_distance += 1./cur_distance;
                        }
                        //==================i-1=j-1=k-1============
                        if (pair_d.first !=0 && 
                            pair_w.first !=0 && 
                            pair_h.first !=0)
                        {
                            cur_ind = (pair_d.second-1)+
                                      (pair_w.second-1)*K+
                                      (pair_h.second-1)*K*M;
                            cur_distance = pow(
                             (dd-pair_d.first)*(dd-pair_d.first)+
                             (dd-pair_w.first)*(dd-pair_w.first)+
                             (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                            this->inds[idx].push_back(cur_ind);
                            this->weights[idx].push_back(
                                               1./cur_distance);
                            total_distance += 1./cur_distance;
                        }
                        //=========================================
                        for (int m=0; m<this->weights[idx].size(); m++)
                            this->weights[idx][m] = 
                                  this->weights[idx][m]/total_distance;
                    }
                }
            }
        }
        return;
    }

    template<typename T>
    void Advection3d_reg<T>::iteration(T* data, T* data_new, int coef, bool boundary)
    {
        //std::swap(data_new, data);
        for (int i = 0; i < N*M*K; i++)
        {
            data_new[i] = 0.0;
            for (int j = 0; j < inds[i].size(); j++)
            {
                if (weights[i][j] == 2.0 && boundary)
                    data_new[i] = 1.0;
                else if (weights[i][j] == 2.0)
                    data_new[i] = data[inds[i][j]];
                else
                    data_new[i] += data[inds[i][j]]*weights[i][j];
            }
        }
        return;
    }

    template<typename T>
    Advection3d_reg<T>::~Advection3d_reg()
    {
        delete [] weights;
        delete [] inds;
        delete [] vels;
    }

    template class Advection3d_reg<double>;
    template class Advection3d_reg<float>;
}

