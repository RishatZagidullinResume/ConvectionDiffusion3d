#include "advection_tet.h"

double W(double v)
{
    return v < 0.225 ? 1. - pow(v/0.225, 0.4) : 0;
}

template<typename F>
F sigma(F const &r, F const &dr, int const &PML, int const &size)
{
    if (fabs(r)>dr*(size-PML))
    {
        return pow((fabs(r)-(size-PML)*dr)/(PML*dr),
               2) * 3.0*log(10.0)*10.0/(PML*dr);
    }
    else return 0;
}

template<typename F>
F limiter(F const &r_factor, F const & n_plus, F const &n_minus)
{
    F phi_superbee = std::max(std::min(n_minus*r_factor, (F)1.0),
                              std::min(r_factor, n_plus));
    phi_superbee = std::max(phi_superbee, (F)0.0);
    return phi_superbee;
}

template<typename F>
void find_gradients_eigen(std::vector<std::vector<F>> &res,
                          std::vector<F> const &u,
                          std::vector<int> const &num_dums,
                          Vector3d<F> const * const cell_centroids,
                          std::vector<int> const * const neighbor_ids,
                          std::vector<bool> const &sphere_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < u.size(); i++)
    {
        int s = neighbor_ids[i].size();
        Eigen::VectorXd T(s-num_dums[i]);
        Eigen::MatrixXd d(s-num_dums[i], 3);
        int val = 0;
        for (int j = 0; j < s; j++)
        {
            if ((neighbor_ids[i][j]!=-1)&&(neighbor_ids[i][j]!=i))
            {
                d(val,0) = (cell_centroids[neighbor_ids[i][j]].x
                          - cell_centroids[i].x);
                d(val,1) = (cell_centroids[neighbor_ids[i][j]].y
                          - cell_centroids[i].y);
                d(val,2) = (cell_centroids[neighbor_ids[i][j]].z
                          - cell_centroids[i].z);
                F w = 1.0/pow(d(val,0)*d(val,0)
                             +d(val,1)*d(val,1)
                             +d(val,2)*d(val,2),
                          0.5);
                d(val,0) *= w;
                d(val,1) *= w;
                d(val,2) *= w;
                T[val] = (u[neighbor_ids[i][j]] - u[i])*w;
                if (abs(T[val]) < 1e-4) T[val] = 0.0;
                val++;
            }
        }
        //Eigen::Vector3d output = mats[i].bdcSvd(Eigen::ComputeThinU 
        //                     | Eigen::ComputeThinV).solve(T);
        //Eigen::Vector3d output = mats[i].householderQr().solve(T);
        Eigen::Vector3d output =
            (d.transpose() * d)
            .llt().solve(d.transpose() * T);
        res[0][i] = output(0);
        res[1][i] = output(1);
        res[2][i] = output(2);
    }
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::initialize(std::vector<F> &val, Vector3d<F> const  &source)
{
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        F x_sq = (g.cell_centroids[i].x -source.x) * 
                 (g.cell_centroids[i].x-source.x);
        F y_sq = (g.cell_centroids[i].y - source.y) * 
                 (g.cell_centroids[i].y-source.y);
        F z_sq = (g.cell_centroids[i].z - source.z) * 
                 (g.cell_centroids[i].z-source.z);
        val[i] = 1.0*exp(-4.0*(z_sq+y_sq+x_sq));
    }
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::interpolate_to_faces(
                          std::vector<F> const &centered, 
                          std::vector<std::vector<F>> &faced)
{
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (g.cell_neighbor_ids[i*4+j]!=-1)
            {
                F result_p_x = g.beta1_p[i*4+j] * 
                               centered[g.cell_stencil_ids[i]
                               [g.sorted_ids[i*4+j][g.a_p[i*4+j]]]]
                             + g.beta2_p[i*4+j] * 
                               centered[g.cell_stencil_ids[i]
                               [g.sorted_ids[i*4+j][g.b_p[i*4+j]]]]
                            + g.beta3_p[i*4+j] * 
                               centered[g.cell_stencil_ids[i]
                               [g.sorted_ids[i*4+j][g.c_p[i*4+j]]]];

                if (fabs(g.beta1_p[i*4+j] + g.beta2_p[i*4+j]
                       + g.beta3_p[i*4+j] - 1.0) > 1e-5)
                {
                    result_p_x = 0.0;
                }
                F result_m_x = g.beta1_m[i*4+j] * 
                               centered[g.cell_stencil_ids[i]
                               [g.sorted_ids[i*4+j][g.a_m[i*4+j]]]]
                             + g.beta2_m[i*4+j] * 
                               centered[g.cell_stencil_ids[i]
                               [g.sorted_ids[i*4+j][g.b_m[i*4+j]]]]
                             + g.beta3_m[i*4+j] * 
                               centered[g.cell_stencil_ids[i]
                               [g.sorted_ids[i*4+j][g.c_m[i*4+j]]]];
                if (fabs(g.beta1_m[i*4+j] + g.beta2_m[i*4+j]
                       + g.beta3_m[i*4+j] - 1.0) > 1e-5)
                {
                    result_m_x = 0.0;
                }
                F slope_plus_x = (result_p_x - centered[i])
                                 /g.dist_inter_p[i*4+j];
                F slope_minus_x = (centered[i] - result_m_x)
                                 /g.dist_inter_m[i*4+j];
                F n_plus = g.dist_inter_p[i*4+j]
                                 /g.dist_center_face[i*4+j];
                F n_minus = g.dist_inter_m[i*4+j]
                                 /g.dist_center_face[i*4+j];
                F small_cor = 1e-8;
                faced[i][j] = -(centered[i] + slope_plus_x * 
                         limiter((slope_minus_x+small_cor)/
                                 (slope_plus_x+small_cor),
                                 n_plus, n_minus) * 
                                 g.dist_center_face[i*4+j]);
                if (result_m_x == 0.0 || result_p_x == 0.0)
                    faced[i][j] = 0.0;
                if (std::isnan(faced[i][j]))
                {
                    std::cout << i << " is nan" << std::endl;
                    std::cout << g.beta1_m[i*4+j] << " ";
                    std::cout << g.beta2_m[i*4+j] << " ";
                    std::cout << g.beta3_m[i*4+j] << std::endl;
                    std::cout << g.beta1_p[i*4+j] << " ";
                    std::cout << g.beta2_p[i*4+j] << " ";
                    std::cout << g.beta3_p[i*4+j] << std::endl;
                }
            }
            else
                faced[i][j] = 0.0;
        }
    }
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::iteration_advection(
              std::vector<std::vector<F>> & vel, std::vector<F> & val)
{                       
    for (int j = 0; j < 3; j++)
        interpolate_to_faces(vel[j], temp_face_vec[j]);
    interpolate_to_faces(val, temp_face);  
    iteration_muscl(val, temp_face, temp_face_vec);
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::iteration_diffusion(std::vector<F> & val)
{
    find_gradients_eigen(temp_x, val, num_dums, g.cell_centroids, 
                         g.cell_stencil_ids, g.sphere_cells);
    iteration_advection(temp_x, val);
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::poisson_solver(
                std::vector<F> & val, std::vector<F> & rhs)
{
    Eigen::VectorXd T(size);
    #pragma omp parallel for
    for (int i = 0; i < size ; i++)
    {
        if (!g.sphere_cells[i])
        {
            T[i] = -rhs[i];
            //if (denominators[i] == 0.0) T[i] = rhs[i];
            //else T[i] = rhs[i]/denominators[i];
        }
        else T[i] = 0.0;
    }
    Eigen::LeastSquaresConjugateGradient<SpMat> cg;
    cg.compute(mat);
    cg.setTolerance(0.01);
    cg.setMaxIterations(500);
    std::cout << "starting the solver...\n";
    Eigen::VectorXd output = cg.solve(T);
    std::cout << "iterations:	" << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    #pragma omp parallel for
    for (int i = 0; i < size ; i++)
    {
        if (!g.sphere_cells[i])
        {
            val[i] = output[i];
        } 
    }
    
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::iteration_muscl(std::vector<F> &f, 
                std::vector<std::vector<F>> const &f_face,
                std::vector<std::vector<std::vector<F>>> const &vel)
{
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        F temp = 0.0;
        bool boundary = false;
        for (int j = 0; j < 4; j++)
        {
            if (g.cell_neighbor_ids[i*4+j]==-1)
            {
                boundary = true;
                break;
            }
            F sc_pr = (vel[0][i][j]*g.normals[i*4+j].x 
                    + vel[1][i][j]*g.normals[i*4+j].y 
                    + vel[2][i][j]*g.normals[i*4+j].z);
            int neigh_id = -1;
            for (int m = 0; m < 4; m ++)
	    {
	        if (g.cell_neighbor_ids[g.cell_neighbor_ids[i*4+j]*4+m] 
	            == i)
                {
                    neigh_id = m;
                    break;
                }
            }
            if(neigh_id != -1) 
                 temp -= g.tr_areas[i*4+j]*
                 (std::max(sc_pr, (F) 0.0) * f_face[i][j] 
                  + std::min(sc_pr, (F) 0.0) * 
                  f_face[g.cell_neighbor_ids[i*4+j]][neigh_id]);
            else
            {
                std::cout << "THIS HAPPENED\n";
                temp -= g.tr_areas[i*4+j] * f_face[i][j];
            }
        }
        temp/=g.tet_volume[i];
        if (std::isnan(f[i]))
        {
            std::cout << i << " is nan" << std::endl;
            std::cout << g.cell_centroids[i].x << " " 
                      << g.cell_centroids[i].y << " " 
                      << g.cell_centroids[i].z << std::endl;
        }
        if (f[i] > 50)
        {
            std::cout << i << " " << f[i] << " f is bigger than 50\n";
            std::cout << g.cell_centroids[i].x << " "
                      << g.cell_centroids[i].y << " " 
                      << g.cell_centroids[i].z << std::endl;
        }
        if (!g.boundary_cells[i] && !boundary)
            f[i] -= dt*(temp + f[i]*sigma(g.cell_centroids[i].z+3.0,
                                           (F) 0.05, 20, 100));
    }
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::init_mats()
{
    mat = SpMat(size, size);
    std::vector<Triplet> tripletList;
    tripletList.reserve(100*size);
    for (int i = 0; i < size; i++)
    {
        Eigen::MatrixXd B_inv(3, 3);
        Eigen::VectorXd o(3);
        Eigen::Vector3d xi(g.cell_centroids[i].x, 
                           g.cell_centroids[i].y,
                           g.cell_centroids[i].z);
        int s = g.cell_stencil_ids[i].size();
        for (int j = 0; j < s; j++)
        {
            I new_i = g.cell_stencil_ids[i][j];
            if ((new_i!=-1) && (new_i!=i))
            {
                Eigen::Vector3d xj(g.cell_centroids[new_i].x, 
                                   g.cell_centroids[new_i].y, 
                                   g.cell_centroids[new_i].z);
                F psi = W((xj-xi).norm());
                    
                int s_2 = g.cell_stencil_ids[new_i].size();
                F omega = 0.0;
                for (int k = 0; k < s_2; k++)
                {
                    I fin_i = g.cell_stencil_ids[new_i][k];
                    if ((fin_i!=-1) && (fin_i!=new_i))
                    {
                        Eigen::Vector3d xk(
                                   g.cell_centroids[fin_i].x, 
                                   g.cell_centroids[fin_i].y, 
                                   g.cell_centroids[fin_i].z);
                        omega+= W((xj-xk).norm());
                    }
                }
                psi /= omega;
                o += psi * (xj-xi);
                B_inv +=  psi * (xj-xi)*(xj-xi).transpose();
            }  
        }
        F denominator = 0;
        for (int j = 0; j < s; j++)
        {
            I new_i = g.cell_stencil_ids[i][j];
            if ((new_i!=-1) && (new_i!=i))
            {
                Eigen::Vector3d xj(g.cell_centroids[new_i].x, 
                                   g.cell_centroids[new_i].y, 
                                   g.cell_centroids[new_i].z);
                F term = 1.-(xj-xi).transpose()*(B_inv.inverse()*o);
                F psi = W((xj-xi).norm());
                int s_2 = g.cell_stencil_ids[new_i].size();
                F omega = 0.0; 
                for (int k = 0; k < s_2; k++)
                {
                    I fin_i = g.cell_stencil_ids[new_i][k];
                    if ((fin_i!=-1) && (fin_i!=new_i))
                    {
                        Eigen::Vector3d xk(
                                   g.cell_centroids[fin_i].x, 
                                   g.cell_centroids[fin_i].y, 
                                   g.cell_centroids[fin_i].z);
                        omega+= W((xj-xk).norm());
                    }
                }
                psi /= omega;
                denominator += psi * (xj-xi).norm() *
                               (xj-xi).norm() * term;
            }
        }
        F ith_val = 0;
        for (int j = 0; j < s; j++)
        {
            I new_i = g.cell_stencil_ids[i][j];
            if ((new_i!=-1) && (new_i!=i))
            {
                Eigen::Vector3d xj(g.cell_centroids[new_i].x, 
                                   g.cell_centroids[new_i].y, 
                                   g.cell_centroids[new_i].z);
                F term = 1 - (xj-xi).transpose()*(B_inv.inverse()*o);
                F psi = W((xj-xi).norm());
                int s_2 = g.cell_stencil_ids[new_i].size();
                F omega = 0.0; 
                for (int k = 0; k < s_2; k++)
                {
                    I fin_i = g.cell_stencil_ids[new_i][k];
                    if ((fin_i!=-1) && (fin_i!=new_i))
                    {
                        Eigen::Vector3d xk(
                                   g.cell_centroids[fin_i].x, 
                                   g.cell_centroids[fin_i].y, 
                                   g.cell_centroids[fin_i].z);
                        omega+= W((xj-xk).norm());
                    }
                }
                psi /= omega;
                F value = 6*psi*term/denominator;
                ith_val -= value;
                tripletList.push_back(Triplet(i, new_i, value));
            }
        }
        tripletList.push_back(Triplet(i, i, ith_val));
    }

    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    for (int i = 0; i < size; i++)
    {
        int s = g.cell_stencil_ids[i].size();
        int num_dummies = 0;
        for (int j = 0; j < s; j++) 
            if ( (g.cell_stencil_ids[i][j]==-1) ||
                (g.cell_stencil_ids[i][j]==i) )
                num_dummies++;
        num_dums.push_back(num_dummies);
    }
    return;
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::init_help_vectors()
{
    size = g.size;
    f.resize(size);
    initialize(f, Vector3d<F>{0.0,0.0,-3.0});
    u.resize(3);
    temp_x.resize(3);

    for (int i = 0; i < 3; i++) 
    {
        u[i] = std::vector<F>(size);
        temp_x[i] = std::vector<F>(size);
    }
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        u[2][i] = -1.0;
    }
    temp_face.resize(size);
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        temp_face[i] = std::vector<F>(4);
    temp_face_vec.resize(3);
    for (int i = 0; i < 3; i++)
    {
        temp_face_vec[i] = std::vector<std::vector<F>>(size);
        #pragma omp parallel for
        for (int j = 0; j < size; j++)
            temp_face_vec[i][j] = std::vector<F>(4);
    }
    return;
}

template<typename F, typename I>
SpaceHeter3d_tet<F,I>::~SpaceHeter3d_tet()
{
    std::cout << "destructor\n";
}

template class SpaceHeter3d_tet<double, int>;

