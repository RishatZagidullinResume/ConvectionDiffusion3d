#pragma once
#include "../geometry/vector3d.h"
#include "../geometry/geometry_tet.h"
#include "../geometry/projection_reg.h"
#include "../solvers/advection_tet.h"
#include "../solvers/advection_reg.h"

class AnimationData_tet {
public:
    float * colors_mesh;
    float * colors_proj;
    std::vector<int> * inds;
    projection_reg &proj;
    SpaceHeter3d_tet<double, int> &eqn;
    AnimationData_tet(projection_reg & proj, 
        SpaceHeter3d_tet<double, int> & eqn) : proj(proj), eqn(eqn)
    {
       init_color_data();
    }
    void update_proj_colors();
    void update_mesh_colors();
    ~AnimationData_tet();
private:
    void init_color_data();
};

void AnimationData_tet::update_proj_colors()
{
    #pragma omp parallel for
    for (int i = 0; i < proj.max_v; i++)
    {
        colors_proj[i] = 0.0;
        for (int j = 0; j < inds[i].size(); j++)
        {
            colors_proj[i] += eqn.f[inds[i][j]];
        }
        if (inds[i].size() > 0)
            colors_proj[i] /= inds[i].size();
    }
}

void AnimationData_tet::update_mesh_colors()
{
    int c = 0;
    for (auto it(eqn.g.tet.finite_vertices_begin());
              it!=eqn.g.tet.finite_vertices_end(); it++)
    {
        colors_mesh[c] = 0.0;
        std::vector<Cell_handle> finite_incident_cells;
        eqn.g.tet.finite_incident_cells(it, 
                    std::back_inserter(finite_incident_cells));
        for (int i = 0; i < finite_incident_cells.size(); i++)
            colors_mesh[c] += 
              eqn.f[finite_incident_cells[i]->subdomain_index()-1];
        colors_mesh[c] /= ((float) finite_incident_cells.size());
        c++;
    }
}

void AnimationData_tet::init_color_data()
{
    inds = new std::vector<int> [proj.max_v];
    colors_mesh = new float[eqn.g.tet.number_of_vertices()];
    colors_proj = new float[proj.max_v];

    for (int ind=0; ind<eqn.g.tet.number_of_finite_cells(); ind++)
    {
        float val = 0.16;
        float x = eqn.g.cell_centroids[ind].x;
        float y = eqn.g.cell_centroids[ind].y;
        float z = eqn.g.cell_centroids[ind].z;
        if (x*x < val)
        {
            int ind_y = (int) ((y-proj.grid.left_w)/proj.grid.dd);
            int ind_z = (int) ((z-proj.grid.left_h)/proj.grid.dd);
            int proj_ind = ind_z*proj.grid.w+ind_y;
            inds[proj_ind].push_back(ind);
        }
        if (z*z < val)
        {
            int ind_x = (int) ((x-proj.grid.left_d)/proj.grid.dd);
            int ind_y = (int) ((y-proj.grid.left_w)/proj.grid.dd);
            int offset = proj.grid.w*proj.grid.h;
            int proj_ind = ind_y*proj.grid.d+ind_x+offset;
            inds[proj_ind].push_back(ind);   
        }
        if (y*y < val)
        {
            int ind_x = (int) ((x-proj.grid.left_d)/proj.grid.dd);
            int ind_z = (int) ((z-proj.grid.left_h)/proj.grid.dd);
            int offset = proj.grid.w*proj.grid.h + 
                         proj.grid.w*proj.grid.d;
            int proj_ind = ind_z*proj.grid.d+ind_x+offset;
            inds[proj_ind].push_back(ind);
        }
    }
    return;
}

AnimationData_tet::~AnimationData_tet()
{
    delete [] colors_mesh;
    delete [] colors_proj;
    delete [] inds;
}
