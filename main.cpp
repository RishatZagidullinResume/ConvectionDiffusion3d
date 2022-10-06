#include "glad/glad.h"
#include "utils/utils.cuh"
#include "utils/utils.h"
#include "solvers/advection_tet.h"
#include "solvers/advection_reg.h"
#include "solvers/diffusion_reg.h"
#include "geometry/geometry_tet.h"
#include "geometry/projection_reg.h"
#include "gl/gl_viewer.h"	
#include "gl/cuda_gl.cuh"
#include "gl/animation_tet.h"
#include "gl/animation_reg.h"
#include <time.h>
#include <sys/time.h>

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char ** argv)
{
    //sphere_to_off();
    //return 0;
    double start = get_wall_time();

    short SCR_WIDTH = 1600;
    short SCR_HEIGHT = 1200;

    //true is to hide window
    custom_gl::gl_viewer viewer(SCR_WIDTH, SCR_HEIGHT, true);
    geometry3d<double, int> geom("./offs/sphere.off",
                                 "./offs/boundary.off");

    geometry_tet grid(geom.tet, &sphere);

    grid_params proj_grid = grid_params(20,20,50,
                                        -1.0,-1.0,-3.0,0.1);
    projection_reg proj(proj_grid);

    double dt = 0.001;
    SpaceHeter3d_tet<double, int> eqn(geom, dt);

    viewer.buffer_shape_vbo(grid.vert_data_size*4, grid.vertices);
    viewer.buffer_shape_ebo(grid.index_data_size,
                            grid.indices_to_draw.data());

    viewer.buffer_vbo(proj.max_v*4, proj.vertices);
    viewer.buffer_ebo(proj.max_i*6, proj.indices);
    viewer.preset();

    custom_gl::cuda_gl shape_interop(&viewer.VBO2,
                               grid.vert_data_size);
    custom_gl::cuda_gl interop(&viewer.VBO, proj.max_v);
    
    AnimationData_tet colors(proj, eqn);
    int time = 0;
    int TIME_MAX = 8001;
    std::cout <<"Preprocessing time: "<<get_wall_time()-start<<"\n";
    start = get_wall_time();
    while (!glfwWindowShouldClose(viewer.window))
    {
        if (time % 100 == 0)
        {
            colors.update_mesh_colors();
            shape_interop.update_gpu_data(colors.colors_mesh);

            colors.update_proj_colors();
            interop.update_gpu_data(colors.colors_proj);
            viewer.view(proj.max_i*6,
                        grid.index_data_size, time);
        }

        eqn.iteration_advection(eqn.u, eqn.f);

        if (time==TIME_MAX)
        {
            time=0;
            break;
        }

        time++;
        printProgress((double)time/TIME_MAX);
    }
    std::cout <<"\nComputation time: "<< get_wall_time()-start<<"\n";
    return 0;
}



