#pragma once
#include "vector3d.h"
#include <cmath>
#include <vector>

#include <CGAL/centroid.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/centroid.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/make_mesh_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef CGAL::Point_3<K> Point_3;
typedef C3t3::Cell_handle Cell_handle;
typedef Tr::Finite_cells_iterator Cell_iterator;

using namespace CGAL::parameters;

/* cell is a tetrahedron
 *
 *
 */

template<typename F, typename I>
class geometry3d {
public:
    Tr tet;
    I size;
    //inner cell has 4 facet neighbors
    I * cell_neighbor_ids;
    //cell's stencil neighbors share a facet, line or point
    std::vector<I> * cell_stencil_ids;
    //cell has 4 facets, here are it's areas
    F * tr_areas;
    //cell has a volume
    F * tet_volume;
    //normals point out from 4 facets
    Vector3d<F> * normals;
    //facet centroids
    Vector3d<F> * facet_centroids;
    //cell centroids
    Vector3d<F> * cell_centroids;
    //boundary cells that don't take part in the numerical scheme
    std::vector<bool> boundary_cells;
    std::vector<bool> sphere_cells;

    //geometry parameters for the FVM in 3D
    F * beta1_p;
    F * beta2_p;
    F * beta3_p;
    F * beta1_m;
    F * beta2_m;
    F * beta3_m;

    I * a_p;
    I * b_p;
    I * c_p;
    I * a_m;
    I * b_m;
    I * c_m;

    F * dist_center_face;
    F * dist_inter_p;
    F * dist_inter_m;

    std::vector<I> * sorted_ids;

    geometry3d(std::string smin, std::string smout);
    ~geometry3d();
    void init_geometry(Tr const &tet);
    void init_fvm();
};

class geometry_tet {
public:
    Tr & tet;
    float * vertices;
    unsigned int * indices;
    int vert_data_size;
    int index_data_size;
    std::vector<unsigned int> indices_to_draw;
    float (*shape)(float, float, float);
    geometry_tet(Tr & tet, float (*shape)(float, float, float)) : tet(tet), shape(shape)
    {
        init_data();
    }
    ~geometry_tet();
    void out_to_file(std::string vfile, std::string ffile);
private:
    void init_data();
};
