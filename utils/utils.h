#pragma once
#include <sstream>
#include <vector>
#include <string>
#include <iterator>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>

typedef CGAL::Surface_mesh_default_triangulation_3 Tr2d;
typedef CGAL::Complex_2_in_triangulation_3<Tr2d> C2t3;
typedef Tr2d::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;


template<typename T>
std::vector<T> split(const std::string& line)
{
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is),
                          std::istream_iterator<T>());
}

//===================================================================
//===================================================================
//===================================================================

float sphere(float x, float y, float z)
{
    float x2=x*x, y2=y*y, z2=(z+1.7)*(z+1.7);
    return x2+y2+z2-0.09;
}

float cylinder (float x, float y, float z)
{
    float x2 = x * x;
    float y2 = y * y;
    return (x2 + y2 - 1.);
}

//this is hardcoded to params volcano_to_off
float volcano(float x, float y, float z)
{
    float x2 = x * x;
    float y2 = y * y;
    float r = 0.0;
    if (z > 0.0 && z<=1.0) r = 1 - z/2.;
    else if (z<= 0.0 && z>=-1.0) r = 2 - (z+1);
    else return 1.0;
    return (x2 + y2 - r*r);
}

float thin_foil (float x, float y, float z)
{
    float result = 0.0;
    if ( x<=0.15 && fabs(z)-1.0<0.1 && fabs(y)-1.0<0.1) result = 0.0;
    else result = 1.0;
    return result;
}

//===================================================================
//===================================================================
//===================================================================

float sphere_tet(Point_3 p)
{
    return sphere(p.x(), p.y(), p.z());
    //float x2=p.x()*p.x(), y2=p.y()*p.y(), z2=(p.z()-1.5)*(p.z()-1.5);
    //return x2+y2+z2-0.3;
}

void sphere_to_off()
{
    Tr2d tr;
    C2t3 c2t3 (tr);
    Surface_3 surface(sphere_tet, Sphere_3(CGAL::ORIGIN, 9.));
    CGAL::Surface_mesh_default_criteria_3<Tr2d> criteria(30., 0.02, 0.02);
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    std::ofstream out("offs/sphere.off");
    CGAL::output_surface_facets_to_off(out, c2t3);
    out.close();
}

void cylinder_to_off()
{
    std::ofstream out;
    out.open("offs/cylinder.off");
    int max_t = 40;
    int max_h = 20;
    out << "OFF\n"<< max_t*max_h+2 << " " << 2*(max_h-1)*max_t+2*max_t << " 0\n\n";

    for (int i = 0; i < max_h; i++)
    {
        for (int t = 0; t < max_t; t++)
        {
            out << cos(t*2*3.1416/max_t) << " " << sin(t*2*3.1416/max_t) << " " << i/(max_h-1.)*2. - 1 << std::endl;
        }
    }
    out << "0.0 -1.0 0.0\n0.0 1.0 0.0\n";
    for (int i = 0; i < max_h-1; i++)
    {
        int offset = i*max_t;
        for (int t = 0; t < max_t; t++)
        {
            out << "3  " << t+offset << " " << (t+1)%max_t+offset << " " << t+max_t+offset << "\n";
            out << "3  " << t+max_t+offset << " " << (t+1)%max_t+offset << " " << (t+1)%max_t+max_t+offset << "\n";
        }
    }
    for (int i = 0; i < max_t; i++)
    {
        out << "3  " << (i+1)%max_t << " " << i << " " << max_t*max_h << "\n";
        out << "3  " << max_t*(max_h-1)+i << " " << max_t*(max_h-1)+(i+1)%max_t << " " << max_t*max_h+1 << "\n";
    }
    out.close();
}

void volcano_to_off()
{
    std::ofstream out;
    out.open("offs/volcano.off");
    int max_t = 40;
    int max_h = 20;
    out << "OFF\n"<< max_t*max_h+2 << " " << 2*(max_h-1)*max_t+2*max_t << " 0\n\n";

    //the base is half of the total size and has wider radiuses
    for (int i = 0; i < max_h/2; i++)
    {
        for (int t = 0; t < max_t; t++)
        {
            //radius from 2 to 1 base starts at -1 (-1)
            float r = 2-i/(max_h/2-1.);
            out << r*cos(t*2*3.1416/max_t) << " " << r*sin(t*2*3.1416/max_t) << " " << i/(max_h-1.)*2. - 1 << std::endl;
        }
    }
    //the tip is half of the total size
    //some computations to make seamless transition
    for (int i = max_h/2; i < max_h; i++)
    {
        for (int t = 0; t < max_t; t++)
        {
            //radius from 1 to 0.5 top ends at 1 (-1 + 2)
            float r = 2 - max_h/2/(max_h/2-1.) - (i-max_h/2)/(max_h/2-1.)*0.5;
            out << r*cos(t*2*3.1416/max_t) << " " << i/(max_h-1.)*2. - 1 << " " << r*sin(t*2*3.1416/max_t) << std::endl;
        }
    }

    out << "0.0 -1.0 0.0\n0.0 1.0 0.0\n";
    for (int i = 0; i < max_h-1; i++)
    {
        int offset = i*max_t;
        for (int t = 0; t < max_t; t++)
        {
            out << "3  " << t+offset << " " << (t+1)%max_t+offset << " " << t+max_t+offset << "\n";
            out << "3  " << t+max_t+offset << " " << (t+1)%max_t+offset << " " << (t+1)%max_t+max_t+offset << "\n";
        }
    }
    for (int i = 0; i < max_t; i++)
    {
        out  << "3  " << (i+1)%max_t << " " << i << " " << max_t*max_h << "\n";
        out  << "3  " << max_t*(max_h-1)+i << " " << max_t*(max_h-1)+(i+1)%max_t << " " << max_t*max_h+1 << "\n";
    }
    out.close();
}

void thin_foil_to_off()
{
    std::ofstream out;
    out.open("offs/thin_foil.off");
    int max_y = 50;
    int max_x = 30;
    out << "OFF\n"<< max_x*max_y*2 << " " << 2*(max_x-1)*(max_y-1)+2*(max_x-1)+2*(max_y-1) << " 0\n\n";

    for (int i = 0; i < max_y; i++)
    {
        for (int j = 0; j < max_x; j++)
        {
            out << 0.0  << " " << j/(max_x-1.)*4. - 2 << " " << i/(max_y-1.)*4. - 2 << std::endl;
        }
    }
    for (int i = 0; i < max_y; i++)
    {
        for (int j = 0; j < max_x; j++)
        {
            out << 0.1  << " " << j/(max_x-1.)*4. - 2 << " " << i/(max_y-1.)*4. - 2 << std::endl;
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        int offset = i*max_x;
        for (int j = 0; j < max_x-1; j++)
        {
            out << "3  " << j+offset << " " << j+1+offset << " " << j+max_x+offset << "\n";
            out << "3  " << j+max_x+offset << " " << j+1+offset << " " << j+1+max_x+offset << "\n";
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        int offset = i*max_x+max_x*max_y;
        for (int j = 0; j < max_x-1; j++)
        {
            out << "3  " << j+offset << " " << j+1+offset << " " << j+max_x+offset << "\n";
            out << "3  " << j+max_x+offset << " " << j+1+offset << " " << j+1+max_x+offset << "\n";
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        out << "3  " << i*max_x << " " << i*max_x+max_x << " " << i*max_x+max_x*max_y << "\n";
        out << "3  " << i*max_x+max_x*max_y << " " << i*max_x+max_x << " " << i*max_x+max_x*max_y+max_x << "\n";
        out << "3  " << i*max_x+max_x-1 << " " << i*max_x+max_x+max_x-1 << " " << i*max_x+max_x-1+max_x*max_y << "\n";
        out << "3  " << i*max_x+max_x-1+max_x*max_y << " " << i*max_x+max_x+max_x-1  << " " << i*max_x+max_x-1+max_x*max_y+max_x << "\n";
    }
    for (int j = 0; j < max_x-1; j++)
    {
        out << "3  " << j << " " << j+1 << " " << j+max_x*max_y << "\n";
        out << "3  " << j+max_x*max_y << " " << j+1 << " " << j+max_x*max_y+1 << "\n";
        out << "3  " << j+max_x*(max_y-1) << " " << j+1+max_x*(max_y-1) << " " << j+max_x*(max_y-1)+max_x*max_y << "\n";
        out << "3  " << j+max_x*(max_y-1)+max_x*max_y << " " << j+1+max_x*(max_y-1) << " " << j+max_x*(max_y-1)+max_x*max_y+1 << "\n";
    }
    out.close();
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
