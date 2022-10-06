# Convection-Diffusion solver library in 3D using finite volume and finite difference approaches.

## Code description

This is a 3D library for solving convection-diffusion equation. You can choose to use finite difference method with the regular computational domain and finite volume method for tetrahedralized computational domain.

* **geometry** folder contains classes for generation of computational domains, classes for creating projections for visualization and a `Vector3d` structure. Computational domains can be regular and tetrahedralized (using **CGAL**). Some functions are given to have special geometry constrains (sphere and cylinder).

* **gl** folder contains classes for animation (using **OpenGL**) and direct GPU data conversion between **CUDA** and **OpenGL** (to render animations with Nvidia graphic cards use `__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia` environment variable settings when launching the executable and have a look at `cuda_gl_interop.h` header file provided by Nvidia). Also the folder contains shaders for OpenGL. You can do real time rendering without dumping any data or you can swith the OpenGL window off (third argument in `custom_gl::gl_viewer` constructor). **glad** folder contains util files for OpenGL loader.

* **solvers** folder contains numerical solvers: finite difference method for advection, diffusion equations; finite volume method for advection, diffusion and Poisson equations.

* **utils** folder contains functions for generating `off` files. They are used by **CGAL** to generate geometry constraints.

## Demo

* Here's an example of the tetrahedralized computational domain:

![data_0.jpg](/data_0.jpg)

* Here's the solution to the advection equation using finite volume method. We also add projections to see the solutions as well as keeping the constrained sphere:

![image.gif](/image.gif)


