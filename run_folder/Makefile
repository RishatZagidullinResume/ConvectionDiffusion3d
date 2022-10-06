GCC = g++
NVCC = nvcc

LIB  = -lgomp -lblas -lm -lgmp -lmpfr -lm -ldl -lglfw

GCC_COPT = -fopenmp

VPATH = ../:../glad:../utils:../solvers:../gl:../geometry

SOURCES_SOL = advection_reg.cpp advection_tet.cpp geometry_tet.cpp diffusion_reg.cpp vector3d.cpp coagulation.cpp
SOURCES_GL =  glad.c
SOURCES_MAIN = main.cpp
SOURCES_CU = cuda_gl.cu

OBJECTS_SOL = $(SOURCES_SOL:%.cpp=$(OBJDIR)/%.o)
OBJECTS_GL = $(SOURCES_GL:%.c=$(OBJDIR)/%.o)
OBJECTS_MAIN = $(SOURCES_MAIN:%.cpp=$(OBJDIR)/%.o)
OBJECTS_CU = $(SOURCES_CU:%.cu=$(OBJDIR)/%_cu.o)

OBJDIR = build

all: tet.exe

tet.exe: $(OBJECTS_MAIN) $(OBJECTS_GL) $(OBJECTS_SOL) $(OBJECTS_CU)
	$(NVCC) $^ -o $@ $(LIB)

$(OBJDIR)/%.o: %.cpp
	$(GCC) -c -I/usr/local/cuda-11.6/include -I/usr/include/eigen3 $(GCC_COPT) $< -o $@

$(OBJDIR)/%.o: %.c
	$(GCC) -c $< -o $@

$(OBJDIR)/%_cu.o: %.cu
	$(NVCC) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o *.exe
