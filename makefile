ifeq ($(origin CXX),default)
CXX := $(firstword $(wildcard /opt/local/stow/gcc-11.2.0/bin/g++ /opt/local/stow/gcc-15.1.0/bin/g++) g++)
endif
NVCC ?= $(firstword $(wildcard /opt/common/cuda/cuda-11.8.0/bin/nvcc /opt/common/cuda/cuda-12.4.1/bin/nvcc /opt/common/cuda/cuda-12.8.1/bin/nvcc /opt/common/cuda/cuda-13.1.1/bin/nvcc) nvcc)
CUDAHOSTCXX ?= $(CXX)

CXXFLAGS ?= -std=gnu++2a -O3 -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?=
OMPFLAGS ?= -fopenmp
CUDA_ARCH ?= sm_70
NVCCFLAGS ?= -std=c++17 -O3 -arch=$(CUDA_ARCH) -ccbin $(CUDAHOSTCXX) --expt-relaxed-constexpr -Xcompiler -Wall,-Wextra

comma := ,
CXX_ROOT := $(patsubst %/bin/g++,%,$(CXX))
CXX_RPATH := $(if $(wildcard $(CXX_ROOT)/lib64/libstdc++.so.6),-Wl$(comma)-rpath$(comma)$(CXX_ROOT)/lib64)
CUDAHOST_ROOT := $(patsubst %/bin/g++,%,$(CUDAHOSTCXX))
CUDAHOST_RPATH := $(if $(wildcard $(CUDAHOST_ROOT)/lib64/libstdc++.so.6),-Xlinker -rpath -Xlinker $(CUDAHOST_ROOT)/lib64)

OMP_TARGET := bin/solve_omp
OMP_SOURCES := src/solve_omp.cpp
OMP_OBJECTS := $(OMP_SOURCES:src/%.cpp=bin/%.o)

CUDA_TARGET := bin/solve_cuda
CUDA_SOURCES := src/solve_cuda.cu
CUDA_OBJECTS := $(CUDA_SOURCES:src/%.cu=bin/%.cuda.o)

DEPS := $(OMP_OBJECTS:.o=.d) $(CUDA_OBJECTS:.o=.d)

.PHONY: all clean cuda omp run run-cuda

all: omp cuda

omp: $(OMP_TARGET)

cuda: $(CUDA_TARGET)

$(OMP_TARGET): $(OMP_OBJECTS) | bin
	$(CXX) $(LDFLAGS) $(CXX_RPATH) $^ $(LDLIBS) $(OMPFLAGS) -o $@

bin/%.o: src/%.cpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OMPFLAGS) -MMD -MP -c $< -o $@

$(CUDA_TARGET): $(CUDA_OBJECTS) | bin
	$(NVCC) $(NVCCFLAGS) $(CUDAHOST_RPATH) $^ -o $@

bin/%.cuda.o: src/%.cu | bin
	$(NVCC) $(NVCCFLAGS) -MMD -MP -c $< -o $@

bin:
	mkdir -p $@

run: $(OMP_TARGET)
	@if [ -z "$(A_MAX)" ]; then echo "Usage: make run A_MAX=<limit>"; exit 2; fi
	./$(OMP_TARGET) $(A_MAX)

run-cuda: $(CUDA_TARGET)
	@if [ -z "$(A_MAX)" ]; then echo "Usage: make run-cuda A_MAX=<limit> [MAX_SOLUTIONS=<limit>]"; exit 2; fi
	./$(CUDA_TARGET) $(A_MAX) $(MAX_SOLUTIONS)

clean:
	find bin -type f ! -name .gitkeep -delete

-include $(DEPS)
