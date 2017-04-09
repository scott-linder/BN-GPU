GENCODE ?= arch=compute_30,code=sm_30

NVCC ?= nvcc
NVCFLAGS += -gencode $(GENCODE)

CXX ?= g++
CPPFLAGS ?= -std=c++11

.PHONY: all
all: ordergraph_gpp ordergraph_gpu

ordergraph_gpp: ordergraph.cc data.h
	$(CXX) $(CPPFLAGS) -o $@ $<

ordergraph_gpu: ordergraph.cu ordergraph_kernel.cu data.h
	$(NVCC) $(NVCFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f ordergraph_gpu ordergraph_gpp
