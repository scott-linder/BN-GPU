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

node_modules:
	mkdir -p node_modules

node_modules/webppl-json: node_modules
	cd node_modules; git clone https://github.com/stuhlmueller/webppl-json

data.json: data.wppl node_modules/webppl-json
	webppl data.wppl --require webppl-json

data.h: data.json
	./data.py >$@

.PHONY: clean
clean:
	rm -f ordergraph{_gpp,gpu} data.{json,h}
