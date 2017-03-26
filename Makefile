GENCODE ?= arch=compute_30,code=sm_30

CC = nvcc
CFLAGS += -gencode $(GENCODE)

ordergraph: ordergraph.cu ordergraph_kernel.cu data.cu
	$(CC) $(CFLAGS) -o $@ $<

node_modules:
	mkdir -p node_modules

node_modules/webppl-json: node_modules
	cd node_modules; git clone https://github.com/stuhlmueller/webppl-json

data.json: data.wppl node_modules/webppl-json
	webppl data.wppl --require webppl-json

data.cu: data.json
	./data.py >data.cu

.PHONY: clean
clean:
	rm -f ordergraph data.{json,cu}
