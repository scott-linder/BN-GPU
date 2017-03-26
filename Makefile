GENCODE ?= arch=compute_30,code=sm_30

CC = nvcc
CFLAGS += -gencode $(GENCODE)

ordergraph: ordergraph.cu
	$(CC) $(CFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f ordergraph
