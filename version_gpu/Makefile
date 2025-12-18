
# Makefile for MCML GPU Version

NVCC = nvcc
CFLAGS = -O3 -lm
ARCH = -arch=sm_75 # Adjusted to Turing

all: mcml_gpu

mcml_gpu: mcml_gpu.cu mcmlio_gpu.c mcml_gpu.h
	$(NVCC) $(ARCH) -o mcml_gpu mcml_gpu.cu mcmlio_gpu.c $(CFLAGS)

clean:
	rm -f mcml_gpu
