TARGETS = static_amr
all: $(TARGETS)
                                 
clean: 
	rm -f *.o *~ $(TARGETS)
	rm -f *.txt *.out denytempo*.dat Eytempo*.dat den*.dat Erms*.dat diffu*.dat ionfr*.dat panim/*.dat panim/*.png canim/*.dat canim/*.png canim2/*.dat canim2/*.png canim3/*.dat canim3/*.png panimE/*.png panimE/*.dat canimE/*.dat canimE/*.png anim/*.png anim/*.dat animE/*.png animE/*.dat bothdenEx/*.txt  bothdenEy/*.txt parentgrid/*.txt childgrid/*.txt canim/*.dat canim/*.png canimE/*.dat canimE/*.png *.csv
static_amr: xyfdtd.o static_amr.o 
	nvcc -o static_amr xyfdtd.o static_amr.o -lm

static_amr.o: static_amr.cu
	nvcc -c static_amr.cu

xyfdtd.o: xyfdtd.cu
	nvcc -c xyfdtd.cu
