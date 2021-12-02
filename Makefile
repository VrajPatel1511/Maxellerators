objects = main.o particle.o v3.o

all: $(objects)
    nvcc -arch=sm_20 $(objects) -o app

%.o: %.cpp
    nvcc -x cu -arch=sm_20 -I. -dc $< -o $@

clean:
    rm -f *.o app
                                 
clear: 
	rm -f *.o *~ $(TARGETS)
	rm -f *.txt *.out denytempo*.dat Eytempo*.dat den*.dat Erms*.dat diffu*.dat ionfr*.dat panim/*.dat panim/*.png canim/*.dat canim/*.png canim2/*.dat canim2/*.png canim3/*.dat canim3/*.png panimE/*.png panimE/*.dat canimE/*.dat canimE/*.png anim/*.png anim/*.dat animE/*.png animE/*.dat bothdenEx/*.txt  bothdenEy/*.txt parentgrid/*.txt childgrid/*.txt canim/*.dat canim/*.png canimE/*.dat canimE/*.png *.csv

