COMPILER = ifort
COMPOPTS = -O3 -unroll -fno-alias -fPIC -vec-report0 -mcmodel=medium -i-dynamic -i8
LINKOPTS = -L/usr/local/lib -lrfftw -lfftw -lm

res_to_dat.o: res_to_dat.f fft.o\
	        header grid_def
	$(COMPILER) $(COMPOPTS) res_to_dat.f -o res_to_dat \
	fft.o $(LINKOPTS)  -L/uc/local/lib -lrfftw

fft.o:  fft.f header grid_def
	$(COMPILER) $(COMPOPTS) -c fft.f

clean:
	rm -f *.o fort.* *~ diablo core