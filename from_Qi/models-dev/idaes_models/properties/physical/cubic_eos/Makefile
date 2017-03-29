
ASL=solvers
ASL_BUILD=$(wildcard $(ASL)/sys.*)

CC = gcc
RM = rm
CP = cp

CFLAGS = -Wall -Wextra -O3 -I$(ASL_BUILD)
LDFLAGS = -shared -lm

ALL: phys_prop.so test_prop

test_prop: test_prop.c
	$(CC) $(CFLAGS) -Wl,--no-as-needed -ldl test_prop.c -o test_prop

phys_prop.o: phys_prop.c $(ASL) $(ADOL)
	$(CC) -c $(CFLAGS) -fPIC phys_prop.c -o phys_prop.o

phys_prop.so: phys_prop.o
	$(CC) $(LDFLAGS) phys_prop.o -o phys_prop.so

$(ASL):
	if ! test -e solvers.tar; then curl 'http://netlib.sandia.gov/cgi-bin/netlib/netlibfiles.tar?filename=netlib/ampl/solvers' -o solvers.tar; fi
	tar -xf solvers.tar
	cd $(ASL) \
	    && chmod 755 configure* \
	    && bash ./configure \
	    && make \
	      
clean:
	rm -f *.o
	rm -f *.so
	rm -f test_prop
	rm -rf $(ASL)
	


