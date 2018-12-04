# Makefile

FC = 		mpiifort  # MPI - Fortran 90 compiler
CC = 		mpiicc   # C compiler
OPT = 		-O3 -fpp -qopt-report=0
#-vec-report0 # -ipo -noalign -fp-model source -w
MMOD = 		-module
MODULEDIR = 	./Modules/
F90FLAGS =	$(MMOD) $(MODULEDIR)

# patsubst = path substitution, wildcard = tout
# si le repertoire contient Fn1.c Fn2.c Fn3.f90 et Fn4.f90 alors OBJ = "Fn1.o Fn2.o Fn3.o Fn4.o"
OBJ = 		$(patsubst %.c,%.o, $(wildcard *.c))  $(patsubst %.f90,%.o, $(wildcard *.f90)) 
OBJ2 = 		$(patsubst %.f,%.o, $(wildcard BLAS/*.f)) 
OBJ_MOD = 	$(patsubst %.c,%.o, $(wildcard $(MODULEDIR)*.c)) $(patsubst %.f90,%.o, $(wildcard $(MODULEDIR)*.f90))

all : spec.x

Modules/Domain.o :	$(MODULEDIR)Element.o $(MODULEDIR)Face.o $(MODULEDIR)Source.o $(MODULEDIR)Receiver.o \
                        $(MODULEDIR)Vertex.o $(MODULEDIR)TimeParam.o $(MODULEDIR)Subdomain.o $(MODULEDIR)Edge.o \
                        $(MODULEDIR)Adjoint.o $(MODULEDIR)Comm.o $(MODULEDIR)Mirror.o $(MODULEDIR)External_Source.o
Modules/Element.o :	$(MODULEDIR)Simu.o
Modules/Face.o :	$(MODULEDIR)Simu.o
Modules/Edge.o :	$(MODULEDIR)Simu.o
Modules/Vertex.o :	$(MODULEDIR)Simu.o
Modules/Source.o :	$(MODULEDIR)angles.o
Modules/nrutil.o :	$(MODULEDIR)nrtype.o
Modules/module_spline.o :	$(MODULEDIR)nrutil.o $(MODULEDIR)nrtype.o
Modules/module_sort.o : 	$(MODULEDIR)nrutil.o $(MODULEDIR)nrtype.o
Modules/init_cond.o :	$(MODULEDIR)module_spline.o
Modules/module_ellipticity.o : $(MODULEDIR)module_A3d.o $(MODULEDIR)module_spline.o $(MODULEDIR)funaro.o
Modules/module_A3d.o :	$(MODULEDIR)def_gparam.o $(MODULEDIR)earth_modele.o \
                        $(MODULEDIR)module_spline.o $(MODULEDIR)spl_A3d.o


spec.x : make_module make_BLAS $(OBJ)
	$(FC) $(F90FLAGS) $(OPT) -o ../bin_bup/$@ $(OBJ) $(OBJ2) $(OBJ_MOD)

make_module : 	$(OBJ_MOD)

make_BLAS :	$(OBJ2)

%.o : %.f90 
	$(FC) $(F90FLAGS) $(OPT) -c $< -o $@

%.o : %.c
	$(CC) -O2 -w -c $< -o $@


clean: clean_modules clean_BLAS cleanx
	rm -f *.o

clean_modules:
	rm -rf  $(MODULEDIR)/*.o $(MODULEDIR)/*.mod

clean_BLAS:
	rm -f BLAS/*.o

cleanx:
	rm -f ../bin_bup/*.x
