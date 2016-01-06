#USE THIS MAKEFILE TO BUILD NON LINEAR EQUATION SOLVER PACKAGE 
#with LF95
#MPI = 1 # comment out this line if MPI is not used
#OPEN_MP = 1  #comment out this line if Open Mp is not used
#AUTO_PARA = 1  # comment out this line if auto parallelization is not desired

SHELL     = /bin/sh
NAME      = non_lin_test 


OBJECTS   =  nedriver.o nefn.o neinck.o nemodel.o nemodelfac.o \
             nestop0.o broyfac.o dogdriver.o dogstep.o fdjac.o \
             hooksteporig.o hookdriver_orig.o  jacrotate.o \
              linesearchorig.o qrupdate_mod.o com.o copy.o terminate.o \
	     machine_epsd.o broyunfac.o check_convergence.o \
	     rsolve.o condest.o sparse_jac.o \
	     choldecomp.o cholsolve.o euclidnorm.o \
	     trust_region.o lsolve.o sdot.o  \
	     rosenbrock.o jac_rosen.o powell.o jac_powell.o \
             main.o allocate.o  helicalValley.o \
	     jac_helicalValley.o broyden.o jac_broyden.o

#following is for debug only:
#LIBRARIES = /usr/local/pgi3.2/linux86/lib/liblapack.a \
	/usr/local/pgi3.2/linux86/lib/libblas.a

#for a21p machine using rpm lapack and blas requires linking to lg2c:
LIBRARIES      = /usr/lib/liblapack.a /usr/lib/libblas.a  \
	/usr/lib/gcc-lib/i586-mandrake-linux/2.95.3/libg2c.a \

LIBRARIES	=/usr/lib/liblapack.a /usr/lib/libblas.a  /usr/lib/libg2c.a 	/usr/lib/gcc-lib/i586-suse-linux/3.3.5/libgcc.a
LIBRARIES       = -L/usr/local/atlas/Linux_P4/lib  -llapack -lcblas -lf77blas -latlas \
			 /usr/lib/gcc/i586-suse-linux/4.0.2/libgcc.a
LIBRARIES       = -L/usr/local/lf9562/lib/  -llapackmt -lblasmt  \
			 /usr/lib/gcc/i586-suse-linux/4.0.2/libgcc.a
DEBUG     = -g -Mbounds
OPTIMIZE  = -fast
LISTING   = -Mlist
SPECIAL   = -Msave -byteswapio -r8
COMP      =/usr/local/mpi/bin/mpif90 
OPTION = $(DEBUG)
#OPTION = $(OPTIMIZE)
ifndef OPEN_MP 
    OMP =
#    FTN90LIB = /usr/local/lib/libpgftnrtl.a \
		/usr/local/lf9555/lib/libfj9e6.a \
		/usr/local/lf9555/lib/libfj9f6.a \
		/usr/local/lf9555/lib/libfj9i6.a 

#FTN90LIB = /usr/local/lib/libpgftnrtl.a  /usr/local/lib/libpgthread.a \
                                                 /usr/local/lib/libpgmp.a

else
    OMP = -mp
    FTN90LIB = /usr/local/lib/libpgftnrtl.a  /usr/local/lib/libpgthread.a \
                                                 /usr/local/lib/libpgmp.a
endif
ifndef AUTO_PARA
    AUTP =
else
    AUTP = -Mconcur -Minfo
    FTN90LIB = /usr/local/lib/libpgftnrtl.a  /usr/local/lib/libpgthread.a \
                                                 /usr/local/lib/libpgmp.a
endif
ifndef  MPI
  MSG = nodlines
  FC = $(COMP) 
  CXX = gcc 
  INC = ./
else
  MSG = dlines
# compiler specs, for Portland Group Compiler:
  FC= /usr/local/mpich-1.2.1/bin/mpif90 
  CXX= /usr/local/mpich-1.2.1/bin/mpicc 
  INC=/usr/local/mpich-1.2.1/include 
endif

#COMPILE77   = pgf77 -c $(OPTIMIZE) $(LISTING) $(SPECIAL)
COMPILE   = $(FC) -I$(INC) -L$(LIBS)  -M$(MSG) $(OMP) $(AUTP) -c \
	    	 $(OPTION) $(LISTING) $(SPECIAL)
LOAD      = $(FC)   -o $(NAME) # -m
ARCHIVE   = ar rv

PROTECT   = chmod 755
DELETE    = rm -f
RENAME    = echo 
.SUFFIXES:

all: $(NAME)
	@echo
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"
	@echo
	@echo make of $(NAME) in `pwd` on `hostname` completed

$(NAME):           $(OBJECTS)  
	$(LOAD)    $(OBJECTS)   $(LIBRARIES) $(FTN90LIB)
	$(PROTECT) $(NAME)


broyunfac.o:         broyunfac.f com.o copy.o terminate.o \
	            machine_epsd.o
	$(COMPILE) broyunfac.f

nedriver.o:         nedriver.f com.o copy.o terminate.o \
	            machine_epsd.o check_convergence.o nestop0.o \
		    allocate.o
	$(COMPILE) nedriver.f
	$(RENAME)  nedriver.lst nedriver.l
nefn.o:         nefn.f         com.o
	$(COMPILE) nefn.f
	$(RENAME)  nefn.lst nefn.l
neinck.o:         neinck.f   
	$(COMPILE) neinck.f
	$(RENAME)  neinck.lst neinck.l
nemodel.o:         nemodel.f  
	$(COMPILE) nemodel.f
	$(RENAME)  nemodel.lst nemodel.l
nemodelfac.o:         nemodelfac.f  
	$(COMPILE) nemodelfac.f
	$(RENAME)  nemodelfac.lst nemodelfac.l
nestop0.o:         nestop0.f 
	$(COMPILE) nestop0.f
	$(RENAME)  nestop0.lst nestop0.l
check_convergence.o:         check_convergence.f 
	$(COMPILE) check_convergence.f
	$(RENAME)  check_convergence.lst check_convergence.l
condest.o:         condest.f 
	$(COMPILE) condest.f
	$(RENAME)  condest.lst condest.l
choldecomp.o:         choldecomp.f 
	$(COMPILE) choldecomp.f
	$(RENAME)  choldecomp.lst choldecomp.l
cholsolve.o:         cholsolve.f 
	$(COMPILE) cholsolve.f
	$(RENAME)  cholsolve.lst cholsolve.l
broyfac.o:         broyfac.f 
	$(COMPILE) broyfac.f
	$(RENAME)  broyfac.lst broyfac.l
dogdriver.o:         dogdriver.f  
	$(COMPILE) dogdriver.f
	$(RENAME)  dogdriver.lst dogdriver.l
dogstep.o:         dogstep.f  
	$(COMPILE) dogstep.f
	$(RENAME)  dogstep.lst dogstep.l
euclidnorm.o:         euclidnorm.f  
	$(COMPILE) euclidnorm.f
	$(RENAME)  euclidnorm.lst euclidnorm.l
fdjac.o:         fdjac.f 
	$(COMPILE) fdjac.f
	$(RENAME)  fdjac.lst fdjac.l
sparse_jac.o:         sparse_jac.f 
	$(COMPILE) sparse_jac.f
	$(RENAME)  sparse_jac.lst sparse_jac.l
hooksteporig.o:         hooksteporig.f   
	$(COMPILE) hooksteporig.f
	$(RENAME)  hooksteporig.lst hooksteporig.l
hookdriver_orig.o:         hookdriver_orig.f   
	$(COMPILE) hookdriver_orig.f
	$(RENAME)  hookdriver_orig.lst hookdriver_orig.l
jacrotate.o:         jacrotate.f   
	$(COMPILE) jacrotate.f
	$(RENAME)  jacrotate.lst jacrotate.l
linesearchorig.o:         linesearchorig.f   
	$(COMPILE) linesearchorig.f
	$(RENAME)  linesearchorig.lst linesearchorig.l
lsolve.o:         lsolve.f   
	$(COMPILE) lsolve.f
	$(RENAME)  lsolve.lst lsolve.l
qrupdate_mod.o:         qrupdate_mod.f   
	$(COMPILE) qrupdate_mod.f
	$(RENAME)  qrupdate_mod.lst qrupdate_mod.l
rsolve.o:         rsolve.f 
	$(COMPILE) rsolve.f
	$(RENAME)  rsolve.lst rsolve.l
machine_epsd.o:         machine_epsd.f   
	$(COMPILE) machine_epsd.f
	$(RENAME)  machine_epsd.lst machine_epsd.l

main.o:         main.f   allocate.o
	$(COMPILE) main.f
	$(RENAME)  main.lst main.l

trust_region.o:         trust_region.f   allocate.o
	$(COMPILE) trust_region.f
	$(RENAME)  trust_region.lst trust_region.l

com.o:         com.f90   
	$(COMPILE) com.f90
	$(RENAME)  com.lst com.l

copy.o:         copy.f90   
	$(COMPILE) copy.f90
	$(RENAME)  copy.lst copy.l

sdot.o:         sdot.f   
	$(COMPILE) sdot.f
	$(RENAME)  sdot.lst sdot.l
rosenbrock.o:         rosenbrock.f  
	$(COMPILE) rosenbrock.f
	$(RENAME)  rosenbrock.lst rosenbrock.l
jac_rosen.o:         jac_rosen.f 
	$(COMPILE) jac_rosen.f
	$(RENAME)  jac_rosen.lst jac_rosen.l
powell.o:         powell.f  
	$(COMPILE) powell.f
	$(RENAME)  powell.lst powell.l
jac_powell.o:         jac_powell.f 
	$(COMPILE) jac_powell.f
	$(RENAME)  jac_powell.lst jac_powell.l

helicalValley.o:         helicalValley.f  
	$(COMPILE) helicalValley.f
	$(RENAME)  helicalValley.lst helicalValley.l
jac_helicalValley.o:         jac_helicalValley.f 
	$(COMPILE) jac_helicalValley.f
	$(RENAME)  jac_helicalValley.lst jac_helicalValley.l

broyden.o:         broyden.f  
	$(COMPILE) broyden.f
	$(RENAME)  broyden.lst broyden.l
jac_broyden.o:         jac_broyden.f 
	$(COMPILE) jac_broyden.f
	$(RENAME)  jac_broyden.lst jac_broyden.l


terminate.o:         terminate.f90   
	$(COMPILE) terminate.f90
	$(RENAME)  terminate.lst terminate.l 

allocate.o:         allocate.f90   
	$(COMPILE) allocate.f90
	$(RENAME)  allocate.lst allocate.l
clean:
	$(DELETE) $(NAME) *.a *.o
