# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
   FC = f90                    # default compiler on this machine
#  FC = gfortran               # a compiler for testing purposes
endif
#
#option to choose level of optimization:  make debug=1
ifndef debug
    FFLAGS =  -O2             # debug not specified, so optimize
  else                           
#   FFLAGS =  -q  -C          # debug=1  (basic gtortran debug option)
#    FFLAGS =  -g     -u       # debug=1  (basic debug option)
    ifeq ($(debug),2)
        FFLAGS =  -C -g        # debug=2  (higher-level debug option)
    endif
#                 can add more options for  debug=3 , etc., as desired
  endif
#
# as usual, list the objects
#
OBJECTS = dpotfit.o masses16.o readata.o tvsort.o readpot.o writepot.o vgen.o quadCORR.o dampF.o AFdiag.o mkpredict.o dyidpj.o dvirdp.o dvacdp.o vgenp.o diffstats.o prepott.o mappar.o alf.o scecor.o cdjoel.o schrq.o nllssrr.o fununc.o gpround.o PSfunc.o uncbv.o

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o tdpot.x

#--------------------------------------------------------------
# Shell Form from Sean Mcleod,  13 February 2008
#--------------------------------------------------------------
