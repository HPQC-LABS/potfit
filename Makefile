# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
#  FC = ifort                  # old alternate compiler on this machine
   FC = f90                    # default compiler on this machine
#  FC = gfortran               # a compiler for testing purposes
endif
#
#option to choose level of optimization:  make DEBUG=1
ifndef DEBUG
    FFLAGS =  -O2 -u            # DEBUG not specified, so optimize
  else                           
    FFLAGS =  -g -u             # DEBUG=1  (basic DEBUG option)
    ifeq ($(DEBUG),2)           
        FFLAGS =  -C -g -u       # DEBUG=2  (higher-level DEBUG option)
    endif
#                 can add more options for  DEBUG=3 , etc., as desired
  endif
#
# as usual, list the objects
#
OBJECTS = masses.o dpotfit.o readata.o tvsort.o readpot.o writepot.o vgen.o dampF.o AFdiag.o mkpredict.o dyidpj.o dvirdp.o vgenp.o diffstats.o prepott.o mappar.o alf.o scecor.o cdjoel.o schrq.o nllssrr.o fununc.o gpround.o PSfunc.o uncbv.o dpsidp.o 

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o tdpot.x

#--------------------------------------------------------------
# Shell Form from Sean Mcleod,  13 February 2008
#--------------------------------------------------------------
