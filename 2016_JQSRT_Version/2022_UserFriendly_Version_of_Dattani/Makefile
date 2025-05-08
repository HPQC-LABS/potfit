# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
#  FC = f90                    # default compiler on this machine
   FC = gfortran               # a compiler for testing purposes
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
OBJECTS = 1-s2.0-S0022407316301078-mmc2.o

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o pot.x

#--------------------------------------------------------------
# Shell Form from Sean Mcleod,  13 February 2008
#--------------------------------------------------------------
# To run the code, execute:
# ./pot.x < input.5

# input.5 must have lines 2-3 like the following:
# './2016_04_28.4'             % Name of PotFit input data file
# 'output'                     % Writefile"
