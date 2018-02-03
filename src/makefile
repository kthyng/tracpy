###############################################################
### TEMPLATE FOR Mac OS X 10.7 with gfortran from fink ###
###############################################################

#F95COMPILER        = "g95"
F95COMPILER        = "gfortran"

#MYCFG            = /usr/local/mysql/bin/mysql_config
#MYI_FLAGS        = `$(MYCFG) --cflags`
#MYL_FLAGS        = `$(MYCFG) --libs`

# LIB_DIR           = -L/sw/lib -L/sw/opt/netcdf7/lib
# INC_DIR           = -I/sw/include -I/sw/opt/netcdf7/include \
#                     -I/usr/local/mysql/include

# LNK_FLAGS         = -lnetcdf -lnetcdff

#================================================================

# VPATH = src:projects/$(PROJECT)
# vpath %.o tmp

ifeq ($(F95COMPILER),"gfortran")
	FF_FLAGS         = -m64 -c -x f95-cpp-input -fconvert=big-endian -gdwarf-2 -fbounds-check -cpp -fPIC
	# F90_FLAGS        =-fno-underscoring
#	FF               = gfortran $(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)
	# FF               = /sw/bin/gfortran $(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)
	FF               = gfortran $(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)

endif
# ifeq ($(F95COMPILER),"g95")
# 	FF_FLAGS = -c -cpp -fendian=big
# 	F90_FLAGS        = -O3 -C  -g  -fno-underscoring
# 	FF               = /Applications/fort/g95/bin/i386-apple-darwin8.11.1-g95 $(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)
# endif
CC                = gcc -O  $(INC_DIR)

COMPUTER = $(shell uname -n)

ifeq ($(COMPUTER),rainier)
	F2PY = f2py
else ifeq ($(findstring hafen,$(COMPUTER)),hafen)
	F2PY = f2py
else
	F2PY = f2py
endif

objects           = pos.o cross.o calc_dxyz.o calc_time.o loop_pos.o \
					vertvel.o turb.o diffusion.o
f2py_source       = outdir/step.f95
#f2py_source       = step.f95
MODULENAME		  = tracmass


runtracmass : $(objects)
	$(FF)  $(MYI_FLAGS) -o runtracmass $(objects) $(LNK_FLAGS) $(MYL_FLAGS)

%.o : %.f95
	$(FF) $(FF_FLAGS) $(ORM_FLAGS) $(PROJECT_FLAG) $(CASE_FLAG) $(ARG_FLAGS)  $< -o $@

$(objects) :

f2py : $(objects)
	$(FF) -E $(ARG_FLAGS) -x f95-cpp-input step.f95 -o outdir/step.f95
	$(F2PY) $(objects) -c $(f2py_source) -m $(MODULENAME)

.PHONY : clean
clean :
	-rm $(objects)
	-rm $(MODULENAME).so
