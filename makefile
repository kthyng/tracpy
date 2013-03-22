###############################################################
### TEMPLATE FOR Mac OS X 10.7 with gfortran from fink ###
###############################################################
# PROJECT	          = tes
# # possible PROJECTS i.e. GCMs: rco, baltix, occ, orc, sim, for, ifs, tes, gomoos
# CASE              = $(PROJECT)

#F95COMPILER        = "g95"
F95COMPILER        = "gfortran"

#================================================================

# PROJMAKE           := $(wildcard projects/$(PROJECT)/Makefile)
# ifneq ($(strip $(PROJMAKE)),)
# 	include projects/$(PROJECT)/Makefile
# endif

# PROJECT_FLAG      = -DPROJECT_NAME=\'$(PROJECT)\'
# CASE_FLAG         = -DCASE_NAME=\'$(CASE)\'
   fl01   = -Dtimestep       # Analytical time scheme used to solve the differential Eqs.
#  fl01   = -Dtimeanalyt    # Time steps with analytical stationary scheme 
#   fl02   = -Dregulardt      # Regular time steps to be used with -Dtimestep
#------------------------------------------------------------------------
#  fl03   = -Dtempsalt       # Include temp and salt
#  fl04   = -Dtextwrite      # Write results to textfile
#  fl05   = -Dbinwrite       # Write results to binaryfile
#  fl06   = -Dmysqlwrite     # Write results to mysql
#------------------------------------------------------------------------
#  fl07  = -Dstreamxy       # Calculates the barotropic stream function.
#  fl08  = -Dstreamr        #    vertical stream function, z=density
#  fl09  = -Dstreamts       #    vertical stream function, z=???
#  fl10  = -Dstreamv        #    vertical stream function, z=depth
#  fl11  = -Drerun          # Stores the Lagrangian stream functions as a
                            # function of the end positions that has been
                            # calculated in an identical previous run.    
#   fl12  = -Dinitxyt        # Start trajectories at given positions and times
#------------------------------------------------------------------------
  fl13  = -Dtwodim         # Turn off vertical velocities.
#  fl14  = -Dfull_wflux     # Use a full 3D wflux field.
#  fl15  = -Dexplicit_w     # Use a given vertical velocity.
#------------------------------------------------------------------------
#   fl16  = -Dvarbottombox   # Variable bottom box to match actual depth
#   fl17  = -Dfreesurface    # Variable bottom box to match actual depth
#   fl18  = -Dzvec1D         # Cell depths defined as vector (for z-coord?)
#  fl18  = -Dzgrid3D        # Cell depths defined as 3D grid (for sigma)
#  fl18  = -Dzgrid3Dt       # Cell depths 3D and time interp. (for atm)
#------------------------------------------------------------------------
#  fl20  = -Dselect         # Select only one trajectory (for debugging)
#  fl21  = -Dtracer         # Stores a simulated tracer
#  fl22  = -Dsediment       # Sediment code developed for RCO
#------------------------------------------------------------------------
#  fl23  = -Dturb           # Adds subgrid turbulent velocities 
#  fl24  = -Ddiffusion      # Adds a diffusion on trajectory
#  fl25  = -Danisodiffusion # Adds an anisotropic diffusion on trajectory
#  fl26  = -Dcoarse         # Adds a diffusion on trajectory
#========================================================================

ARG_FLAGS=  \
$(fl01)$(fl02)$(fl03)$(fl04)$(fl05)$(fl06)$(fl07)$(fl08)$(fl09)$(fl10)\
$(fl11)$(fl12)$(fl13)$(fl14)$(fl15)$(fl16)$(fl17)$(fl18)$(fl19)$(fl20)\
$(fl21)$(fl22)$(fl23)$(fl24)$(fl25)$(fl26)$(fl27)$(fl28)$(fl29)$(fl30)\

# ARG_FLAGS         = -D

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
	F2PY = f2py-2.7
else ifeq ($(findstring hafen,$(COMPUTER)),hafen)
	F2PY = f2py2.7
else
	F2PY = f2py
endif

objects           = pos.o cross.o calc_dxyz.o calc_time.o loop_pos.o
f2py_source       = step.f95
MODULENAME		  = tracmass
                    
# objects           = modules.o seed.o  savepsi.o loop_pos.o init_seed.o \
#                     sw_stat.o loop.o time_subs.o getfile.o \
#                     vertvel.o coord.o cross.o init_par.o \
#                     interp.o interp2.o pos.o \
#                     sw_seck.o sw_pres.o sw_dens0.o \
#                     writepsi.o writetracer.o turb.o main.o \
# 		            setupgrid.o readfield.o diffusion.o \
#objwdir = $(patsubst %,tmp/%,$(objects))
#jacket.o

runtracmass : $(objects)
	$(FF)  $(MYI_FLAGS) -o runtracmass $(objects) $(LNK_FLAGS) $(MYL_FLAGS)

%.o : %.f95
	$(FF) $(FF_FLAGS) $(ORM_FLAGS) $(PROJECT_FLAG) $(CASE_FLAG) $(ARG_FLAGS)  $< -o $@

$(objects) : 

f2py : $(objects)
	$(F2PY) $(objects) -c $(f2py_source) -m $(MODULENAME)
#readfield.o:  projects/$(PROJECT)/readfield.f95
#	$(FF) $(FF_FLAGS) $(ORM_FLAGS) projects/$(PROJECT)/readfield.f95 -o tmp/$@

#stat.o:  $(PROJECT)/stat.f95
#	$(FF) $(FF_FLAGS) $(ORM_FLAGS) $(PROJECT)/stat.f95

# jacket.o : ../mysql/jacket.c
# 	$(CC)  -c ../mysql/jacket.c

#main.o : main.f95 
#	$(FF) $(FF_FLAGS) $(ORM_FLAGS) main.f95

.PHONY : clean
clean :
	-rm $(objects)
	-rm $(MODULENAME).so
