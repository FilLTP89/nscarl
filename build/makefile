#### PATH TO HDF5 AND MPI LIBRARY AND INCLUDES - to be modified by user (to be completed)

LIBHDF5     = -L ~/LOCAL/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -ldl -lm
INCLUDEHDF5 = -I ~/LOCAL/include

LIBMPI      = -L/gpfs/opt/compilers/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/lib -lmpi -lmpi_dbg -lmpi_mt -lmpigf -lmpi_ilp64
INCLUDEMPI  = -I/gpfs/opt/compilers/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/include

LIBMPI      = -L/gpfs/opt/compilers/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/lib -lmpi -lmpi_dbg -lmpi_mt -lmpigf -lmpi_ilp64
INCLUDEMPI  = -I/gpfs/opt/compilers/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/include

#LIBFFTW     = -L ~/LOCAL/lib -lfftw3_mpi -lfftw3 -lfftw3_threads
#INCLUDEFFTW = -I ~/LOCAL/include

EXEC  = merge3dFields.exe
#EXEC2 = statistics.exe
#LIB1 = libRF.a
FC   = mpiifort
FFLAGS = -g -O1 -cpp #-DMAKELOG

#### FROM THIS POINT ON THE USER SHOULD NOT MODIFY
RM := rm -rf

#Specifying the sources "f90" and "f"
SRCS = $(wildcard *.f90 *.f) 

#Create a ".f90" for each source
#OBJS = $(SRCS:.f90=.o) #SYNTAX NOT WORKING

RF_OBJ += ./merge3dFields.o

LIB_OBJ += ./sem_hdf5.o
#./spectra_RF.o \
#./writeResultFile_RF.o \
#./mesh_RF.o \
#./type_MESH.o \
#./type_RF.o \
#./topography_RF.o \
#./readUNV_RF.o \
#./localization_RF.o \
#./type_inputRF.o \
#./sample_RF.o \
#./calls_RF.o \
#./ranlib.o \
#./rnglib.o \
#./ipmpar.o \
#./spmpar.o \
#./cumnor.o \
#./devlpl.o \
#./stvaln.o \
#./dinvnr.o \
#./cdfnor.o

#STAT_OBJ += ./main_Stat.o \
#./statistics_RF.o \
#./type_STAT.o

COMM_OBJ += ./constants_RF.o \
./constants.o \
./hdf5_RF.o

#./type_TEST.o \
#./test_func_RF.o \

#LIBS = $(LIBHDF5) $(LIBMPI) $(LIBFFTW)
#INCLUDE = $(INCLUDEHDF5) $(INCLUDEMPI) $(INCLUDEFFTW)
LIBS = $(LIBHDF5) $(LIBMPI) 
INCLUDE = $(INCLUDEHDF5) $(INCLUDEMPI) 

#Dependencies
merge3dFields.o 	 : constants_RF.o \
				hdf5_RF.o			  	

sem_hdf5.o           : constants_RF.o
# Making all the ".o" from the ".f90"
%.o: ../SRC_EXE/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '
		
%.o: ../SRC_LIB/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '

# All Target
all: merge3dFields
#all: randomField

# Tool invocations
merge3dFields: $(RF_OBJ) $(COMM_OBJ)
		@echo 'Building target: $@'
		@echo 'Invoking: Fortran Linker'
		$(FC) -o $(EXEC) $(FFLAGS) $(RF_OBJ) $(COMM_OBJ) $(INCLUDE) $(LIBS)
		@echo 'Finished building target: $@'
		@echo ' '	
	
#statistics: $(STAT_OBJ) $(COMM_OBJ)
#		@echo 'Building target: $@'
#		@echo 'Invoking: Fortran Linker'
#		$(FC) -o $(EXEC2) $(FFLAGS) $(STAT_OBJ) $(COMM_OBJ) $(INCLUDE) $(LIBS)
#		@echo 'Finished building target: $@'
#		@echo ' '

#RFlib: $(LIB_OBJ) $(COMM_OBJ)
#		@echo  'Building target: $@'
#		@echo  'Invoking: Fortran Linker'
#		ar rc  $(LIB1) $(LIB_OBJ) $(COMM_OBJ) $(LIBS)
#		@echo  'Finished building target: $@'
#		@echo ' '
#
#RFexe: $(RF_OBJ) $(COMM_OBJ)
#		@echo 'Building target: $@'
#		@echo 'Invoking: Fortran Linker'
#		$(FC) -o RFexe.exe $(FFLAGS) /home/carvalhol/Projects/RANDOM_FIELD/randomField/SRC_EXE/main_RandomField.f90  $(INCLUDE) $(LIBS) libRF.a
#		@echo 'Finished building target: $@'
#		@echo ' '	
	

# Other Targets
clean:
		-$(RM) $(LIB1) $(EXEC) $(RF_OBJ) $(COMM_OBJ) *.mod
		-@echo 'Finished excluding'

.PHONY: all clean dependents
.SECONDARY:
