SHELL=	/bin/sh

include $(ESMFMKFILE)

CMD= global_chgres.exe

OBJS = chgres.o \
       atmosphere.o \
       input_data.o \
       surface.o \
       model_grid.o \
       program_setup.o \
       write_data.o \
       search_util.o \
       static_data.o

$(CMD):	$(OBJS)
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -o $(CMD) $(OBJS) $(W3NCO_LIBd) $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKRPATHS) $(ESMF_F90ESMFLINKLIBS) -g -traceback
	mv $(CMD) ../exec

model_grid.o:  program_setup.o model_grid.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c model_grid.F90

program_setup.o:  program_setup.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c program_setup.f90

atmosphere.o:  program_setup.o model_grid.o input_data.o atmosphere.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c atmosphere.F90

chgres.o: atmosphere.o model_grid.o program_setup.o surface.o chgres.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c chgres.F90

write_data.o: input_data.o atmosphere.o  model_grid.o program_setup.o surface.o static_data.o write_data.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c write_data.F90

input_data.o: program_setup.o model_grid.o input_data.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c input_data.F90

surface.o: search_util.o model_grid.o input_data.o program_setup.o static_data.o surface.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c surface.F90

search_util.o: search_util.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c search_util.f90

static_data.o: model_grid.o program_setup.o static_data.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c static_data.F90

clean:
	rm -f *.o *.mod ${CMD} ../exec/${CMD}
