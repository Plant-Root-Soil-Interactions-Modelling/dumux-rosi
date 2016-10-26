include ../Makefile.defs

PREFIX = ../

sources_f = typedef.f90 sparsematrix.f90 Modules.f90 SoluteRoot.f90 Input.f90 Orthofem.f90 Water.f90 Doussan.f90 Solute.f90  Sink.f90 Environmental.f90 Output.f90 RootTyp.f90 RootGrowth.f90 PlantGrowth.f90 SolveRootDirect.f90 SolutePartrace.f90 RSWMS34.f90 
objects_f = $(addprefix $(PREFIX)$(BUILDDIR)/, ${sources_f:.f90=.o})

sources_c = RootTyp1boite.c RootTyp1_interface.c 
objects_c = $(addprefix $(PREFIX)$(BUILDDIR)/, ${sources_c:.c=.o})

.SUFFIXES: .f90 .o .c

all: $(PREFIX)$(BUILDDIR) $(PREFIX)$(LIBRSWMS)

$(PREFIX)$(BUILDDIR):
	mkdir -p $@

$(PREFIX)$(BUILDDIR)/%.o: %.f90
	@echo "compiling f90 file" $<
	cd $(PREFIX)$(BUILDDIR); $(F90C) $(F90FLAGS) -I../$(SRC_RSWMS)/ -c ../$(SRC_RSWMS)/$< -o $(@F)

$(PREFIX)$(BUILDDIR)/%.o: %.c
	@echo "compiling c file" $<
	cd $(PREFIX)$(BUILDDIR) ; $(CC) $(CFLAGS) -c ../$(SRC_RSWMS)/$< -o $(@F)

$(PREFIX)$(LIBRSWMS):  $(objects_f) $(objects_c)
	cd $(PREFIX)$(BUILDDIR) ; rm -f $(@F)
	cd $(PREFIX)$(BUILDDIR) ; ar -rs $(@F) $(^F)

clean:
	rm -rf $(PREFIX)$(BUILDDIR) 


RootTyp1boite.o: RootTyp1boite.h
RootTyp1_interface.o: RootTyp1boite.h RootTyp1_interface.h
