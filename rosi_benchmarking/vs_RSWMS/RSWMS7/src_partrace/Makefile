include ../Makefile.defs

PREFIX = ../

sources_cp = partrace.cpp elements.cpp elementsZ.cpp elementsXYZ.cpp \
	  elementsRectilinear.cpp particlelist.cpp particle.cpp io.cpp \
	  particleclassSorp.cpp decay.cpp uptake.cpp interact.cpp parallel.cpp \
	  particleclass.cpp random.cpp runtime.cpp event.cpp \
	  distdata.cpp partrace_interface.cpp
objects_cp = $(addprefix $(PREFIX)$(BUILDDIR)/, ${sources_cp:.cpp=.o})

.SUFFIXES: .o .cpp

all: $(PREFIX)$(BUILDDIR) $(PREFIX)$(LIBPARTRACE)

$(PREFIX)$(BUILDDIR):
	mkdir -p $@

$(PREFIX)$(BUILDDIR)/%.o: %.cpp
	@echo "compiling cpp file" $<
	cd $(PREFIX)$(BUILDDIR) ; $(CXX) $(CXXFLAGS) -c ../$(SRC_PARTRACE)/$< -o $(@F)

$(PREFIX)$(LIBPARTRACE): $(objects_cp)
	cd $(PREFIX)$(BUILDDIR) ; rm -f $(@F)
	cd $(PREFIX)$(BUILDDIR) ; ar -rs $(@F) $(^F)

clean:
	rm -rf $(PREFIX)$(BUILDDIR) 


partrace.o: partrace.hpp runtime.hpp particlelist.hpp particleclass.hpp
partrace.o: particle.hpp definitions.hpp elements.hpp dispersion.hpp
partrace.o: distdata.hpp parallel.hpp random.hpp interact.hpp
elements.o: elements.hpp definitions.hpp dispersion.hpp distdata.hpp
elements.o: particleclass.hpp particle.hpp random.hpp partrace.hpp
elements.o: runtime.hpp particlelist.hpp parallel.hpp event.hpp
elementsZ.o: elements.hpp definitions.hpp dispersion.hpp distdata.hpp
elementsZ.o: particleclass.hpp particle.hpp random.hpp partrace.hpp
elementsZ.o: runtime.hpp particlelist.hpp parallel.hpp
elementsXYZ.o: elements.hpp definitions.hpp dispersion.hpp distdata.hpp
elementsXYZ.o: particleclass.hpp particle.hpp random.hpp partrace.hpp
elementsXYZ.o: runtime.hpp particlelist.hpp
elementsRectilinear.o: elements.hpp definitions.hpp dispersion.hpp
elementsRectilinear.o: distdata.hpp particleclass.hpp particle.hpp random.hpp
elementsRectilinear.o: partrace.hpp runtime.hpp particlelist.hpp parallel.hpp
particlelist.o: particlelist.hpp particleclass.hpp particle.hpp
particlelist.o: definitions.hpp elements.hpp dispersion.hpp distdata.hpp
particlelist.o: runtime.hpp partrace.hpp parallel.hpp decay.hpp uptake.hpp
particle.o: particle.hpp definitions.hpp
io.o: io.hpp definitions.hpp partrace.hpp runtime.hpp particlelist.hpp
io.o: particleclass.hpp particle.hpp elements.hpp dispersion.hpp distdata.hpp
io.o: parallel.hpp interact.hpp event.hpp random.hpp
particleclassSorp.o: particleclass.hpp particle.hpp definitions.hpp
particleclassSorp.o: elements.hpp dispersion.hpp distdata.hpp random.hpp
particleclassSorp.o: event.hpp partrace.hpp runtime.hpp particlelist.hpp
particleclassSorp.o: parallel.hpp
decay.o: decay.hpp runtime.hpp elements.hpp definitions.hpp dispersion.hpp
decay.o: distdata.hpp partrace.hpp particlelist.hpp particleclass.hpp
decay.o: particle.hpp random.hpp
uptake.o: uptake.hpp runtime.hpp elements.hpp definitions.hpp dispersion.hpp
uptake.o: distdata.hpp partrace.hpp particlelist.hpp particleclass.hpp
uptake.o: particle.hpp random.hpp
interact.o: interact.hpp elements.hpp definitions.hpp dispersion.hpp
interact.o: distdata.hpp random.hpp partrace.hpp runtime.hpp particlelist.hpp
interact.o: particleclass.hpp particle.hpp parallel.hpp
parallel.o: parallel.hpp 
particleclass.o: particleclass.hpp particle.hpp definitions.hpp elements.hpp
particleclass.o: dispersion.hpp distdata.hpp random.hpp event.hpp
particleclass.o: partrace.hpp runtime.hpp particlelist.hpp parallel.hpp
particleclass.o: decay.hpp uptake.hpp
random.o: random.hpp
runtime.o: runtime.hpp event.hpp partrace.hpp particlelist.hpp
runtime.o: particleclass.hpp particle.hpp definitions.hpp elements.hpp
runtime.o: dispersion.hpp distdata.hpp parallel.hpp
event.o: event.hpp
distdata.o: distdata.hpp parallel.hpp
partrace_interface.o: partrace_interface.hpp parallel.hpp distdata.hpp
