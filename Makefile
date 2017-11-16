PKGNAME  := fcHiggs
SRCDIR   := src
BINDIR   := bin
LIBDIR   := lib
CXXFLAGS := -g -O3 -m64 -march=native -Wall -Wextra -std=c++14 -pedantic -I$(SRCDIR) \
	$(CXXFLAGS)
LDFLAGS  := -O3 -m64 $(LDFLAGS)
LIBS     :=
AR       := ar crs
MKDIR    := mkdir -p
RM       := rm -f

# Targets
EXE    := $(BINDIR)/pph $(BINDIR)/pphb_neutral \
	$(BINDIR)/ppht_charged $(BINDIR)/pphb_charged \
	$(BINDIR)/hdecay
EXESRC := $(patsubst $(BINDIR)/%,$(SRCDIR)/%.cc,$(EXE))
EXEOBJ := $(EXESRC:.cc=.o)
LIB    := $(LIBDIR)/lib$(PKGNAME).a
LIBSRC := $(filter-out $(EXESRC),$(wildcard $(SRCDIR)/*.cc))
LIBOBJ := $(LIBSRC:.cc=.o)

# LHAPDF (http://lhapdf.hepforge.org/)
CXXFLAGS += -I$(shell lhapdf-config --incdir)
LDFLAGS  += -Wl,-rpath,$(shell lhapdf-config --libdir)
LIBS     += -L$(shell lhapdf-config --libdir) -lLHAPDF

.PHONY: all build clean

all: $(EXE)

$(BINDIR)/%: $(SRCDIR)/%.o build $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $< -L$(LIBDIR) -l$(PKGNAME) $(LIBS)

$(LIB): CXXFLAGS += -fPIC
$(LIB): $(LIBOBJ)
	$(AR) $@ $^
	ranlib $@

build:
	$(MKDIR) $(LIBDIR)
	$(MKDIR) $(BINDIR)

clean::
	$(RM) $(EXEOBJ) $(LIBOBJ)
	$(RM) $(EXE) $(LIB)
	$(RM) -r $(BINDIR) $(LIBDIR)
