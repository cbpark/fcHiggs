SRCDIR   := src
BINDIR   := bin
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++14 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LDFLAGS  := -O2 $(LDFLAGS)
LIBS     :=
AR       := ar crs
MKDIR    := mkdir -p
RM       := rm -f

# Targets
EXE    := $(BINDIR)/pph
EXESRC := $(patsubst $(BINDIR)/%,$(SRCDIR)/%.cc,$(EXE))
EXEOBJ := $(EXESRC:.cc=.o)
LIBSRC := $(filter-out $(EXESRC),$(wildcard $(SRCDIR)/*.cc))
LIBOBJ := $(LIBSRC:.cc=.o)

.PHONY: all build clean

all: $(EXE)

$(BINDIR)/%: $(SRCDIR)/%.o build $(LIBOBJ)
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBOBJ) $(LIBS)

$(LIBOBJ): CXXFLAGS += -fPIC

build:
	$(MKDIR) $(BINDIR)

clean::
	$(RM) $(EXEOBJ) $(LIBOBJ) $(EXE)
	$(RM) -r $(BINDIR)
