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
OBJSRC := $(filter-out $(EXESRC),$(wildcard $(SRCDIR)/*.cc))
OBJ    := $(OBJSRC:.cc=.o)

.PHONY: all build clean

all: $(EXE)

$(BINDIR)/%: $(SRCDIR)/%.o build $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $< $(OBJ) $(LIBS)

$(LIBOBJ): CXXFLAGS += -fPIC

build:
	$(MKDIR) $(BINDIR)

clean::
	$(RM) $(EXEOBJ) $(OBJ) $(EXE)
	$(RM) -r $(BINDIR)
