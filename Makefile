IDIR=include
CXXFLAGS=-I$(IDIR) -Wall -O3 -std=c++14

_DEPS = utils.hh vec.hh matrix.hh
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

ODIR=obj
_OBJ = utils.o vec.o matrix.o simplex.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

SDIR=src


all: simplex.exe

simplex.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

$(ODIR)/%.o: $(SDIR)/%.cc $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o simplex.exe