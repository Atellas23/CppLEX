IDIR=include
CXXFLAGS=-I$(IDIR) -Wall -O3 -std=c++14 -D_GLIBCXX_DEBUG

_DEPS = utils.hh vec.hh matrix.hh
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

ODIR=obj
_OBJ = utils.o vec.o matrix.o simplex.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

TDIR=test
_TOBJ = utils.o vec.o matrix.o 
TOBJ = $(patsubst %,$(ODIR)/%,$(_TOBJ))

SDIR=src


all: simplex.exe

simplex.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

$(ODIR)/%.o: $(SDIR)/%.cc $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

.PHONY: clean test $(TDIR)/test_%.exe

$(TDIR)/test_%.exe: $(TDIR)/test_%.cc $(TOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

all-tests := $(addsuffix .exe, $(basename $(wildcard $(TDIR)/test_*.cc)))

test: $(all-tests)
	$(patsubst %, %;, $^)

clean:
	rm -f $(ODIR)/*.o simplex.exe $(TDIR)/*.exe