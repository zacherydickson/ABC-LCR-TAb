override CCFLAGS += -std=c++14
CC = g++ $(CCFLAGS)
HOST := $(shell hostname)
SDIR = src
ODIR = obj
HODIR = $(ODIR)/$(HOST)
BDIR = bin
HBDIR = $(BDIR)/$(HOST)
PROG = run_abc
HEADERS = $(wildcard $(SDIR)/*.hpp)
OBJECTS = $(addsuffix .o,$(basename $(subst $(SDIR),$(HODIR),$(HEADERS))))
STANDALONES = $(wildcard $(SDIR)/*.h)
BOOSTLIB=/home/zac/applications/boost/current

all: $(HBDIR)/$(PROG)

$(HBDIR)/$(PROG): $(SDIR)/$(PROG).cpp $(STANDALONES) $(OBJECTS)
	@mkdir -p $(HBDIR)
	$(CC) -pthread $(HODIR)/* $(SDIR)/$(PROG).cpp -o $(HBDIR)/$(PROG)

$(HODIR)/Distributions.o : $(SDIR)/Distributions.cpp $(SDIR)/Distributions.hpp
	@mkdir -p $(HODIR)
	$(CC) -c $(SDIR)/Distributions.cpp -o $(HODIR)/Distributions.o -I$(BOOSTLIB)

$(HODIR)/%.o : $(SDIR)/%.cpp $(SDIR)/%.hpp
	@mkdir -p $(HODIR)
	$(CC) -c -o $(HODIR)/$*.o $< 

.PHONY: clean
clean:
	@rm -rf $(HBDIR) $(HODIR)

.PHONY: distclean
distclean:
	@rm -rf $(BDIR) $(ODIR)
	@find . -maxdepth 1 -type f ! -name Makefile ! -name .gitignore -exec rm -f {} +
