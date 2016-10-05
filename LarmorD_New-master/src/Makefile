#Sean M. Law
#Aaron T. Frank

include ../makefile.defs

OBJS = libmoletools.$(LIBEXT)

all:	$(BINDIR) libmoletools.$(LIBEXT) larmord larmord_extractor

clean:
	rm -rf *.o
	rm -rf *stackdump
	rm -rf ../bin/*
	rm -rf *~

$(BINDIR):
	mkdir -p $(BINDIR)

libmoletools.$(LIBEXT):
	$(MAKE) -C ../lib libmoletools.$(LIBEXT)

larmord: larmord.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LIBFLAGS) $(INCFLAGS) -o $(BINDIR)/$@ $< $(LINKFLAGS)

larmord_extractor: larmord_extractor.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LIBFLAGS) $(INCFLAGS) -o $(BINDIR)/$@ $< $(LINKFLAGS)

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(INCFLAGS) -c $< 
