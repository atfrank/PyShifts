include ./makefile.defs

all: libmoletools.$(LIBEXT) programs

libmoletools.$(LIBEXT):
	$(UNAME)
	$(MAKE) -C lib 

programs:
	$(MAKE) -C src

clean:
	$(MAKE) -C lib clean
	$(MAKE) -C src clean
	rm -rf *~
	rm -rf *stackdump
