HOMEL=$(PWD)
TOP=../bin
## binary for install  
# DEBUG_OPTIONS= -g -pg
DEBUG_OPTIONS= -g -pg
BIN=$(HOMEL)/../bin
PERLDIR=$(HOMEL)/perlsrc

NLIB=$(HOMEL)/admixlib/nicklib.a
IDIR=$(HOMEL)/admixinclude
VPATH=.:nicksrc
DEBUG_OPTIONS= -g -pg

CFLAGS= -c -g -p -I$(IDIR) -Wimplicit

OBJ=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o
TWTAB=\"$(HOMEL)/admixtables/twtable\"

statsubs.o:     tables	nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o nicksrc/statsubs.c


X1=neanderthal_dating
X1O=neanderthal_dating.o  mcio.o  admutils.o  egsubs.o h2d.o  exclude.o badpairs.o

X2=grabpars  
X2O=grabpars.o 

 
IDIR=$(HOMEL)/admixinclude  
##ISTR= -I$(HOME)/include  -I/opt/nag/libC/include
##NAGL= -L/opt/nag/libC/lib -lnagc_nag -lm -pthread

##nagqp.o:		nagqp.c
#	$(CC)  $(CFLAGS) -I/opt/nag/libC/include -o nagqp.o nagqp.c


$(X1): nicklib	$(X1O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(X1) $(X1O) $(NLIB)

$(X3): nicklib  $(X3O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(X3) $(X3O) $(NLIB)

$(X2): nicklib	$(X2O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(X2) $(X2O) $(NLIB)


PROGS=neanderthal_dating
all: nicklib $(PROGS)

install:	all	tables 
	mkdir -p $(TOP)
	cp $(PROGS) $(TOP)

uninstall:
	rm -f $(NLIB) tables
	rm -f $(TOP)/*

libnick.a:	dirs tables  $(OBJ)
	ar -r libnick.a $(OBJ)
	ranlib libnick.a

nicklib:	dirs  tables  libnick.a 
	cp libnick.a  $(NLIB)

tables:    
	echo "tables made"  > tables
	cp twtable  $(HOMEL)/admixtables
	
dirs:	
	mkdir -p  $(HOMEL)/admixlib
	mkdir -p  $(HOMEL)/admixtables
	mkdir -p  $(HOMEL)/admixinclude
	mkdir -p  $(BIN)
	cp  *.h  $(IDIR)
	cp  nicksrc/*.h  $(IDIR)

clean: 
	rm -f *.o 
	rm -f *junk*
	rm -f core
	rm -f libnick.a
	rm -f $(PROGS)
	rm -f nicksrc/*.o

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/admixlib 
	rm -rf $(HOMEL)/admixtables 
	rm -rf $(HOMEL)/admixinclude 

NAGL= -L/opt/nag/libC/lib -lnagc_nag -lm -pthread


