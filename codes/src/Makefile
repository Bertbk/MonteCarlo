#---------------------------------------------------------------------
BINDIR  = ../bin
LDIR    = 
IDIR    = 
LIBS    = 

# Directives
#OMP=-DNO_OPENMP
#OMP=-D_OPENMP
OBJ     = Message.o Point.o Main_Monte.o
#---------------------------------------------------------------------

all: clean checkdir Main_Monte

clean:
	rm -f $(OBJ)

checkdir:
	@(test -d "$(BINDIR)"  || mkdir -p $(BINDIR) )

Main_Monte : $(OBJ)
	g++ -o $(BINDIR)/MonteCarlo $(OBJ) $(IDIR) $(LDIR) $(LIBS)

Main_Monte.o : Main_Monte.cpp
	g++ -o $@ -c Main_Monte.cpp $(IDIR)

Point.o : Point.cpp
	g++ -o $@ -c Point.cpp $(IDIR)

Message.o: Message.cpp
	g++ -o $@ -c Message.cpp $(IDIR)
