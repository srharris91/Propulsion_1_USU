comp = g++
flags = -g -Wall -O3

OBJS = Prandtl_Meyer.o CylindricalPort.o RP.o Solver.o Project2.o
EXEC = ../Project2.exec

all: exec run

exec: $(OBJS)
	@echo ''
	@echo 'Linking $(EXEC) from $(OBJS)'
	@echo ''
	@$(comp) $(flags) $(OBJS) -o $(EXEC)

run:
	./$(EXEC) > answer.txt

clean: 
	@rm -f *.o *.exec *.mod *.dat *.txt $(EXEC)






.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	@echo 'Compiling [$@] from [$<] using $(comp)'
	@$(comp) $(flags) $*.cpp -c -o $*.o

