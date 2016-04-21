comp = g++
flags = -g -Wall

OBJS = X.o Project1.o
EXEC = ../Project1.exec

all: exec run

exec: $(OBJS)
	@echo ''
	@echo 'Linking $(EXEC) from $(OBJS)'
	@echo ''
	@$(comp) $(flags) $(OBJS) -o $(EXEC)

run:
	./$(EXEC)

clean: 
	@rm -f *.o *.exec *.mod *.dat *.txt $(EXEC)






.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	@echo 'Compiling [$@] from [$<] using $(comp)'
	@$(comp) $(flags) $*.cpp -c -o $*.o

