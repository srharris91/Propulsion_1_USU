VAR=`cat input.txt`
all: CEA plot

CEA:
	@cd CEAexec; ./FCEA2 < ../input.txt
plot:
	@cd Plot_CEA; make
