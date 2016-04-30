VAR=`cat input.txt`
all:
	cd CEAexec; ./FCEA2 < ../input.txt
	cp CEAexec/$(VAR).plt Plot_CEA/
	cd Plot_CEA; make
	rm Plot_CEA/$(VAR).plt
	#vi CEAexec/$(VAR).out
