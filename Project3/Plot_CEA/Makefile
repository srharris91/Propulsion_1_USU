file=../input.txt
VAR=`cat ${file}`
all:
	@cp ../CEAexec/$(VAR).plt ./
	@python plot_stuff.py `cat ${file}`.plt
	@rm $(VAR).plt
