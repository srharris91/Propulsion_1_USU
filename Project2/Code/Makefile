all: cpp python
cpp:
	cd cpp_stuff; make

python:
	cp cpp_stuff/answer.txt python_stuff/
	cd python_stuff; python plot_stuff.py

clean:
	cd cpp_stuff; make clean

