CC = g++
CFLAGS = -std=c++11

all:clean hungarian_lib cAccel_lib

clean:
	-rm CAccel*.so libhungarian.so

hungarian_lib:
	g++ Hungarian.cpp -fPIC -shared -o libhungarian.so

cAccel_lib:
	python setup.py build_ext --inplace
