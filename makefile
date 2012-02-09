tvresize: makefile tvresize.cpp
	g++46 -g -pedantic -Wall -Wextra -O3 -o tvresize tvresize.cpp $$(libpng-config --cflags --ldflags) -I/usr/local/include
