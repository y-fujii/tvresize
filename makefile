tvresize: makefile tvresize.cpp
	clang++ -g -pedantic -Wall -Wextra -O3 -o tvresize tvresize.cpp $$(libpng-config --cflags --ldflags)
