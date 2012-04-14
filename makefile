tvresize: makefile tvresize.cpp
	$(CXX) -pedantic -Wall -Wextra -O3 -o tvresize tvresize.cpp $$(libpng-config --cflags --ldflags) -I/usr/local/include
