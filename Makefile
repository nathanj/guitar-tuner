CC     = gcc `pkg-config --cflags sdl`
CFLAGS = -g
LDLIBS = -lfftw3 -lasound `pkg-config --libs sdl`

dft: dft.o

