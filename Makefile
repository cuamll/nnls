SHELL   = /bin/sh
CC      = gfortran
FLAGS   = -std=f2018 -ffree-form -Wall -Werror -pedantic
LDFLAGS = -llapack -lblas
DFLAGS  = -g
RFLAGS  = -O2
SOURCES = nnls.f

# if SHARED = 1 compile nnls module as a shared library
# if SHARED = 0 compile standalone program to test/debug
SHARED = 0
ifeq ($(SHARED), 1)
	CFLAGS += -fPIC
	LDFLAGS += -shared
	TARGET = libnnls.so
else
	SOURCES += main.f
	TARGET = test
endif
OBJECTS = $(SOURCES:.f=.o)

%.o: %.f
	$(CC) $(FLAGS) $(CFLAGS) $(RFLAGS) -c $<

.PHONY: clean

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(FLAGS) $(CFLAGS) $(DFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)
