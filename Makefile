CXX=g++
CFLAGS=-I.
LDFLAGS=-L. -ltinyxml2 -lstdc++ -lm
DEPS = tinyxml2.h
OBJ = sketcher.o 

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

sketcher: $(OBJ)
	gcc -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f *.o sketcher
