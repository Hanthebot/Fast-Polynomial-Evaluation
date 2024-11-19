CC     = g++
CFLAGS = -std=c++20 -Wall -O3 -o
LFLAGS = -std=c++20 -Wall -O3 -c

span_test: span_test.cpp nd_vector.o
	$(CC) $(CFLAGS) span_test span_test.cpp nd_vector.o

nd_vector_test: nd_vector_test.cpp nd_vector.o
	$(CC) $(CFLAGS) nd_vector_test nd_vector_test.cpp nd_vector.o

nd_vector_set_test: nd_vector_set_test.cpp nd_vector.o
	$(CC) $(CFLAGS) nd_vector_set_test nd_vector_set_test.cpp nd_vector.o

nd_vector.o: nd_vector.cpp nd_vector.h
	$(CC) $(LFLAGS) nd_vector.cpp

clean:
	del *.o

