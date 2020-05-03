demo : example.cpp crn.o
	g++ example.cpp crn.o -o demo
crn.o : crn.cpp crn.h
	g++ -c crn.cpp
clean:
	rm demo crn.o
