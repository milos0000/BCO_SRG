program : bco_srg.o jacobi_eigenvalue.o normal.o rand_el_unordered_set.o spectral_distance.o write_to_file.o 
	g++ bco_srg.o jacobi_eigenvalue.o normal.o rand_el_unordered_set.o spectral_distance.o write_to_file.o -o program

write_to_file.o : write_to_file.cpp 
	g++ -c write_to_file.cpp

bco_srg.o : bco_srg.cpp bco_srg.hpp
	g++ -c bco_srg.cpp

jacobi_eigenvalue.o : jacobi_eigenvalue.cpp jacobi_eigenvalue.hpp
	g++ -c jacobi_eigenvalue.cpp 

normal.o : normal.cpp
	g++ -c normal.cpp 

rand_el_unordered_set.o : rand_el_unordered_set.cpp 
	g++ -c rand_el_unordered_set.cpp 

spectral_distance.o : spectral_distance.cpp 
	g++ -c spectral_distance.cpp 

clean: 
	rm -f *.o *.gch program