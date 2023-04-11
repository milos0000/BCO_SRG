program : bco_srg.o jacobi_eigenvalue.o normal.o rand_el_unordered_set.o spectral_distance.o write_to_file.o 
	g++ bco_srg.o jacobi_eigenvalue.o normal.o rand_el_unordered_set.o spectral_distance.o write_to_file.o -o program

write_to_file.o : src/write_to_file.cpp 
	g++ -c src/write_to_file.cpp

bco_srg.o : src/bco_srg.cpp include/bco_srg.hpp
	g++ -c src/bco_srg.cpp

jacobi_eigenvalue.o : src/jacobi_eigenvalue.cpp include/jacobi_eigenvalue.hpp
	g++ -c src/jacobi_eigenvalue.cpp 

normal.o : src/normal.cpp
	g++ -c src/normal.cpp 

rand_el_unordered_set.o : src/rand_el_unordered_set.cpp 
	g++ -c src/rand_el_unordered_set.cpp 

spectral_distance.o : src/spectral_distance.cpp 
	g++ -c src/spectral_distance.cpp 

clean: 
	rm -f *.o *.gch program