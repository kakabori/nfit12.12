include ./make.inc

all:
	(cd .; \
	cd ./lib; make -f make.blas; \
	make -f make.subplex; \
	make -f make.misc; \
	make -f make.minimization; \
	make -f make.stat; \
	make)

test: test.cpp
	$(CXX) test.cpp -lmydll

clean:
	rm -f lib/*.o obj/*.o core *~ a.out *.so;


