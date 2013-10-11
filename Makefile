all: BetheGams_mex.mexmaci64 makePotential_mex.mexmaci64 makePotential_mex.mexmaci64 submodularMAP_mex.mexmaci64

submodularMAP_mex.mexmaci64: submodularMAP_mex.cpp
	mex -largeArrayDims -g submodularMAP_mex.cpp MinSum.cpp BetheApprox.cpp graph.cpp maxflow.cpp -o submodularMAP_mex.mexmaci64

BetheGams_mex.mexmaci64: BetheGams_mex.cpp
	mex -largeArrayDims -O BetheGams_mex.cpp MinSum.cpp BetheApprox.cpp graph.cpp maxflow.cpp -o BetheGams_mex.mexmaci64

makePotential_mex.mexmaci64: makePotential_mex.cpp
	mex -largeArrayDims -O makePotential_mex.cpp -o makePotential_mex.mexmaci64

clean:
	rm *.o *.mexmaci64

