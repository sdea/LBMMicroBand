#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
typedef double T;

// Define descriptor for lattice
#define DES descriptors::AdvectionDiffusionD2Q5Descriptor

// Store here paramters
struct Parameters {

     // Lattice dimensions
     plint nx;
     plint ny;

     //Relaxation time
     T omega;

     // Real time
     T realTime, dt;
};

void assigParam (Parameters &param) {

     param.nx = 598;
     param.ny = 400;

     param.omega = 1.;

     param.realTime = 100; // [s]
     param.dt = 1E-04; // [s]

}

// Setup simulation
void setupSim (MultiBlockLattice2D<T, DES> &adLattice,
               OnLatticeAdvectionDiffusionBoundaryCondition2D<T, DES> *adBC,
               Parameters &param) {
    
    // Lattice dimeunsions
    plint nx = adLattice.getNx();
    plint ny = adLattice.getNy();
    
    // Parameters for AD
    Array<T, 2> adU0(0.,0.);
    T adRho = 0.;                      // Init at 0.5 to fast convergence

    pcout << "SETUP SIMULATION...\n" << std::endl;

    // Set here boundary conditions
    // On all bottom for now
    Box2D bottom(0,nx-1,  0,0);

    adBC->addTemperatureBoundary1N(bottom, adLattice); 
    setBoundaryDensity(adLattice, bottom, 1.);
    initializeAtEquilibrium(adLattice, adLattice.getBoundingBox(), adRho, adU0);

    // Init lattice: Execute the data processor before streaming
    adLattice.initialize();
    delete adBC;
}

void writeVTK(MultiBlockLattice2D<T, DES>& lattice,
              plint iter)
{    
    T dx = 1.;
    T dt = 1.;
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeDensity(lattice), "Concentration", dx/dt);
}

int main(int argc, char* argv []) {

	// Init palabos simulation
    	plbInit(&argc, &argv);
	pcout << "Hello palabos" << std::endl;

	// Call struct and parameters
	Parameters param;
	assigParam(param);

	// Dimensions of grid (LU)
	const plint nx = param.nx;
	const plint ny = param.ny;
	const T omega = param.omega;

	// Instantiate lattice
	MultiBlockLattice2D<T, DES> lattice(nx, ny, new AdvectionDiffusionBGKdynamics<T, DES>(omega));
	
	// Setup simulation with BC
	setupSim(lattice, createLocalAdvectionDiffusionBoundaryCondition2D<T,DES>(), param);	

	// Print here 2D vtk 
        writeVTK(lattice, 0);	
	plint numIter =  (int)(param.realTime/param.dt);
	pcout << "Simulating  " << numIter << " iterations" <<std::endl;
	for(plint iTT=0; iTT < numIter; ++iTT) {

	    if(iTT%10 == 0) {
	       
		    pcout << "iTT:  " << iTT << std::endl;
            }
    	    
            // Print output	    
	    if(iTT%1000 == 0) {
	       
		    pcout << "Writng VTK..." << std::endl; 
		    writeVTK(lattice, iTT);
            }	
	    
	    
	    lattice.collideAndStream();
	
	}

}

