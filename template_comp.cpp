#include <mpi.h>
#include <iostream>
using std::cout;
using std::endl;


int main(int argc, char **argv)
{

    // Initialize the MPI environment.
    MPI_Init(&argc, &argv);

    // Initialize A

    // Use MPI_BARRIER then call MPI_WTIME on root process

    // Iterate 10 times

    // Compute verification values, send to root process

    // Call MPI_WTIME on root process

    // Print out elapsed time and verification

    // Finalize MPI
    MPI_Finalize();
    return 0;
}


double f(double x)
{
    
}