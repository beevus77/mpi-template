#include <mpi.h>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;


// Function definition matching specs
double f(double x)
{
    double y = 2*x;
    for (int i = 1; i <= 10; i++) {
        y += x*cos(y+i) / pow(1.5,i);
    }
    return y;
}


int main(int argc, char **argv)
{

    // Initialize the MPI environment, get number of processors and processor ID
    MPI_Init(&argc, &argv);
    int P, ID;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);

    // TODO: Initialize A

    // Use MPI_BARRIER then call MPI_WTIME on root process
    MPI_Barrier(MPI_COMM_WORLD);
    double start, end;
    if (ID == 0){
        start = MPI_Wtime();
    }

    // Iterate 10 times, FIXME: not actual iteration
    double val1 = 0;
    for (int i = 0; i < 10; i++) {
        val1 += f(ID + i);
    }
    double val2 = pow(val1, 2);

    // Compute verification values, send to root process, FIXME: verification values not correct
    double verification1, verification2;
    MPI_Reduce(&val1, &verification1, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&val2, &verification2, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Call MPI_WTIME on root process, print elapsed time and verification
    if (ID == 0) {
        end = MPI_Wtime();
        cout << end-start << " seconds" << endl;
        cout << "Sum of entries of A (verification 1):\t\t" << verification1 << endl;
        cout << "Sum of squares of entries of A (verification 2):\t" << verification2 << endl;
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}