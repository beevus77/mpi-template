#include <mpi.h>
#include <cmath>
#include <iostream>
using std::cout;
using std::min;
using std::endl;
using std::string;
using std::stoi;


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


    // Initialize A (horizontal strips)
    // TODO for speedup: flip matrix if m > n (tall)
    int m = stoi(argv[1]), n = stoi(argv[2]);
    int numrows = ceil(m/P), numcols=n;
    if (ID < (m % int(ceil(m/P))) ) {
        numrows++;
    }
    double myA[numrows][numcols];
    for (int i = 0; i < numrows; i++) {
        for (int j = 0; j < n; j++) {
            double val_i = ID*ceil(m/P) + i;

            // increment val_i based on ID
            if (ID < (m % int(ceil(m/P))) ) {
                val_i += ID;
            } else {
                val_i += m % int(ceil(m/P));
            }

            // set myA
            myA[i][j] = j*sin(val_i) + val_i*cos(j) + sqrt(val_i+j+1);
        }
    }


    // Use MPI_BARRIER then call MPI_WTIME on root process
    MPI_Barrier(MPI_COMM_WORLD);
    double start, end;
    if (ID == 0){
        start = MPI_Wtime();
    }

    // Iterate 10 times
    // TODO

    // Compute verification values, send to root process, FIXME: verification values not correct
    double vals[2] = {0,0}; // first entry: sum of all entries, second entry: sum of squares
    for (int i = 0; i < numrows; i++) {
        for (int j = 0; j < numcols; j++) {
            vals[0] += myA[i][j];
            vals[1] += pow(myA[i][j],2);
        }
    }
    double verifications[2];
    MPI_Reduce(&vals, &verifications, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Call MPI_WTIME on root process, print elapsed time and verification
    if (ID == 0) {
        end = MPI_Wtime();
        cout << end-start << " seconds" << endl;
        cout << "Sum of entries of A (verification 1):\t\t" << verifications[0] << endl;
        cout << "Sum of squares of entries of A (verification 2):\t" << verifications[1] << endl;
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}