#include <iostream>
#include <cmath>
using std::cout;
using std::min;
using std::max;
using std::endl;
using std::string;
using std::stoi;


double f(double x)
{
    double y = 2*x;
    for (int i = 1; i <= 10; i++) {
        y += x*cos(y+i) / pow(1.5,i);
    }
    return y;
}


void testSum(int ID)
{
    int m = 20, n = 20, P=4;
    int numrows = ceil(m/P), numcols=n;
    if (ID < (m % int(ceil(m/P))) ) {
        numrows++;
    }
    double myA[numrows][numcols];
    for (int i = 0; i < numrows; i++) {
        for (int j = 0; j < n; j++) {
            double val_i = ID*ceil(m/P) + i;

            if (ID < (m % int(ceil(m/P))) ) {
                val_i += ID;
            } else {
                val_i += m % int(ceil(m/P));
            }

            myA[i][j] = j*sin(val_i) + val_i*cos(j) + sqrt(val_i+j+1);
        }
    }

    double vals[2] = {0,0};
    for (int i = 0; i < numrows; i++) {
        for (int j = 0; j < numcols; j++) {
            vals[0] += myA[i][j];
            vals[1] += pow(myA[i][j],2);
        }
    }
    cout << vals[0] << endl;
    cout << vals[1] << endl;
}


void testBorderLogic(int ID)
{
    // Initialize myA
    int numrows = 5, numcols = 5;
    double myA[numrows][numcols];
    for (int i = 0; i < numrows; i++) {
        for (int j = 0; j < numcols; j++) {
            myA[i][j] = i*j;
        }
    }

    // initialize myNewA
    double myNewA[numrows][numcols];
    for (int i = 0; i < numrows; i++) {
        for (int j = 0; j < numcols; j++) {

            // Unchanged along border
            if (j == 0 or (ID == 0 and i == 0) or j == numcols-1) {
                myNewA[i][j] = myA[i][j];
            } else {
                myNewA[i][j] = -1;
            }
        }
    }

    for (int i = 0; i < numrows; i++) {
        for (int j = 0; j < numcols; j++) {
            cout << myNewA[i][j] << endl;
        }
    }
}


void testMaxMin()
{
    cout << max(3,4) << endl;
    cout << max(-1,-2) << endl;
    cout << max(-25, min(30, 5)) << endl;
}


int main(int argc, char **argv)
{
    // cout << f(1.) << endl;
    // cout << min(3.4, ceil(3.3)) << endl;
    // cout << 20 % int(ceil(20./3)) << endl;
    // cout << argv[1] << endl;

    testSum(stoi(argv[1]));

    testBorderLogic(stoi(argv[1]));

    testMaxMin();
    
}