#include <iostream>
#include <cmath>
using std::cout;
using std::min;
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


int main(int argc, char **argv)
{
    // cout << f(1.) << endl;
    // cout << min(3.4, ceil(3.3)) << endl;
    // cout << 20 % int(ceil(20./3)) << endl;
    // cout << argv[1] << endl;

    int m = 20, n = 20, P=4, ID = stoi(argv[1]);
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
            // cout << "myA[" << i << "][" << j << "]\t" << myA[i][j] << endl;
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