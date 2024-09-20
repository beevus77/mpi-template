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

    // int m = 20, n = 20, P=3, ID = stoi(argv[1]);
    // int numrows = ceil(m/P), numcols=n;
    // if (ID < (m % int(ceil(m/P))) ) {
    //     numrows++;
    // }
    // double myA[numrows][numcols];
    // for (int i = 0; i < numrows; i++) {
    //     for (int j = 0; j < n; j++) {
    //         double val_i = ID*ceil(m/P) + i;

    //         if (ID < (m % int(ceil(m/P))) ) {
    //             val_i += ID;
    //         } else {
    //             val_i += m % int(ceil(m/P));
    //         }

    //         myA[i][j] = j*sin(val_i) + val_i*cos(j) + sqrt(val_i+j+1);
    //         cout << "myA[" << i << "][" << j << "]\t" << myA[i][j] << endl;
    //     }
    // }
    // cout << 3*sin(15) + 15*cos(3) + sqrt(15+3+1) << endl;
    // cout << numrows << endl;
    // cout << ceil(m/P) << endl;
    // cout << m % int(ceil(m/P)) << endl;

    double vals[2] = {0,1};
    cout << vals[0] << endl;
}