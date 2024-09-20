#include <iostream>
#include <cmath>
using std::cout;
using std::endl;
using std::min;
using std::vector;


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

    int m = 20, n = 20, P=3, ID = 2;
    vector<vector<double>> myA;
    if (ID == P-1) {
        double myA[m % int(ceil(m/P))][n];
    } else {
        double myA[int(ceil(m/P))][n];
    }
    for (int i = ID*int(ceil(m/P)); i < min(m, (ID+1)*int(ceil(m/P))); i++) {
        for (int j = 0; j < n; j++) {
            myA[i][j] = j*sin(i) + i*cos(j) + sqrt(i+j+1);
        }
    }
    cout << myA << endl;
}