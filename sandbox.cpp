#include <iostream>
#include <cmath>
using std::cout;
using std::endl;


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
    cout << f(1.) << endl;
}