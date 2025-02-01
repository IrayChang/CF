#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

double maxi(double x, double y)
{
    return (x > y) ? x : y;
}

/* gBM
 Bullet formula: gBM: S_t = S_0 * exp((r- 0.5 * vol^2)*t + vol * W_t) where W_t ~ N(0,t)
 Iterative/recursive formula: S_t2 = S_t1 * exp((r- 0.5 * vol^2)*(t2-t1)) + vol * dWt)
 where dWt ~ N(0,t2-t1).
 Note: (non overlapping) Brownian increments are indepedents! */


double BarrierOptionPricer (double s0, double k, double T, double r, double vol, double B, double n)
{
    int i, j;
    double dt = (double) 1/365;
    double minSt; // mere initialisation
    int m = (int) floor(T/dt); //m = number of days till maturity
    cout << m << endl;
    double x = 1;
    double st;
    double forPayoff = 0;
    for (i = 1; i <= n; i++) // a simulaiton counter
    {
        st = s0; // resetting st back to the spot
        minSt = s0;
        for (j = 1; j <= m; j++) // a time step counter
        {
            st *= exp((r - 0.5 * vol * vol) * dt + vol * x);
            if (st < minSt)
                minSt = st;
        }
        if (minSt > B)
            forPayoff += maxi(k-st,0); //put's payoff
    }
    return exp(-r*T) * forPayoff / n ; //risk neutral growth factor
}

int main()
{
    cout << BarrierOptionPricer(100, 110*exp(0.05*1), 1, 0.05, 0, 90, 100) << endl;
}
