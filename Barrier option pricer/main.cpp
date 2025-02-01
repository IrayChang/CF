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

double Phi(double x) //Numeric recipe for the N(0,1) CDF
{
    int neg = (x < 0);
    if (neg) x *= -1;
    double k(1 / (1 + 0.2316419 * x));
    double y = ((((1.330274429 * k - 1.821255978) * k + 1.781477937) * k - 0.356563782) * k + 0.319381530) * k;
    y = 1.0 - 0.398942280401 * exp(-0.5 * x * x) * y;
    return (1 - neg) * y + neg * (1 - y);
}

double PhiInverse(double p)
{
    double a, b, m, fOfm, diff;
    int i; const int CAP = 100; //'cap' on the num of iterations
    const double TOL = 0.0001; //basis point, the allowable range of error
    a = -10; //left boundary
    b = 10; //right boundary
    i = 0; //iterations
    
    do // use binary search method to approximate the probability
    {
        m = 0.5 * (a + b);
        fOfm = Phi(m); //CDF value
        diff = fOfm - p;
        
        if (diff<0) //(fOfm < p)
            a = m;
        else
            b = m;
        diff = fabs(diff); // absolute value of diff
        i++;//You can write it anywhere within this block/loop!
 
    } while (diff > TOL && i <= CAP);
    
    if (i >= CAP)
        cout << "CAP has been reached" << endl;
    cout << "Btw, it took " << i << " iterations." << endl;
    return m;
}


int main()
{
    cout << BarrierOptionPricer(100, 110*exp(0.05*1), 1, 0.05, 0, 90, 100) << endl;
}
