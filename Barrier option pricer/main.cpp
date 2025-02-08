#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

double maxi(double x, double y)
{
    return (x > y) ? x : y;
}

double BarrierOptionPricer (double s0, double k, double r, double vol, double T, double B, double n);
double Phi(double x);
double PhiInverse(double p);
double EuroPut(double s0, double k, double r, double vol, double T);
double EuroCall(double s0, double k, double r, double vol, double T);

void SensitivityAnalysis_barrier(double s0, double k, double r, double vol, double T, int n);

int main()
{
    cout << "Checking put-call parity: "
    << EuroCall(100, 100, 0.05, 0.2, 1) - EuroPut(100, 100, 0.05, 0.2, 1) << " v.s. "
    << 100 - 100 * exp(-0.05 * 1) << endl;
    
    srand(5);
    SensitivityAnalysis_barrier(100, 100, 0.05, 0.2, 1, 1000);

    cout << "The barrier option price is " << BarrierOptionPricer(100, 110 * exp(0.05*1), 0.05, 0, 1, 90, 100) << endl;
    cout << EuroPut(100, 100, 0.05, 0.2, 1) << endl; // theoretical value when barrier is 0
    cout << EuroCall(100, 100, 0.05, 0.2, 1) << endl;
    
    if (PhiInverse(0.5) == 0)
        cout << "PhiInverse works correctly" << endl;
    /*
    cout << Phi(0) << endl; // surely 0.5
    cout << PhiInverse(0.5) << endl; // surely 0
    */
}


/* gBM
 Bullet formula: gBM: S_t = S_0 * exp((r- 0.5 * vol^2)*t + vol * W_t) where W_t ~ N(0,t)
 Only focuses on now and maturity, no time dependent
  
 Iterative/recursive formula: S_t2 = S_t1 * exp((r- 0.5 * vol^2)*(t2-t1)) + vol * dWt)
 where dWt ~ N(0,t2-t1).
 Observe the process among time path, spliting time into time steps, delta
 
 Note: (non overlapping) Brownian increments are indepedents!
 Note: Wiener process is a continous formula of random walk, where t_1 - t_0 = 0,
       it is a stochastic process
 */


/* barrier option
 Payoff = max(k - St, 0), if I(minSt > B), where I is the indicator function
 
 The stock price, St, follows a Geometric Brownian Motion, gBM
 
 */


/*
 Phi and PhiInverse:
 In the gBM assumption, dW_t ~ N(0, dt),
 to simulate standard normal variables, u ~ U(0, 1)
 and by inversing the cummulative normal distribution, to obtain standard normal samples for simulating the Wiener process.
 
 
 Note: standardization, z = (x - u) / s.d, where z ~ N(0, 1)
       Reverse process, x = z * s.d,
       Then, we obtain a normally distributed variable
 
       If u ~ uniform distribution, we can trasform it into standard normal distribution
 */


double BarrierOptionPricer (double s0, double k, double r, double vol, double T, double B, double n)
{
    int i, j;
    double dt = (double) 1/365;
    double minSt; // mere initialisation
    int m = (int) floor(T/dt); //m = number of days till maturity
    // cout << m << endl;
    double x = 1;
    double st;
    double forPayoff = 0;
    double u,z;
    
    for (i = 1; i <= n; i++) // a simulaiton counter, i is a index but n is a number
    {
        st = s0; // resetting st back to the spot
        minSt = s0; //recursive formula, to check condition in each time step, such as daily or monthly basis
        for (j = 1; j <= m; j++) // a time step counter
        {
            do
            {
                u = (double) rand() / RAND_MAX;
            } while (u==0| u==1);
            
            z = PhiInverse(u);
            x = z * sqrt(dt);
            
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
    
    // use binary search method to approximate the probability
    // a six-sigma event
    
    a = -10; //left boundary
    b = 10; //right boundary
    i = 0; //iterations
    
    do
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
    {
        cout << "CAP has been reached" << endl;
    }
    
    // cout << "Btw, it took " << i << " iterations." << endl;

    return m;
}

double EuroPut(double s0, double k, double r, double vol, double T)
{
    double d1, d2;
    double VolSqrt_T = (vol * sqrt(T));
    d1 = (log(s0 / k) + (r + 0.5 * vol * vol) * T) / VolSqrt_T;
    d2 = d1 - VolSqrt_T;
    
    return Phi(-d2) * k * exp(-r * T) - Phi(-d1) * s0;
}

double EuroCall(double s0, double k, double r, double vol, double T)
{
    double d1, d2;
    double VolSqrt_T = (vol * sqrt(T));
    d1 = (log(s0 / k) + (r + 0.5 * vol * vol) * T) / VolSqrt_T;
    d2 = d1 - VolSqrt_T;
    
    return Phi(d1) * s0 - Phi(d2) * k * exp(-r * T);
}

void SensitivityAnalysis_barrier (double s0, double k, double r, double vol, double T, int n)
{
    int i;
    int B = 0;
    int NumberOfScenerio = 11;
    int s0_over_NumberOfScenerio_Minus1 = s0 / (NumberOfScenerio - 1);
    
    for (i = 0; i < NumberOfScenerio; i++)
    {
        srand(5);
        cout << "Scenerio " << i << ": note's price = " << BarrierOptionPricer(s0, k, r, vol, T, B, n) << endl;
        B += s0_over_NumberOfScenerio_Minus1;
    }
}
