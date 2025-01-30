#include <iostream>
#include <cmath>
using namespace std;

double btPricer(double s0, double k, double r, double vol, double T, int n) //definition
{
    double dt = T/n,
    nu = r - 0.5 * pow(vol,2),
    Detla = sqrt(pow(vol,2)*dt + pow(nu,2)*pow(dt,2)),
    p = 0.5 + 0.5 * nu * dt/Detla;

    double * St = new double [n+1]; //DMA, 動態分配一個 double 變數
    double * V = new double [n+1];

    int i, j; // Declaration of two characteristic: i & j

    // At maturity (i.e at j == n):
    // Xt[0] = x0 (i.e the Ln(the spot s0)) - n * Delta
    St[0] = s0 * exp(-n * Detla); // s0 gone downhill n times
    for(i=1; i<=n; i++){
        St[i] = St[i-1] * exp(2 * Detla);
    }
    for(i=0; i<=n; i++){
        V[i] = max(St[i] - k, 0.0); // call's payoff
    }

    // what happened after maturity
    // Travelling backward-in-time
    for(j=n-1; j>0; j--){
        for(i=0; i<=j; i++){
            V[i] = exp(-r*dt) * (p * V[i+1] + (1-p) * V[i+1] );
        }
    }
    double value = V[0];
    delete[] St; //freeing up the memeory, the opposite of new, 釋放記憶體
    delete[] V;
    return value; // dummy returm
}

// this is a comment

/* GBM
Binomial trees : determine the paremeter first, Moment Matching */

// call option
int main()
{
    cout << "Example europian call option's price: "
    << btPricer(200, 100, 0.05, 0.2, 1, 100) << endl;
    return 0;
}
