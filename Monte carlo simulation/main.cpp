#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

/*
 
 The defualt time, tau, follows an exponential distribution
 
 To simulate the exponential distribution, we introduce the monte carlo to simulate the
 credit risk default time, to try to estimate the probability of an event
 
 STEPs:
 Generate U ~ U[0,1]
 Compute tau
 Average over many simulations estimated default time
 
*/

double expectedDefaultTime(double lambda, int n){
    int i;
    double u, tauTimeslambda, summation = 0;
    for(i=1; i<=n; i++){
        u = (double) rand()/RAND_MAX; // why use uniform distribution?
        if(u == 1)
            u = 0.99999;
        tauTimeslambda = -log(1-u); // tau = -ln(1-U) / lambda
        
        /* but for the computional efficiency, we times lambda in loop first, and divide lambda outisde the loop.  */
        
        summation += tauTimeslambda;
    }
    return summation / lambda / n;
}

double expectedDefaultTime_refined(double lambda, int n){
    int i;
    double u, tauTimeslambda, summation = 0;
    for(i=1; i<=n; i++){
        u = (double) rand()/RAND_MAX;
        if(u == 0)
            u = 0.0000001;
        tauTimeslambda = -log(u); // tau = -ln(1-U) / lambda
        summation += tauTimeslambda;
    }
    return summation / lambda / n;
}

void SensitivityAnalysis()
{
    int i, n = 1000;
    double McVal = 0, TheoVal = 0, prevMcVal, prevTheoVal;
    double lambda = 0.00001;
    for (i = 0 ; i < 5; i++)
    {
        prevMcVal = McVal;
        prevTheoVal = TheoVal;
        srand(5); //to elimate cross-sample noise
        McVal = expectedDefaultTime_refined(lambda, n);
        TheoVal = 1/lambda; //Expectation of hazard rate
        lambda *= 10;
        if (i != 0)
        {
            cout << "The empirical ratio is " << McVal / prevMcVal << endl;
            cout << "The theo ratio is " << TheoVal / prevTheoVal << endl;
        }
    }
}

int main(){
    SensitivityAnalysis();
    srand(5);
    double lamda = 0.01/(1-0.04);
    cout << "Expected default time is " << expectedDefaultTime(lamda, 1000) << endl;
    unsigned seed;
    seed = (unsigned)time(NULL);
    srand(seed);
    cout << "Expected default time is " << expectedDefaultTime_refined(lamda, 1000) << endl;
    return 0;
}
