#include<iostream>
#include<cmath>
#include<algorithm>

using namespace std;
  
double fdmEs_callOption(double s0, double k, double T, double vol, double r); //declaration
int main()
{
    cout << fdmEs_callOption(100,100,1,0.2,0.05) << endl;
    return 0;
}
 
double fdmEs_callOption(double s0, double k, double T, double vol,
    double r) { //definintion
    int i, j, iLower;
    double value;
    const int I = 100; //how many underlying's steps there are
    const double DELTA_S = 2 * k / I;
    /* A simple choice. But it has an implication: the cap on S is now effectively
     chosen to be: 2k.  */

    iLower = (int) floor(s0 / DELTA_S); // First available oppo for giving iLower a value

    double delta_t = 1 / (vol * vol * I * I); /* to do with the

    maths behind the explicit scheme ('stability')*/

    int J = (int) floor(T / delta_t);

    delta_t = T / J; //now delta_t is final !!
 
    double* S; S = new double[I + 1];

    double* Deltas; Deltas = new double[I + 1];

    double* Gammas; Gammas = new double[I + 1];

    double* Vprev; Vprev = new double[I + 1];

    double* Vcurr; Vcurr = new double[I + 1];
 
    // At maturity (i.e. for j==0):

    for (i = 0; i <= I; i++) {

        S[i] = i * DELTA_S; // Won't change !

        Vcurr[i] = max(S[i] - k, 0.0);
        /* Will keep getting

                updated in the looping below, except for

                the lower & upper boundaries

                (ie for i==0 & i==I respectively)

                which will NOT be changed !! */

    }
 
    //Now for j>0: the Big Loop

    for (j = 1; j <= J; j++) {

        for (i = 0; i <= I; i++)

            Vprev[i] = Vcurr[i];

        //Now deal with the unknows (of which there I-1)

        for (i = 1; i <= I - 1; i++) {
            
            // 0 <= i <= I, so we have I+1 values, for i = 0 and i = T, we have I-1 unknowns.

            //Implementing eqn 12 (2 preps, plus eq 12 itself)

            Deltas[i] = (Vprev[i + 1] - Vprev[i - 1]) / (2 * DELTA_S);

            Gammas[i] = (Vprev[i + 1] - 2 * Vprev[i] + Vprev[i - 1]) / (DELTA_S * DELTA_S);

            Vcurr[i] = Vprev[i] + delta_t *

                (0.5 * vol * vol * S[i] * S[i] * Gammas[i]

                + r * S[i] * Deltas[i] - r * Vprev[i]);// eqn 12 itself!

        } //end of the i loop

    }//end of the j loop
 
    // Last available oppo for giving iLower a value. But we've done it!

    value = Vcurr[iLower]; // Could be improved!

    delete[] S; delete[] Deltas; delete[] Gammas; delete[] Vprev;

    delete[] Vcurr;

    return value;

}
 
