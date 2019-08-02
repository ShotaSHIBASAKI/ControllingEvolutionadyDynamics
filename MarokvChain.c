//
//  main.c
//  Dynamic programming with Markov chian
//  long term optimal introduction rates of the sensitive and resistant cooperators, m1 and m2, respectively
//  Created by Shota on 19.11.18.
//  Copyright © 2018 ShotaS. All rights reserved.
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define num_state 14//in model 0 number of  state is 14;model1 7
#define stepsize 100//stepsize of m1 and m2


void TransitionMatrix(int model, double m1, double m2, double mu1, double mu2, double P[num_state][num_state]);
void InnerProd(double pi[num_state], double P[num_state][num_state]);
void StationaryDistribution(int model, double m1, double m2, double mu1,double mu2, double stationary[num_state]);
double Productivity(double prod, double phi[num_state], double pi[num_state]);
char fname[256];


/*main function*/
int main(void) {
    int t;// index for time step
    int i;// index for state
    int j;// index for m1
    int k;// index for m2
    int T_end = 2500;//end of the simulation with discrete time step
    double norm=stepsize;//normalization for m1 and m2
    int counter = 0;
    // define the parameter values
    int model = 0; /*model type. 0: Ergodic. 1: Non-ergodic. See diagram for more detail.*/
    double mu1=0.01;//mu1=0.1 for test, 0.01 for true running
    double mu2=0.00;// Assum[ption no mutation in resistance
    double phi[num_state];//productivity at each state
    
    //define vector of optimal m1 and m2
    double opt_m1_list[T_end];
    double opt_m2_list[T_end];
    
    double prob_list[num_state][(stepsize + 1) * (stepsize + 1)];//define the prob. vector list; state x stepsize**2
    
    double product_list[stepsize+1][stepsize+1];//productivity given value of m1 and m2
    
    double opt_m1;//optimal value of m1 at time step t
    double opt_m2;//optimal value of m2 at time step t
    double opt_prod;//optimal productivity at time step t
    
    for(j=0;j<=stepsize;j++)
    {
        for(k=0;k<=stepsize;k++)
            product_list[j][k]=0.0;
    }
    double prod_exp;//expected productvity
    if (model==0)
    {
        mu2 = 0.0;
        /*define the productivity at each state
         This values can be calculated by another simulation,; however we assume these values are given here*/
        //phi>0 as phi:= alpha(T-Tin)
        //Ergodic MC as in Fig. 5A in the main text
        phi[0] = 0.4;//staying sCo
        phi[1] = 0.2;//mutation sCo -> sCh
        phi[2] = 0.0;//staying sCh
        phi[3] = 0.15;//introduction sCh -> rCo
        phi[4] = 0.3;//staying rCo
        phi[5] = 0.15;//mutation  rCo -> rCh
        phi[6] = 0.0;//staying rCh
        phi[7] = 0.2;//introduction rCh -> sCo
        phi[8] = 0.35;//introduction sCO -> rCo
        phi[9] = 0.35;//introduction rCo->sCh
        phi[10] = 0.3;//mistake-introduction sCo -> sCo (unnecessary)
        phi[11] = 0.0;//misake-introduction sCh -> sCo (failed)
        phi[12] = 0.2;//mistake-introduction rCo -> rCo (unnecessary)
        phi[13] = 0.0;//mistake-introduction rCh ->rCo (failed)
    }
    else
    {
        /*Non-ergodic MC as in FIg. A.11A in the appendix
         To simplify, ignoring transient states*/
        phi[0] = 0.4;//staying sCo (optimal)
        phi[1] = 0.2;//sCo (not optimal)
        phi[2] = 0.0;//staying sCh
        phi[3] = 0.3;//staying rCO (optimal)
        phi[4] = 0.15;//staying rCo (not optimal)
        phi[5] = 0.0;//staying rCh
        phi[6] = 0.2;//Coexistence of sCO with rCh
        /*Ignore below*/
        phi[7] = 0.0;
        phi[8] = 0.0;
        phi[9] = 0.0;
        phi[10] = 0.0;
        phi[11] = 0.0;
        phi[12] = 0.0;
        phi[13] = 0.0;
       
    }
    
    
    
    
    //define the initial state distribution. pi_i=1 if i=1, otherwise pi_i=0
    double pi[num_state];
    //printf("Set initial conditions for the state prob. \n");
    for(i=0;i<num_state;i++)
    {
        if(i==0)
        {
            pi[i] = 1.0;//initial condition is always sCo
            //printf("pi_%d(0) is %.1f\n",i, pi[i]);
        }
        else
        {
            pi[i] = 0.0;
            //printf("pi_%d(0) is %.1f\n",i, pi[i]);
        }
    }
    
    printf("sum of mu %.2f \n", mu1+mu2);
    /*Analysis 1: check whether prob. distribution converge to the analytcal solution or not*/
    double m1 = 0.3;
    double m2 = 0.2;//values for calculating stationary distribution
    // defining stationary distribution and transition matrix depending on the model to check the converegenec
    double P[num_state][num_state];// transition matrix
    TransitionMatrix(model,  m1, m2, mu1,  mu2, P);//set the transition matrix
    
    //probability vector for a long time step later
    for(t=0;t<T_end;t++)
    {
        InnerProd(pi, P);
    }
    double stationary[num_state];//stationary distribution
    StationaryDistribution(model, m1, m2, mu1, mu2, stationary);//set stationary distribution
    //comparing the analytical results with numerical ones
    double diff=0.0;
    for(i=0;i<num_state;i++)
    {
        if (pi[i]-stationary[i]>0)
        {
            diff+=(pi[i]-stationary[i]);
        }
        else
        {
            diff+=(stationary[i]-pi[i]);
        }
    }
    if(diff < num_state * 0.00001)
    {
        printf("Convereged! Go ahead! \n");
        for (i=0; i<num_state;i++)
        {
            //printf("%d th pi is %.3f and stationary is %.3f. \n", i, pi[i], stationary[i]);
        }
    }
    else
    {
        printf("The stationary distribution is wrong or too short time for converegence. \n");
        for (i=0; i<num_state;i++)
        {
            //printf("%d th pi is %.3f and stationary is %.3f. \n", i, pi[i], stationary[i]);
        }
        //return 1;
    }
    
    
    /*Analysis 2: calculating optimal introduction rates*/
    //define probaility vector list
    for(j=0; j < (stepsize+1) * (stepsize+1); j++)
    {
        for(i=0;i<num_state; i++)
        {
            if(i==0)
            {
                prob_list[i][j]=1.0;
            }
            else
            {
                prob_list[i][j]=0.0;
            }
        }
    }
    for(t=1;t<=T_end;t++)
    {
        opt_m1 = 0.0;//optimal value of m1 at time step t
        opt_m2 = 0.0;//optimal value of m2 at time step t
        opt_prod = -100;//optimal productivity of m1 at time step t
        counter=0;
        for(j=0;j<=stepsize;j++)
        {
            m1 = j / norm;
            //printf("m1 = %.2f and",m1);
            for(k=0;k<=stepsize;k++)
            {
                m2 = k / norm;
                //printf(" m2 = %.2f \n",m2);
                if(m1+m2+mu1+mu2<1)
                {
                    //this is necessary for ergodic s Markov chain
                    TransitionMatrix(model,  m1, m2, mu1,  mu2, P);//set markov chain
                    //printf(" m2 = %.2f. \n",m2);
                    //get the probability vector at t-1
                    for(i=0;i<num_state;i++)
                    {
                        pi[i]=prob_list[i][counter];
                        //printf("pi_%d = %.3f ", i, pi[i]);
                    }
                    //printf("\n");
                    //calculate probability vector at t
                    InnerProd(pi, P);
                    //calculate the expected productivity
                    prod_exp=Productivity(product_list[j][k], phi, pi);
                    //comparing current optimum
                    if(prod_exp>opt_prod)
                    {
                        //printf("Current condition m1 = %.1f and m2 = %.1f is  better because %.3f > %.3f \n", m1, m2, prod_exp, opt_prod);
                        opt_m1=m1;
                        opt_m2=m2;
                        opt_prod=prod_exp;
                    }
                    else
                    {
                        //printf("Current condition m1 = %.1f and m2 = %.1f is not good because %.3f < %.3f \n", m1, m2, prod_exp, opt_prod);
                    }
                    //update the probability state and productivity list
                    for(i=0;i<num_state;i++)
                    {
                        prob_list[i][counter]=pi[i];
                    }
                    product_list[j][k]=prod_exp;
                }
                counter+=1;
            }
        }
       //save the optimal values at time step t. Note t begins from 1
        opt_m1_list[t-1]=opt_m1;
        opt_m2_list[t-1]=opt_m2;
        printf("at t=%d, opt m1 = %.2f and m2 = %.2f and max prod is %.4f. \n", t, opt_m1, opt_m2, opt_prod);
        
    }
    /*Show the distribution at the end*/
    /*TransitionMatrix(model,  opt_m1, opt_m2, mu1,  mu2, P);
    for(i=0;i<num_state;i++)
    {
        if(i==0)
        {
            pi[i] = 1.0;//initial condition is always sCo
            //printf("pi_%d(0) is %.1f\n",i, pi[i]);
        }
        else
        {
            pi[i] = 0.0;
            //printf("pi_%d(0) is %.1f\n",i, pi[i]);
        }
    }
    for(t=0.0;t<=T_end;t++)
    {
        InnerProd(pi, P);
    }
    for(i=0;i<num_state;i++)
    {
        printf("state %d prob is %,3f \n", i, pi[i]);
    }*/
    /*Analysis 3: Calculate the optimal rates from the prediction*/
    opt_m1 = 0.0;
    opt_m2 = 0.0;
    opt_prod = -100;
    /*Notice m1 and m2 cannot be 0 if they give the optimal productivity*/
    for(j=1;j<=stepsize;j++)
    {
        m1 = j / norm;
        for(k =1; k<=stepsize;k++)
        {
            m2 = k / norm;
            if(m1+m2<1-mu1-mu2)
            {
                StationaryDistribution(model, m1, m2, mu1, mu2, stationary);
                prod_exp=Productivity(0.0, phi, stationary);
                if(prod_exp>opt_prod)
                {
                    //printf("Current condition m1 = %.1f and m2 = %.1f is  better because %.3f > %.3f \n", m1, m2, prod_exp, opt_prod);
                    opt_m1=m1;
                    opt_m2=m2;
                    opt_prod=prod_exp;
                }
                else
                {
                    //printf("Current condition m1 = %.1f and m2 = %.1f is not good because %.3f < %.3f \n", m1, m2, prod_exp, opt_prod);
                }
            }
        }
    }
    printf("Long term opt m1 = %.2f and m2 = %.2f \n", opt_m1, opt_m2);
    
    /*save csv file*/
    printf("Start making a csv file. \n");
    FILE *fp;
    if(model==0){
        sprintf(fname,"optDPMC_ergodic.csv");
    }
    else{
        sprintf(fname,"optDPMC_Nonergodic.csv");
    }
    fp=fopen(fname, "w");
    fprintf(fp, "time, m1, m2\n");
    for(t=1;t<=T_end;t++)
    {
        fprintf(fp, "%d, %.3f, %.3f\n",t,opt_m1_list[t-1], opt_m2_list[t-1]);
    }
    fprintf(fp, "infinity, %.3f, %.3f\n",opt_m1, opt_m2);
    fclose(fp);
    printf("finish writing a csv file. \n");
   
    
    return 0;
}






/*Analyze the Markov Chain*/
void TransitionMatrix(int model, double m1, double m2, double mu1, double mu2, double P[num_state][num_state])
{
    //Get the transition matrix
    int i, j;
    //First, set P as a zero-matrix because it is a sparse matrix
    for(i=0;i<num_state;i++)
    {
        for(j=0;j<num_state;j++)
        {
            P[i][j]=0.0;
        }
    }
    if(model==0)
    {
        /*Ergodic*/
        P[0][0] = 1.0 - m1 - m2 - mu1;
        P[7][0] = 1.0;
        P[9][0] = 1.0;
        P[10][0] = 1.0;
        P[0][1] = mu1;
        P[1][2] = 1.0;
        P[2][2] = 1.0 - m1 - m2;
        P[11][2] = 1.0;
        P[2][3] = m2;
        P[3][4] = 1.0;
        P[4][4] = 1.0-m1-m2-mu1;
        P[8][4] = 1.0;
        P[12][4] = 1.0;
        P[4][5] = mu1;
        P[5][6] = 1.0;
        P[6][6] = 1.0 - m1 - m2;
        P[13][6] = 1.0;
        P[6][7] = m1;
        P[0][8] = m2;
        P[4][9] = m1;
        P[0][10] = m1;
        P[2][11] = m1;
        P[4][12] = m2;
        P[6][13] = m2;
    }
    else
    {
        /*Non-ergodic*/
        P[0][0] = 1.0 - m2 - mu1;
        P[1][0] = m1;
        P[0][1] = m2;
        P[1][1] = 1.0 - m1 - mu1;
        P[0][2]= mu1;
        P[1][2]= mu1;
        P[2][2]= 1.0-m2;
        P[6][2]= mu1;
        P[2][3]= m2;
        P[3][3]= 1.0 - m1 - mu1;
        P[4][3]= m2;
        P[6][3]= m2;
        P[3][4]= m1;
        P[4][4]= 1.0 - m2 - mu1;
        P[3][5]= mu1;
        P[4][5]= mu1;
        P[5][5]= 1.0 - m1;
        P[5][6]= m1;
        P[6][6]= 1 -m2 -mu1;
    }
    
}
void InnerProd(double pi[num_state], double P[num_state][num_state])
{
    //innner product for markov chain
    int ir, ic;//index for row and colmun
    double next[num_state];
    for(ic=0;ic<num_state;ic++)
    {
        next[ic] = 0.0;
        for(ir=0;ir<num_state;ir++)
        {
            next[ic]+= pi[ir] * P[ir][ic];
        }
    }
    for (ic=0;ic<num_state;ic++)
    {
        pi[ic] = next[ic];
    }
}
void StationaryDistribution(int model, double m1, double m2, double mu1,double mu2, double stationary[num_state])
{
    //calculating stationry distribution
    double coeff[num_state];//coefficient vector
    double sum_coeff=0.0;
    int i;
    double baseline;//baseline
    if(model==0)
    {
        /*Ergodic*/
        coeff[0] = 1.0;//pi1
        coeff[1] = mu1 * coeff[0];//pi2
        coeff[2] = mu1 / m2;//pi3
        coeff[3] = m2 * coeff[2];//pi4
        coeff[4] = (mu1 + m2) / (mu1 + m1);//pi5
        coeff[5] = mu1 * coeff[4];//pi6
        coeff[6] = mu1 * (mu1 + m2) / m1 /  (mu1 + m1);//pi7
        coeff[7] = m1 * coeff[6];//pi8
        coeff[8] = m2 * coeff[0];//pi9
        coeff[9] = m1 * coeff[4];//pi10
        coeff[10] = m1 * coeff[0];//pi11
        coeff[11] = m1 * coeff[2]; //pi12
        coeff[12] = m2 * coeff[4];//pi13
        coeff[13] = m2 * coeff[6];//pi14
        
        
    }
    else
    {
        /*Non ergodic*/
        coeff[0] = 0.0;//pi1
        coeff[1] = 0.0;//pi2
        coeff[2] = mu1 / m2;//pi3
        coeff[3] = (m2 + mu1) * (m2 + mu1) / (mu1 * (m1 + m2 + mu1));//pi4
        coeff[4] = m1 /(m2 + mu1) * coeff[3];//pi5
        coeff[5] = (m2 + mu1) / m1;//pi6
        coeff[6] = 1.0;//pi7
        coeff[7] = 0.0;
        coeff[8] = 0.0;
        coeff[9] = 0.0;
        coeff[10] = 0.0;
        coeff[11] = 0.0;
        coeff[12] = 0.0;
        coeff[13] = 0.0;
    }
    for(i=0;i<num_state;i++)
    {
        sum_coeff+=coeff[i];
    }
    baseline=1/sum_coeff;
    for(i=0;i<num_state;i++)
    {
        stationary[i] = coeff[i] * baseline;
    }
}
/*Expected productivity*/
double Productivity(double prod, double phi[num_state], double pi[num_state])
{
    int i;
    double next;
    next=prod;
    for(i=0;i<num_state;i++)
    {
        next+=phi[i]*pi[i];
    }
    return next;
}
