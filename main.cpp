#include <iostream>
#include "theorver.h"
#include<iomanip>
#include <conio.h>
#include <string> 
#include <cmath>      
using namespace std;
void start(){
    int n;
    cout<<endl;
    cout<<"Print a number to choose:"<<endl;
    cout<<"1 to find a combination "<<endl;
    cout<<"2 to find a discombination "<<endl;
    cout<<"3 to use Bernulli's formula "<<endl;
    cout<<"4 to use De Moivre Laplace integral theorem "<<endl;
    cout<<"5 to use De Moivre Laplace local theorem "<<endl;
    cout<<"6 to use "<<"qp"<<" function "<<endl;
    cout<<"7 to find a expected value"<<endl;
    cout<<"8 to find a variance"<<endl;
    cout<<"9 to find a standard deviation"<<endl;
    cout<<"0 to end "<<endl;
    cin>>n;   
    if(n==1)
        {
            int a,b;
            cout<<"Add a number of all objects"<<endl;
            cin>>a;
            cout<<"Add a number of the set"<<endl;
            cin>>b;
            cout<<"The combination is "<<comb(a,b)<<endl;
            cout<<"***"<<endl;
            start();
        }
    else if(n==2)
        {
            int a;
            cout<<"Add a number of the set"<<endl;
            cin>>a;
            cout<<"The discombination is "<<discomb(a)<<endl;
            cout<<"***"<<endl;
            start();
        }
    else if(n==3)
        {
            int a,b;
            float p;
            cout<<"Add a number of all objects"<<endl;
            cin>>a;
            cout<<"Add a number of the set"<<endl;
            cin>>b;
            cout<<"Add a probability"<<endl;
            cin>>p;
            cout<<"The combination is "<<bernulli(a,b, p)<<endl;
            cout<<"***"<<endl;
            start();
        }
    else if(n==4)
        {
            double a,b;
            float p;
            cout<<"Add a number of the first border"<<endl;
            cin>>a;
            cout<<"Add a number of the second border"<<endl;
            cin>>b;
            cout<<"Add a number of all objects"<<endl;
            cin>>n;
            cout<<"Add a probability"<<endl;
            cin>>p;
            cout<<"The probability is "<<int_moivre_laplace(a,b,n,p)<<endl;
            cout<<"***"<<endl;
            start();
        }
    else if(n==5)
        {
            int m;
            float p;
            cout<<"Add a number of all objects"<<endl;
            cin>>n;
            cout<<"Add a number of the set"<<endl;
            cin>>m;
            cout<<"Add a probability"<<endl;
            cin>>p;
            cout<<"The probability is "<<loc_moivre_laplace(n,m, p)<<endl;
            cout<<"***"<<endl;
            start();
        }
    else if(n==6)
        {
            double x;
            cout<<"Add a number "<<endl;
            cin>>x;

            cout<<"The meaning is "<<phi(x)<<endl;
            cout<<"***"<<endl;
            start();
        }
    else if(n==7)
        {
            int n;
            string s;
            cout<<"Add a number of elements ";
            cin>>n;
            discred_random_variable drv_1 (n);
            cin>>drv_1;
            cout<<drv_1.expected_value();
            start();
        }
    else if(n==8)
        {
            int n;
            string s;
            cout<<"Add a number of elements ";
            cin>>n;
            discred_random_variable drv_1 (n);
            cin>>drv_1;
            cout<<drv_1.variance();
            start();
        }
    else if(n==9)
        {
            int n;
            string s;
            cout<<"Add a number of elements ";
            cin>>n;
            discred_random_variable drv_1 (n);
            cin.get();
            cin>>drv_1;
            cout<<drv_1.standard_deviation();
            start();
        }
    else if(n==0){cout<<"***Program ended***"<<endl;}
    else {cout<<"Wrong addition"; start();} 
    cout<<endl;
    
}
int main()
{
    start();
    getch();
    return 0;
}