#include <iostream>
#include <cmath>
#include <string> 
#include <cstdlib>
#define _USE_MATH_DEFINES
using namespace std;


int comb(int a, int b)
{
    int i;
    long long c{1},f;
    if (b<a-b)
        {
            for (i=1; i<=b; i++)
            {
                c*=i;
            }
            b=a-b;
        }
    else
        {
            for (i=1; i<=a-b; i++)
            {
                c*=i;
            }
        }
    f=c;
    c=1;
    for (i=b+1; i<=a; i++)
    {
        c*=i;
    }
    c=c/f;
    return c;
}
long long discomb(int n)
{
int f{1};
double s{1};
for (int i{1}; i<=n;i++)
    {
        f*=i;
        s+=(pow(-1,i)/f);
    }
s*=f;
return s;
}
double bernulli(int n, int k, float p)
{
    double b;
    b=comb(n,k)*pow(p,k)*pow(1-p,n-k);
    return b;
}
double phi(double x)
    {
        double phi;
        phi=(erf(x/M_SQRT2))/2;
        return phi;
    }
double int_moivre_laplace(double a, double b, int n, double p)
{
    double ml, np{n*p},npq;
    npq=sqrt((1-p)*np);
    ml=phi((b-np)/npq)-phi((a-np)/npq);
    return ml;
}
double loc_moivre_laplace(int n, int m, double p)
{
    double ml, np{n*p},npq;
    npq=sqrt((1-p)*np);
    ml=((1/sqrt(2*M_PI))*pow(M_E,(-(pow(((m-np)/npq), 2))/2)))/npq;
    return ml;
}
double max_str(string s)
{
    int i,n,l;
    double e, max;
    string ns{""};
    l=size(s);
    n=s.find(" ",0);
    max=stod(s.substr(0,n));
    i=n;
    while (i<l)
    {
        i++;
        if(s[i]!=' ') ns+=s[i];
        if(s[i]==' ') 
        {
            e=stod(ns);
            ns="";
            
        }
        if (max<e) max=e;
    }
    return max;
}



class discred_random_variable
{
private:
  int number;
  double **cells;
  void init ()
  {
    int i,j;
      //destroy ();
    cells = new double *[2];
    for (i = 0; i < 2; i++)
      {
        cells[i] = new double[number];
        for (j = 0; j < number; j++)
          {
            cells[i][j] = 0;
          }
      }
  }
  void destroy ()
        {
        int i;
        if (cells)
        {
        delete cells;
        }
        cells = nullptr;
        }

public:

discred_random_variable (int n)
    {
        number=n;
        init();
    }
void demo ()
    {
        number=3;
        init();
        int j;
        for (j = 0; j < number; j++)
          {
            cells[0][j] = j;
          }
          cells[1][0] = 0.1;
          cells[1][1] = 0.3;
          cells[1][2] = 0.6;
    }

void sort ()
  {
    int i,j,k,c;
    double max{cells[0][0]},tmp1,tmp2;
    for (i = 1; i < number; i++)
        {
          if (max<cells[0][i]) max=cells[0][i];
        }
        
    for (i = 0; i < number; i++)
        {
        for (j = i+1; j < number; j++)
            {
               
                if (cells[0][i]==cells[0][j]) 
                {    
                    cells[1][i]+=cells[1][j];
                    cells[1][j]=-1;
                    cells[0][j]=max+1;
                }
            }
        }
    for (k = number / 2; k > 0; k /= 2)
        for (i = k; i < number; i++)
        {
            tmp1 = cells[0][i];
            tmp2 = cells[1][i];
            for (j = i; j >= k; j -= k)
            {
                if (tmp1 < cells[0][j - k])
                   {
                    cells[0][j] = cells[0][j - k];
                    cells[1][j] = cells[1][j - k];
                   }
                else break;
            }
            cells[0][j] = tmp1;
            cells[1][j] = tmp2;
        }

      for (i=0; i<number; i++)
      {
        if (cells[1][i]<0) {number=i; break;}
      }
    double** a;  
    a= new double *[2];
    a[0] = new double [number] ;
    a[1] = new double [number] ;
    for (i=0; i<number; i++)
      {
        a[0][i]=cells[0][i];
        a[1][i]=cells[1][i];
      }
    cells=nullptr;
    cells=a;
  }
  friend ostream & operator<< (ostream & stream, const discred_random_variable & drv)
  {
    int j;
	for (j = 0; j < drv.number; j++)
	  {
	    stream <<"variable  " << drv.cells[0][j] << " probability " << drv.cells[1][j]<< endl;
	  }
    return stream;
  }

friend istream & operator>> (istream & in, discred_random_variable & drv)
    {
      int i;
      double s{0};
      for (i = 0; i < drv.number; i++)
      {
  
            if (s<=1.00)
            {
                cout<<"Add variable ";
                in>>drv.cells[0][i];
                cout<<"Add probability ";
                in>>drv.cells[1][i];
                s+=drv.cells[1][i];
            }
            else
            {
              cout<<"Error: a sum of the probabilities is over 1. Add again ";
              s=0;
              i=-1;
            }
        }
        return in;
    }
discred_random_variable *operator = (const discred_random_variable & drv)
  {
    destroy ();
    number = drv.number;
    init ();
    int i, j;
    for (i = 0; i < 2; i++)
      {
	for (j = 0; j < number; j++)
	  {
	    cells[i][j] = drv.cells[i][j];
	  }
      }
    return this;
  }
  discred_random_variable operator + (const discred_random_variable & drv)
  {
    int Cn{number*drv.number};
    discred_random_variable C (Cn);   
	int i,j,k{0};
    for (i = 0; i < number; i++)
	  {
	    for (j = 0; j < drv.number; j++)
	      {
            C.cells[0][k] =(cells[0][i] + drv.cells[0][j]);
            C.cells[1][k] =(cells[1][i] * drv.cells[1][j]);
            k++;
	      }
	  }
    C.sort();
    return C;
  }

discred_random_variable operator * (const discred_random_variable & drv)
  {
  int Cn{number*drv.number};
  discred_random_variable C (Cn);   
	int i,j,k{0};
    for (i = 0; i < number; i++)
	  {
	    for (j = 0; j < drv.number; j++)
	      {
            C.cells[0][k] =(cells[0][i] * drv.cells[0][j]);
            C.cells[1][k] =(cells[1][i] * drv.cells[1][j]);
            k++;
	      }
	  }
    C.sort();
    return C;
  }
discred_random_variable operator - (const discred_random_variable & drv)
  {
  int Cn{number*drv.number};
  discred_random_variable C (Cn);   
	int i,j,k{0};
    for (i = 0; i < number; i++)
	  {
	    for (j = 0; j < drv.number; j++)
	      {
            C.cells[0][k] =(cells[0][i] - drv.cells[0][j]);
            C.cells[1][k] =(cells[1][i] * drv.cells[1][j]);
            k++;
	      }
	  }
    C.sort();
    return C;
  }
void get_DRV_by_String(string s)//Не работает пока
  {
    int i{0},j{0},l;
    double x;
    bool b{true};
    string ns{""};
    s.insert(size(s),"  ");
    s.insert(0," ");
    l=size(s);
    for (i=0;i<l; i++)
      {
        if(s[i]!=' ')
        { 
          if (s[i]==',') ns+='.'; 
          else ns+=s[i];
        }
        if((s[i]==' ')and(ns!=""))
        {
          x=stod(ns);
          ns="";
          if(b) 
            {
              cells[0][j]=x; b=false;
            } 
          else 
            {
              cells[1][j]=x;  j++; b=true;
            };
        }
      }
      sort();
  }
double expected_value()
  {int i; double s{0};
  for (i=0; i<number; i++)
    {
      s+=(cells[0][i]*cells[1][i]);
    }
    s/=number;
    return s;
  }
double variance()
  {int i; double s1{0}, s2{0};
  for (i=0; i<number; i++)
    {
      s1+=(cells[0][i]*cells[0][i]*cells[1][i]);
      s2+=(cells[0][i]*cells[1][i]);
    }
    s1=(s1-(s2*s2))/number;
    return s1;
  }
double standard_deviation()
  {int i; double s1{0}, s2{0};
  for (i=0; i<number; i++)
    {
      s1+=(cells[0][i]*cells[0][i]*cells[1][i]);
      s2+=(cells[0][i]*cells[1][i]);
    }
    s1=sqrt((s1-(s2*s2))/number);
    return s1;
  }
};

class function_random_variable
{
private:
  int number;
  double *cells;
  void init ()
  {
    cells = new double[number];
    for (int i{0}; i < number; i++)
      {
        cells[i] = i;
      }
  }
  void destroy()
  {
    delete [] cells;
    cells=nullptr;
  }
public:
  function_random_variable(int n)
  {
    number=n;
    init();
  }
};