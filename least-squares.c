#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// Does a least squares fit given points x, and y, and errors
//  prints out chi^2, int, slope, error in a, error in b 

#define f(i) (slope*x[i] + intercept)
#define sqr(x) ((x)*(x))
#define s dely
#define N 63 //how many data points?
#define n 6 //what degree of polynomial? 

void gaussLS(int m, int n, double a[m][n], double x[n-1]){
    int i,j,k;
    for(i=0;i<m-1;i++){
        //This part does the Partial Pivoting
        for(k=i+1;k<m;k++){
            //If this case has diagonal element that is smaller than any of the terms below it
            if(fabs(a[i][i])<fabs(a[k][i])){
                //Exchange rows
                for(j=0;j<n;j++){                
                    double tempMatrix;
                    tempMatrix=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=tempMatrix;
                }
            }
        }
        //Start the Gauss Elimination
        for(k=i+1;k<m;k++){
            double element=a[k][i]/ a[i][i];
            for(j=0;j<n;j++){
                a[k][j]=a[k][j]-element*a[i][j];
            }
        }  
    }
    //Begin the Back-substitution
    for(i=m-1;i>=0;i--){
        x[i]=a[i][n-1];
        for(j=i+1;j<n-1;j++){
            x[i]=x[i]-a[i][j]*x[j];
        }
        x[i]=x[i]/a[i][i];
    }       
}

void showMatrix(int a, int b, double matrix[a][b]){
    int i,j;
    for(i=0;i<a;i++){
        for(j=0;j<b;j++){
            printf("%lf\t",matrix[i][j]);
        }
        printf("\n");
    } 
}
int main(){
  
  double x[N], y[N], delx[N], dely[N];
  double del, sumofx, sumofy,sumofxx,sumofyy, sumofxy, una, unb;
  double  sumofxxs, sumofxs, sumofss, sumofxys, sumofys;
  double  slope, ave, intercept, chi;
  FILE *datafile;
  datafile=fopen("rowdata.dat","r");
  sumofx=sumofy=sumofxxs=sumofxs=una=unb=0;
  sumofss=sumofxys=sumofys=0;
  

  for(int i=0;i<N;i++)
    fscanf(datafile, "%lf %lf %lf %lf\n", 
	   &x[i],&y[i],&dely[i],&delx[i]);


  for(i=0;i<N;i++){
    s[i]=dely[i];
    sumofx   += x[i];
    sumofy   += y[i]; 
    sumofxx  += x[i]*x[i];
    sumofyy  += y[i]*y[i];
    sumofxy  += x[i]*y[i];
    sumofxxs += (x[i]*x[i])*s[i];
    sumofxs  += x[i]*s[i];
    sumofss  += s[i];
    sumofxys +=(x[i]*y[i])*s[i];
    sumofys  +=(y[i])*s[i];
  }

  del=sumofss*sumofxxs-sumofxs*sumofxs;
  ave=N*sumofxx-sumofx*sumofx;
  slope=(sumofss*sumofxys-sumofxs*sumofys)/del;
  intercept=(sumofxxs*sumofys-sumofxs*sumofxys)/del;
  una=sqrt((1/del)*(sumofxxs));
  unb=sqrt((1/del)*(sumofss));

  for(chi=0,i=0;i<N;i++)
     chi += sqr(y[i] - f(i))*s[i];
    double X[2*n+1];  //this is for the chi squared for higher order we need to matrix it again
    for(i=0;i<=2*n;i++){
        X[i]=0;
        for(j=0;j<N;j++){
            X[i]=X[i]+pow(x[j],i);
        }
    }
    //the normal augmented matrix for our calculations
    double B[n+1][n+2];  
    // the right hand side
    double Y[n+1];      
    for(i=0;i<=n;i++){
        Y[i]=0;
        for(j=0;j<N;j++){
            Y[i]=Y[i]+pow(x[j],i)*y[j];
        }
    }
    for(i=0;i<=n;i++){
        for(j=0;j<=n;j++){
            B[i][j]=X[i+j]; 
        }
    }
    for(i=0;i<=n;i++){
        B[i][n+1]=Y[i];
    }
    double A[n+1];
    printf("The best fit equation:\n");
    showMatrix(n+1,n+2,B);
    gaussLS(n+1,n+2,B,A);
    for(i=0;i<=n;i++){
        printf("%lfx^%d+",A[i],i);//final output of program to be used in the next program. save this value. 
    }

  fclose(datafile);
  printf("chi= %lf\t a= %e.9\t b= %e.9\t una= %lf\t unb= %lf\n", 
	 chi, slope, intercept, unb, una);
}
