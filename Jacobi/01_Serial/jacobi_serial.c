#include <math.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#define NN 4096
#define NM 4096

double A[NN][NM];
double Anew[NN][NM];

struct timeval timerStart;

void StartTimer()
{
    gettimeofday(&timerStart, NULL);
}

// time elapsed in ms
double GetTimer()
{
  struct timeval timerStop, timerElapsed;
  gettimeofday(&timerStop, NULL);
  timersub(&timerStop, &timerStart, &timerElapsed);
  return timerElapsed.tv_sec*1000.0+timerElapsed.tv_usec/1000.0;
}

int main(int argc, char** argv)
{
  const int n = NN;
  const int m = NM;
  const int iter_max = 1000;
  
  const double tol = 1.0e-6;
  double error     = 1.0;
  
  memset(A, 0, n * m * sizeof(double));
  memset(Anew, 0, n * m * sizeof(double));
  
  for (int j = 0; j < n; j++)
    {
      A[j][0]    = 1.0;
      Anew[j][0] = 1.0;
    }
  
  printf("Jacobi relaxation Calculation: %d x %d mesh\n", n, m);
  
  StartTimer();
  int iter = 0;
  
  while ( error > tol && iter < iter_max )
    {
      error = 0.0;
      
      for( int j = 1; j < n-1; j++)
        {   
	  for( int i = 1; i < m-1; i++ )
            {   
	      Anew[j][i] = 0.25 * ( A[j][i+1] + A[j][i-1]
                                    + A[j-1][i] + A[j+1][i]);
	      error = fmax( error, fabs(Anew[j][i] - A[j][i]));
            }
        }
      
      for( int j = 1; j < n-1; j++)
        {   
	  for( int i = 1; i < m-1; i++ )
            {   
	      A[j][i] = Anew[j][i];
            }
        }
      
      if(iter % 100 == 0) printf("%5d, %0.6f\n", iter, error);
      
      iter++;
    }
  
  double runtime = GetTimer();
  
  printf(" total: %f s\n", runtime / 1000);
}
