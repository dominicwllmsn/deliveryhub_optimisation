// hill
// demonstration of hillclimbing as an optimization method

// includes and definitions

#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>

#define PI 3.14159

int f_evals = 0; // keep a count function evaluations - note, global variable

using namespace std;

// functions

double fitness ( double x, double y ) {
  // function to return the 'value' of thing being optimized
  // this is the function which tells your program how good the current solution is
  double f;
  f = cos(x/3) * cos(y/3) * pow(cos(2.*x),2) * pow(cos(2.*y),2);
  f_evals++;
  return f;
}

double random_number ( double lower, double upper, int n ) {
  // function to return a random number between two limits, lower and upper
  // n is the amount of bits to split the range into
  double r;
  r = lower + (rand() % (n + 1) * (1./n) * (upper-lower));
  return r;
}

// main program

int main() {
  // declare variables
  int dx, dy;
  double x, y; // holds the current best values of x and y
  double step = 0.001; // step size to move in x and y
  // note in a more sophisticated program you could vary the step size
  // depending on whether the maximum was reached or not and try to 
  // make the method more efficient and more accurate like that
  double value, oldValue, newValue, maxValue;
  int N = 50;
  double globalValue = 0;
  double glx;
  double gly;

  srand(time(NULL)); // seeds random number generator

  for (int n = 0; n < N; n++)
  {
  
  // pick a starting point at random between -pi/4 and +pi/4
  x = random_number ( -5., 5., 100 ) ;
  y = random_number ( -5., 5., 100 ) ;

  // work out value of function at that point - how 'good' is point x,y?
  value = fitness(x,y);

  // main do-loop to continually 'improve' x, y
  do {
    
    // save the current value
    oldValue = value;
    maxValue = oldValue; // set the maxValue for the local search to be the current value
    
    // now look around the current point to see if there's a better one nearby
    for ( int i = -1; i <= 1; i++ ) {
      for ( int j = -1; j <= 1; j++ ) {
    // this gives 9 points including the current point (when i=0, j=0)
    if ( i==0 && j==0 ) {
    // so maybe you want to miss that one and save some function evaluations
    } else {
      newValue = fitness(x + step * i, y + step * j); // value at neighbouring point
      if ( newValue >= maxValue ) { // is it bigger than maxValue?
        // yes so set maxValue and save point i,j values
        dx = i; 
        dy = j;
        maxValue = newValue;
      }
    }
      }
    }
    
    // update x and y to new point with higher 'value'
    x += step * dx;
    y += step * dy;
    value = maxValue;
    

  } while ( value > oldValue ); // repeat all this while we can find a greater value than the previous one
  
  if (value > globalValue){
    globalValue = value;
    glx = x;
    gly = y;
  }

  }

  cout << "Function evaluations: " << f_evals << "\n"; // write out the total number of function evaluations
  cout << "\"Global\" maximum is at (" <<x << " , " <<y << ") -> " << globalValue;
  return 0;
}