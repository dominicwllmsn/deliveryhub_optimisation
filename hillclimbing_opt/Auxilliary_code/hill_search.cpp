// hill_search
// version of hill that searches randomly within the parameter space and climbs

#include <iostream>
#include <cmath>
#include <ctime>

#define PI 3.14159

using namespace std;

int fevals = 0;

double random_number ( double upper, double lower, int n ) {
  double r;
  r = lower + ( rand() % (n+1) * (1./n) * (upper-lower));
  return r;
}

double fitness ( double x, double y ) {
  double f;
  f = cos(x/3) * cos(y/3) * pow(cos(2.*x),2) * pow(cos(2.*y),2); 
  fevals++;
  return f;
}

int main() {

  double x,y,glx,gly;
  int dx,dy;
  double step = 0.001;
  double value, newValue, oldValue;
  double globalMax = 0;

  srand(time(NULL));


  //  cout << x << " " << y << " -> " << value << "\n";

  for ( int k = 0; k < 50; k++ ) {

  // first pick a random starting point for x and y
  x = random_number ( PI, -PI, 100);
  y = random_number ( PI, -PI, 100);

  // now work out the value of the function at point x,y
  value = fitness(x,y);
    
  do {
  // now look around the current point and see if we can go somewhere where value is higher
  oldValue = value;

  for ( int i = -1; i <= 1; i++ ) {
    for ( int j = -1; j <= 1; j++ ) {
      if ( i == 0 && j == 0 ) {
      } else {
    newValue = fitness( x + step * i, y + step * j );
    if ( newValue >= value ) {
      dx = i;
      dy = j;
      value = newValue;
    }
      }
    }
  }

  // update to new position and new value
  //value = maxValue;
  x += step * dx;
  y += step * dy;

  } while ( value > oldValue );
  
  if(value > globalMax) {
    globalMax = value;
    glx = x;
    gly = y;
  }

  }

  cout << glx << " " << gly << " -> " << globalMax << "\n";

  cout << fevals << "\n";

  return 0;
}