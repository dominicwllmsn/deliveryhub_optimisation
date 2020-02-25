// Final Assignment
// Read in the GBPlaces.csv and use an optimization algorithm over many
// iterations to find the ideal place to build a delivery hub.

// includes and definitions
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>

#define PI 3.14159

int f_evals = 0; // counter for the number of function evaluations

using namespace std;

// Function declarations
double greatCircle(double lat1, double long1, double lat2, double long2);
double random_number ( double lower, double upper, int n );


// Main program

int main() {
// READ IN FILE
    // declare the required arrays
    double population[100];
    double latitude[100];
    double longitude[100];

    // open the input file as an input file stream (ifstream)
    ifstream dataFile( "/Users/dominicwilliamson/Documents/Programming/UoMCourse/C++/UoM_CPP/Week10_11/GBplaces.csv");

    // if the file is open...
    if (dataFile.is_open()) {
        cout << "GBplaces.csv has been opened.\n";
        // define a variable to hold the input line 
        string line;
        int row = -1;

        // while there are still lines to read from input
        while (getline(dataFile, line))
        {   
            // skip the first row (which is just the headings)
            if (row == -1) {
                row++;
            } else {
                // define input stringstream object iss
                // the use of the input stringstream creates a buffer (like cin) 
                // which can be used in the getline() function (to split up a string)
                istringstream iss;
                iss.str(line); // assign to iss the line of the .csv we want to split up
                string field; // create field to store each substring from getline
                int column = 1;

                // while there are substrings in the buffer, we can split up about the ',' delimiter char
                while (getline(iss, field, ','))
                {   
                    // skip the first two columns - we only need the pop, lat and long data
                    if (column < 3) {
                        column++;
                    } else {
                        // define element as a double of the field str and insert into the appropriate array 
                        double element = atof(field.c_str());
                        if (column == 3) {
                            population[row] = element;
                        } else if (column == 4) {
                            latitude[row] = element;
                        } else if (column == 5) {
                            longitude[row] = element;
                        }
                        column++; // move onto the next column, if available
                    }
                }
                // printf("(%5.2f,%5.2f) -> %7.0f \n", latitude[row],longitude[row],population[row]);
                row++; // move onto the next row, if available
            }
        }    

    // now at the file end, so close it
    dataFile.close();
    cout << "Input file read and closed" << endl;
  } else {
    // if the file couldn't be opened, write a message and end main()
    cout << "Unable to open input file" << endl;
    return 1;
  }

  // Find minimum and maximum of long, lat and pop 
  // need long,lat extremes for the random_number function bounds, and pop extremes
  // to normalize the pop array in order to 'weight' each distance
  double minLong = longitude[0];
  double maxLong = longitude[0];
  double minLat = latitude[0];
  double maxLat = latitude[0];
  double minPop = population[0];
  double maxPop = population[0];

  for (int i = 1; i < 100; i++) {
    if (longitude[i] > maxLong)
      maxLong = longitude[i];
    if (longitude[i] < minLong)
      minLong = longitude[i];
    if (latitude[i] > maxLat)
      maxLat = latitude[i];
    if (latitude[i] < minLat)
      minLat = latitude[i];
    if (population[i] > maxPop)
      maxPop = population[i];
    if (population[i] < minPop)
      minPop = population[i]; 
  }

  // cout << minLong << "," << minLat << " ___ " << maxLong << "," << maxLat << endl;


// OPTIMIZATION
  // declare variables
  int dx, dy;
  double x, y; // holds the current best values of long and lat respectively.
  double step = 0.01; // step size to move in x and y
  // note in a more sophisticated program you could vary the step size
  // depending on whether the minimum was reached or not and try to 
  // make the method more efficient and more accurate like that
  double value, oldValue, newValue, minValue;
  int N = 200; // number of iterations of for-loop
  double globalX = 0; // current global minimum of long
  double globalY = 0; // current global minimum of lat
  double globalValue = 1e9; // current global minimum of distance - set to a large number


  srand(time(NULL)); // seeds random number generator

  for (int n = 0; n < N; n++)
  {
  
  // pick a starting point at random between min & max long/lat values
  x = random_number(minLong, maxLong, 100);
  y = random_number(minLat, maxLat, 100);
  value = 0;

  // work out the distance sum at the random location 
  for (int i = 0; i < 100; i++) {
    double distance = greatCircle(y, x, latitude[i], longitude[i]);
    value += distance;
  }
  // cout << "Initial sum: " << value << endl;
  // main do-loop to continually 'improve' x, y
  do {
    
    // save the current value
    oldValue = value;
    minValue = oldValue; // set the minValue for the local search to be the current value
    
    // now look around the current point to see if there's a better one nearby
    for ( int i = -1; i <= 1; i++ ) {
      for ( int j = -1; j <= 1; j++ ) {
    // this gives 9 points including the current point (when i=0, j=0)
    if ( i==0 && j==0 ) {
    // so maybe you want to miss that one and save some function evaluations
    } else {
      // values at neighbouring point
      x = x + step * i;
      y = y + step * j; 
      // calculate the value for the sum of distances at the new coordinates
      newValue = 0;
      for (int n = 0; n < 100; n++) {
        double distance = greatCircle(y, x, latitude[n], longitude[n]);
        newValue += distance;
      }
      if ( newValue <= minValue ) { // is it smaller than minValue?
        // if yes, set minValue and save point i,j values
        dx = i; 
        dy = j;
        minValue = newValue;
      }
    }
      }
    }
    
    // update x and y to new point with smaller 'value'
    x += step * dx;
    y += step * dy;
    value = minValue;
    // cout << "Subsequent sum: " << value << endl;

  } while ( value < oldValue ); // repeat all this while we can find a greater value than the previous one
  
  if (value < globalValue){
    globalValue = value;
    globalX = x;
    globalY = y;
  }

  }

  cout << "Function evaluations: " << f_evals << "\n"; // write out the total number of function evaluations
  cout << "\"Global\" minimum is at (" << globalY << " , " << globalX << ") -> " << globalValue;


  return 0;
}

// Function definitions

double greatCircle(double lat1, double long1, double lat2, double long2) {
        /*
        Inputs: latitudes and longitudes of two locations on Earth
        Output: 'Great Circle' distance between the two locations in meters
        */
        lat1 = lat1*PI/180.;
        long1 = long1*PI/180.;
        lat2 = lat2*PI/180.;
        long2 = long2*PI/180.;
        double dLat = lat2 - lat1;
        double dLong = long2 - long1;
        double R = 6371e3;
        double a = pow(sin(dLat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dLong/2),2);
        double c = 2 * atan2((pow(a,0.5)), pow((1-a),0.5));

        return R * c;
}


double random_number(double lower, double upper, int n) {
  /*
  Inputs: lower and upper bounds for the range of the random number 
  n is the number of segments the range is split into
  Output: pseudo-random number in the specified range 
  */ 
  double r;
  r = lower + (rand() % (n + 1) * (1./n) * (upper-lower));
  return r;
}
