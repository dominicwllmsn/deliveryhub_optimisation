/* 
  Final Assignment 
  Dominic Williamson - ID: 9846595
  version 1.0 : STANDARD PROBLEM
  Program reads in the GBPlaces.csv and use an optimization algorithm over many
  iterations to find the ideal place to build a delivery hub. The Great 
  Circle formula is used to find the distance between two coordinates.
*/

// includes and definitions
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>
#include <vector>
#include <string>

#define PI 3.14159

int f_evals = 0; // counter for the number of function evaluations

using namespace std;

// Function definitions

double random_number(double lower, double upper, int n) {
  /*
  Inputs: lower and upper bounds for the range of the random number, 
  n is the number of segments the range is split into.
  Output: pseudo-random number in the specified range.
  */ 
  return lower + (rand() % (n + 1) * (1./n) * (upper-lower));
}

double greatCircle(double lat1, double long1, double lat2, double long2) {
        /*
        Inputs: latitudes and longitudes of two locations on Earth
        Output: the Great Circle distance between the two locations in kilometers
        */
        f_evals++;
        lat1 = lat1*PI/180.;
        long1 = long1*PI/180.;
        lat2 = lat2*PI/180.;
        long2 = long2*PI/180.;
        double dLat = lat2 - lat1;
        double dLong = long2 - long1;
        double R = 6371;
        double a = pow(sin(dLat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dLong/2),2);
        double c = 2 * atan2((pow(a,0.5)), pow((1-a),0.5));

        return R * c;
}

// Main program

int main() {
// READ IN FILE
    // declare the required vectors
    vector< vector<double> > coords;
    vector<double> temp;

    // open the input file as an input file stream (ifstream)
    ifstream dataFile( "/Users/dominicwilliamson/Documents/Programming/UoMCourse/C++/UoM_CPP/Week10_11/GBplaces.csv");

    // if the file is open...
    if (dataFile.is_open()) {
        cout << "Input file has been opened.\n";
        // define a variable to hold the input line 
        string line;
        // while there are still lines to read from the input
        while (getline(dataFile, line))
        {   
            // skip the first row (which is just the headings)
            if (line.length() > 0 && line.find('%') != 0) {
                // define an input stringstream object called iss
                // N.B. using the input stringstream creates a buffer (like cin) 
                // which can be used in the getline() function to split up a string
                istringstream iss;
                iss.str(line); // assign to iss the line of the .csv we want to split up
                string field; // create field to store each substring from getline
                int column = 1;
                vector<double> temp; // clear 'temp' when reading in a new row

                // while there are substrings in the buffer, we can split up about the ',' delimiter char
                while (getline(iss, field, ','))
                {   
                    // skip the first three columns - we only need the latitude and longitude data
                    if (column < 4) {
                        column++;
                    } else {
                        // define element as a double of the field str and push back onto temp
                        double element = atof(field.c_str());
                        temp.push_back(element);
                        column++; // move onto the next column, if available
                    }
                }
                coords.push_back(temp);
            }
        }    

    // we are now at the end of the file, so close it
    dataFile.close();
    cout << "Input file has been read and closed." << endl;
  } else {
    // if the file couldn't be opened, write a message and return main() 
    // with an error message (i.e. a non-zero return value)
    cout << "Unable to open input file" << endl;
    return 1;
  }


// Find the minimum and maximum of long & lat to use as the bounds for 
// the random_number function 
  double minLat = coords[0][0];
  double maxLat = coords[0][0];
  double minLong = coords[0][1];
  double maxLong = coords[0][1];


  for (int i = 1; i < coords.size(); i++) {
    if (coords[i][0] > maxLat)
      maxLat = coords[i][0];
    if (coords[i][0] < minLat)
      minLat = coords[i][0];
    if (coords[i][1] > maxLong)
      maxLong = coords[i][1];
    if (coords[i][1] < minLong)
      minLong = coords[i][1];
  }


// OPTIMIZATION
  // declare variables
  int dx, dy;
  double x, y; // holds the current best values of long and lat respectively.
  double step = 0.01; // step size to move in x and y
  double value, oldValue, newValue, minValue; // variables required to hold each distance sum
  int N = 100; // number of iterations of the for-loop (i.e. number of random coordinates tested)
  double globalX, globalY; // current global minimum of long and lat respectively.
  double globalValue = 1e10; // current global minimum of distance - initially set to a v.large number

  srand(time(NULL)); // seeds the random number generator

  for (int n = 0; n < N; n++)
  {

  // pick a starting point at a random point between min & max long/lat values
  x = random_number(minLong, maxLong, 100);
  y = random_number(minLat, maxLat, 100);
  value = 0; // 'value' will hold the sum of all of the distances calculated

  // work out the distance sum at the random location 
  for (int i = 0; i < coords.size(); i++) {
    double distance = greatCircle(y, x, coords[i][0], coords[i][1]);
    value += 2*distance;
  }

  // main do-loop to search for a local minimum of the summed distance
  do {

    oldValue = value;     // save the current value
    minValue = oldValue; // set the minValue for the local search to be the current value
    
    // now look around the current point to see if there's a better one nearby
    // this gives 8 points excluding the hub itself (when i=0, j=0)
    for ( int i = -1; i <= 1; i++ ) {
      for ( int j = -1; j <= 1; j++ ) {
    if ( i==0 && j==0 ) {
      continue;
    } else {
      // find values at a neighbouring point
      x = x + step * i;
      y = y + step * j; 

      // calculate the value for the sum of distances at the new coordinates
      newValue = 0;
      for (int m = 0; m < coords.size(); m++) {
        double distance = greatCircle(y, x, coords[m][0], coords[m][1]);
        newValue += 2*distance;
      }

      // if the newValue is smaller than the current minimum, set the new minValue
      // and save the i,j values
      if ( newValue <= minValue ) { 
        dx = i; 
        dy = j;
        minValue = newValue;
      }
    }
      }
    }
    
    // update x and y to new point with smaller summed distance value
    x += step * dx;
    y += step * dy;
    value = minValue;

  } while ( value < oldValue ); // repeat all this while we can find a smaller value than the previous one
  
  // Test to see if the local minimum we just calculated is smaller than the current smallest known minimum
  // if yes: set globalValue to this newly found value, as well as the coordinates of this new minimum.
  if (value < globalValue){
    globalValue = value;
    globalX = x;
    globalY = y;
  }

  }

  // Average distance from hub to each location. Divide by 100 since 100 locations in GBplaces. Divide  
  // by 2 since globalValue is found by calculating the distance from the hub to the location and back.
  double avgDist = (globalValue*0.6214/100)/2;

  // Print the total number of function evaluations and the 'global' minimum 
  cout << "******************************" << endl;
  cout << "Function evaluations: " << f_evals << "\n"; 
  printf("- The delivery hub should be placed at (%.4f, %.4f).\n",globalY, globalX);
  printf("- A hub at these coordinates would be an average distance of %.2f miles away from each \
town/city evaluated.\n", avgDist);

  return 0;
}

