/*
  Final Assignment 
  Dominic Williamson - ID: 9846595
  version 1.1 : TWO HUBS + STEP OPTIMIZATION
  Program reads in the GBPlaces.csv and uses an optimization algorithm over many
  iterations to find the ideal place to build two delivery hubs. 
  The step length for this algorithm is optimized by decreasing it by a factor of 
  ten over several iterations, allowing for more efficient and more accurate
  determination of the global minima. The Great Circle formula is used to find
  the distance between two coordinates.
*/

// includes and definitions
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>
#include <vector>

#define PI 3.14159
#define RADIUS 6371 // Earth's radius in kilometers

int f_evals = 0; // counter for the number of function evaluations

using namespace std;

// Function definitions

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

double fitness(double lat1, double long1, double lat2, double long2, vector< vector<double> > data) {
  /*
    Inputs: latitudes and longitudes of two locations on Earth, a 2D vector of location coordinates
    Output: the (minimum) sum of the distances from the hubs to each location, using the Great Circle
            distance formula, in kilometers
  */

  // Declare the necessary variables to find the Great Circle distance
  double dLat1, dLong1, a1, c1;
  double dLat2, dLong2, a2, c2;
  // Start with f, the sum variable, being equal to 0
  double f = 0;

  // Convert test coordinates to radians
  lat1 = (PI*lat1) / 180;
  long1 = (PI*long1) / 180;
  lat2 = (PI*lat2) / 180;
  long2 = (PI*long2) / 180;

  // Loop over all the latitudes and longitudes in the 2D vector
  for (int j = 0; j < data.size(); j++) {

    // Convert the latitudes and longitudes to radians
    data[j][0] = (PI*data[j][0]) / 180;
    data[j][1] = (PI*data[j][1]) / 180;

    // Do the calculations via the Haversine formula for the first hub
    dLat1 = lat1 - data[j][0];
    dLong1 = long1 - data[j][1];
    a1 = (pow(sin(dLat1 / 2), 2)) + (cos(data[j][0]) * cos(lat1) * pow(sin(dLong1 / 2), 2));
    c1 = 2 * atan2(sqrt(a1), sqrt(1 - a1));

    // Do the calculations via the Haversine formula for the second hub
    dLat2 = lat2 - data[j][0];
    dLong2 = long2 - data[j][1];
    a2 = (pow(sin(dLat2 / 2), 2)) + (cos(data[j][0]) * cos(lat2) * pow(sin(dLong2 / 2), 2));
    c2 = 2 * atan2(sqrt(a2), sqrt(1 - a2));

    if (c1 <= c2) {
      // f is multiplied by 2 because we want to know the distances from hub AND back
      f += 2 * RADIUS * c1;
    } else {
      f += 2 * RADIUS * c2;

    }
  }

  f_evals++;
  return f;

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
// declare variables to hill climb (or descend?) to each minima
  int dx1, dy1, dx2, dy2;
  double value, oldValue, newValue, minValue; // // variables required to hold each distance sum
  double lathub1, longhub1, lathub2, longhub2;
  int N = 50; // number of iterations of the for-loop (i.e. number of random coordinates tested)  
  double globalY1, globalX1, globalY2, globalX2; // current global minimum of each lat/long value
  double globalValue = 1e10; // current global minimum of distance - set initially to a v.large number

  srand(time(NULL)); // seeds the random number generator

  for (int n = 0; n < N; n++) {
  // pick the first hub's starting point at a random point between the min & max long/lat values
  lathub1 = random_number(minLat, maxLat, 100);
  longhub1 = random_number(minLong, maxLong, 100);

  // pick the second hub's starting point at a random point between the min & max long/lat values
  lathub2 = random_number(minLat, maxLat, 100);
  longhub2 = random_number(minLong, maxLong, 100);

  // Using the fitness function we can calculate the total distance travelled
  value = fitness(lathub1, longhub1, lathub2, longhub2, coords);

  double step = 0.1; // step size to move in x and y

  // for loop surrounding the hill climbing method to reduce the step size each iteration
  for (int p = 0; p < 5; p++) {
    step /= pow(10,p);

  // main do-loop to search for a local minimum of the summed distance
  do {

    oldValue = value; // Save the current value
    minValue = oldValue; // set the minValue for the local search to be the current value

    // now look around the current point to see if there's a better one nearby
    // this gives 8 points about each hub excluding the hub itself (when i=0, j=0 or k = 0, l = 0)    
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          for (int l = -1; l <= 1; l++) {

            if (i == 0 && j == 0 && k == 0 && l == 0) {
              continue;
            } else {
              // Value at neighbouring point
              newValue = fitness(lathub1 + step * i, longhub1 + step * j, lathub2 + step * k, longhub2 + step * l, coords);
              // Is the newvalue smaller than the minvalue?
              // If yes, set minvalue to newvalue and save point i,j values
              if (newValue <= minValue) {
                dx1 = i;
                dy1 = j;
                dx2 = k;
                dy2 = l;
                minValue = newValue;
              }
            }
          }
        }
      }
    }

    // Update lathub and longhub to the new point with a smaller value
    lathub1 += step * dx1;
    longhub1 += step * dy1;
    lathub2 += step * dx2;
    longhub2 += step * dy2;
    value = minValue;

  } while (value < oldValue);
  // Repeat all this while we can find a lower value than the previous one
  }

  // Test to see if the local minimum we just calculated is smaller than the current smallest known minimum
  // if yes: set globalValue to this newly found value, as well as the coordinates of this new minimum.
  if (value < globalValue){
    globalValue = value;
    globalY1 = lathub1;
    globalX1 = longhub1;
    globalY2 = lathub2;
    globalX2 = longhub2;
  }

  }
  // Average distance from hub to each location. Divide by 100 since 100 locations in GBplaces. Divide  
  // by 2 since globalValue is found by calculating the distance from the hub to the location and back.
  double avgDist = (globalValue*0.6214/100)/2;

  // Print the total number of function evaluations
  cout << "******************************" << endl;
  cout << "Function evaluations = " << f_evals << "\n";
  // Print the ideal location for each delivery hub 
  printf("- Delivery hub 1 should be placed at (%.4f, %.4f) \n", lathub1, longhub1);
  printf("- Delivery hub 2 should be placed at (%.4f, %.4f) \n", lathub2, longhub2);

  // Print out the average distance from a hub
  printf("- Placing the hubs at these coordinates leads to an average distance\nfrom any location to a \
hub of %.2f miles.\n", avgDist);

  return 0;
}
