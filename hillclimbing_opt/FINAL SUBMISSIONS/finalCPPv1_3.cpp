/*
  Final Assignment 
  Dominic Williamson - ID: 9846595
  version 1.3 : TWO HUBS + STEP OPTIMIZATION + NEAREST CITY 
  Program reads in the GBPlaces.csv and uses an optimization algorithm over many
  iterations to find the ideal place to build two delivery hubs.
  The step length for this algorithm is optimized by decreasing it by a factor of 
  ten over several iterations, allowing for more efficient and more accurate
  determination of the global minima. This version utilises the Great Circle distance
  formula which takes into account the curvature of the Earth's surface/that the planet 
  is approximately spherical. 
  The program also travels from its current location to its nearest neighbour (where
  feasible), creating a more accurate delivery model.
*/


// includes and definitions
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>
#include <vector>

#define PI 3.14159

int f_evals = 0; // counter for the number of function evaluations

using namespace std;

double greatCircle(double lat1, double long1, double lat2, double long2) {
        /*
        Inputs: latitudes and longitudes of two locations on Earth
        Output: the Great Circle distance between the two locations in kilometers
        */
        f_evals++;

        double dLat = lat2 - lat1;
        double dLong = long2 - long1;
        double R = 6371;
        double a = pow(sin(dLat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dLong/2),2);
        double c = 2 * atan2((pow(a,0.5)), pow((1-a),0.5));

        return R * c;
}

double fitness(double lat1, double long1, double lat2, double long2, vector< vector<double> > data) {
  /*
    Inputs: latitudes and longitudes of two locations on Earth, a 2D vector of location coordinates
    Output: the (minimum) sum of the distances from the hubs to each location using the Vincenty 
            formulae
  */

  // Declare necessary variables to use the Great Circle formula 
  double dist1, dist2;
  // Declare necessary variables for the "nearest city" extension
  vector<int> blacklist; // Vector to hold blacklisted towns/cities
  double nearLat; 
  double nearLong;
  int nearIndex;
  double distance;
  double minDist = 1e10; // Current minimum distance from the current location to the test 
                         // - initially set to a very large value

  // Start with f, the sum variable, being equal to 0
  double f = 0;

  lat1 *= PI/180.;
  long1 *= PI/180.;
  lat2 *= PI/180.;
  long2 *= PI/180.;

  // Loop over all the latitudes and longitudes in the matrix
  for (int j = 0; j < data.size(); j++) {

    // if 'j' is in the blacklist, move onto the next one
    if (find(blacklist.begin(), blacklist.end(), j) != blacklist.end())
      continue;

    // reset the necessary variables for each value of j
    double nearLat; 
    double nearLong;
    int nearIndex;
    double distance;
    double minDist = 1e10;

    // iterate over all locations which are NOT in the blacklist and ARE 
    // within 0.5 degrees of 'j' in both latitude and longitude
    for (int k = 0; k < data.size(); k++) {

      if (k == j || find(blacklist.begin(), blacklist.end(), k) != blacklist.end())
        continue;
      if (abs(data[k][0] - data[j][0]) > 0.5 || abs(data[k][1] - data[j][1]) > 0.5)
        continue;

      // find the distance from location 'j' to 'k'
      distance  = greatCircle(data[k][0]*PI/180, data[k][1]*PI/180, data[j][0]*PI/180, data[j][1]*PI/180);

      // if less than minDist, replace and go to the next 'k'. If not, simply go to the next 'k'
      if (distance < minDist) {
        minDist = distance;
        nearLat = data[k][0];
        nearLong = data[k][1];
        nearIndex = k;
      }
    }
    // blacklist j so that the algorithm won't visit it again (since a truck has already been there)
    blacklist.push_back(j);



    // Convert the latitudes and longitudes to radians
    data[j][0] = (PI*data[j][0]) / 180;
    data[j][1] = (PI*data[j][1]) / 180;

    dist1 = greatCircle(lat1, long1, data[j][0], data[j][1]);
    dist2 = greatCircle(lat2, long2, data[j][0], data[j][1]);

    if (minDist < 80) {
      // if minDist < 80 km (50 miles), go to that location and then go back to the hub
      // if not, only travel to j and back to the hub
      // (this conditional means we avoid adding minDist if it is still equal to 1e10)
    dist1 += minDist + greatCircle(lat1, long1, nearLat*PI/180, nearLong*PI/180);
    dist2 += minDist + greatCircle(lat2, long2, nearLat*PI/180, nearLong*PI/180);
    blacklist.push_back(nearIndex);
    } else {
      dist1 *= 2;
      dist2 *= 2;
    }

    if (dist1 <= dist2) {
      f += dist1;
    } else {
      f += dist2;
    }
  }

  f_evals++;
  return f;

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
  int N = 5; // number of iterations of the for-loop (i.e. number of random coordinates tested)  
  double globalY1, globalX1, globalY2, globalX2; // current global minimum of each lat/long value
  double globalValue = 1e10; // current global minimum of distance - set initially to a v.large number

  srand(time(NULL)); // seeds the random number generator

  for (int n = 0; n < N; n++)
  {
  // pick the first hub's starting point at a random point between the min & max long/lat values
  lathub1 = random_number(minLat, maxLat, 100);
  longhub1 = random_number(minLong, maxLong, 100);

  // pick the second hub's starting point at a random point between the min & max long/lat values
  lathub2 = random_number(minLat, maxLat, 100);
  longhub2 = random_number(minLong, maxLong, 100);

  // Using the fitness function we can calculate the total distance travelled
  value = fitness(lathub1, longhub1, lathub2, longhub2, coords);

  double step = 0.1; // step size to move in x and y

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
  double avgDist = (globalValue/1.609/100);

  cout << "******************************" << endl;
  // Print the total number of function evaluations
  cout << "Function evaluations = " << f_evals << "\n";

  // Print the ideal location for each delivery hub 
  printf("- Delivery hub 1 should be placed at (%.4f, %.4f) \n", lathub1, longhub1);
  printf("- Delivery hub 2 should be placed at (%.4f, %.4f) \n", lathub2, longhub2);

  // Print out the total distance travelled with our best solution
  printf("- Placing the hubs at these coordinates leads to an average distance\nfrom any location to a \
hub of %.2f miles.\n", avgDist);

  return 0;
}