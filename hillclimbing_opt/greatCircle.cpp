// Final Assignment
// Read in the GBPlaces.csv and use an optimization algorithm over many
// iterations to find the ideal place to build a delivery hub


// includes and definitions

#include <iostream>
#include <cmath>

#define PI 3.14159

using namespace std;

// functions

double greatCircle(double lat1, double long1, double lat2, double long2) {
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

// main program

int main() {
  cout << "Manchester -> London distance is: " << greatCircle(53.48,-2.23743,51.51,-0.12574)/1609.3;
  return 0;
}