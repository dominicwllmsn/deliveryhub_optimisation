// A program to read in GBplaces.csv

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

int main() {

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
            if (row == -1) {
                row++;
            } else {
                // define input stringstream object iss
                // the use of the input stringstream creates a buffer (like cin or dataFile) 
                // which can be used in the getline() function (to split up a string)
                istringstream iss;
                iss.str(line); // assign to iss the line of the .csv we want to split up
                string field; // create field to store each substring from getline
                int column = 1;

                // while there are substrings we can split up about the ',' delimiter char
                while (getline(iss, field, ','))
                {   
                    if (column < 3) {
                        column++;
                    } else {
                        // define element as double of the field str and append to temp 
                        double element = atof(field.c_str());
                        if (column == 3) {
                            population[row] = element;
                        } else if (column == 4) {
                            latitude[row] = element;
                        } else {
                            longitude[row] = element;
                        }
                        column++;
                    }
                }
                printf("(%5.2f,%5.2f) -> %7.0f \n", latitude[row],longitude[row],population[row]);
                row++;
            }
        }    

    // now at the file end, close it
    dataFile.close();
    cout << "Input file read and closed" << endl;
  } else {
    // if the file couldn't be opened, write a message and end main()
    cout << "Unable to open input file" << endl;
    return 1;
  }
  
  return 0;
}