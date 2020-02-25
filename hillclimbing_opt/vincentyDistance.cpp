// C++ Standard Libraries
#include <cmath>
#include <iostream>
#include <ctime>

#define PI 3.14159

/**
 * Vincenty Distance
*/
double Vincenty_Distance( const double& latitude_01, const double& longitude_01,
                          const double& latitude_02, const double& longitude_02,
                          const double& a,
                          const double& b )
{
    // Flattening
    const double f = (a-b)/a;

    // tan U1
    const double tan_U1 = (1-f) * std::tan(latitude_01);
    const double tan_U2 = (1-f) * std::tan(latitude_02);

    // Longitudinal Distance
    const double cos_U1 = 1 / std::sqrt(1 + tan_U1 * tan_U1);
    const double cos_U2 = 1 / std::sqrt(1 + tan_U2 * tan_U2);
    const double sin_U1 = tan_U1 * cos_U1;
    const double sin_U2 = tan_U2 * cos_U2;

    // Iterate until complete
    const double L = longitude_02 - longitude_01;
    double lambda = L;
    double diff, sigma;
    double cos_alpha_sq, cos_2sigma_m;
    double cos_sigma, sin_sigma;

    while( true ){

        // 
        double sin_lambda = std::sin( lambda );
        double cos_lambda = std::cos( lambda );

        double c1 = (cos_U2 * sin_lambda)*(cos_U2 * sin_lambda);
        double c2 = (cos_U1 * sin_U2);
        double c3 = (sin_U1 * cos_U2 * cos_lambda);


        //  sin sigma
        sin_sigma = std::sqrt( c1 + ( c2 - c3 )*( c2 - c3 ) );

        // cos sigma
        cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lambda;

        // sigma
        sigma = std::atan2( sin_sigma, cos_sigma );

        // sin alpha
        double sin_alpha = (cos_U1 * cos_U2 * sin_lambda)/(sin_sigma);

        // cos^2 alpha
        cos_alpha_sq = 1 - (sin_alpha*sin_alpha);

        // cos^2 2sigmam
        cos_2sigma_m = cos_sigma - (2 * sin_U1 * sin_U2)/(cos_alpha_sq);

        // C
        double C = (f/16.0) * cos_alpha_sq * (4 + f * (4 - 3 * cos_alpha_sq));

        // Update Lambda
        diff = lambda;
        lambda = L + (1-C) * f * sin_alpha * (sigma + C * sin_sigma * ( cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m*cos_2sigma_m)));
        diff = lambda - diff;
        if( std::fabs(diff) < 0.00001 ){ break; }
    }

    // U2
    double u_sq = cos_alpha_sq  * (a*a - b*b)/(b*b);

    // Compute A, B
    double A = 1 + (u_sq/16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));

    double B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))));

    // Sigma
    double cos_2sigma_m_sq = cos_2sigma_m * cos_2sigma_m;
    double delta_sigma = B * sin_sigma * ( cos_2sigma_m + (B/4.0) * (cos_sigma * (-1 * 2 * cos_2sigma_m_sq ) - (B/6.0) * cos_2sigma_m * (-3 + 4 * sin_sigma*sin_sigma) * (-3 + 4 * cos_2sigma_m_sq)));

    // Distance
    double s = b * A * (sigma - delta_sigma);
    return s;
}

double greatCircle(double lat1, double long1, double lat2, double long2) {
        /*
        Inputs: latitudes and longitudes of two locations on Earth
        Output: 'Great Circle' distance between the two locations in meters
        */
        lat1 = lat1;
        long1 = long1;
        lat2 = lat2;
        long2 = long2;
        double dLat = lat2 - lat1;
        double dLong = long2 - long1;
        double R = 6371e3;
        double a = pow(sin(dLat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dLong/2),2);
        double c = 2 * atan2((pow(a,0.5)), pow((1-a),0.5));

        return R * c;
}


/**
 * @brief Main Function
 */
int main( int argc, char* argv[] )
{

    // Set our coordinates
    const double latitude_01  =   -4.25763 * M_PI / 180.0;
    const double longitude_01 = 50.3715 * M_PI / 180.0;

    const double latitude_02 =    1.29834 * M_PI / 180.0;
    const double longitude_02 = 57.1437 * M_PI / 180.0;

    // Set the datum components
    const double a = 6378137.0;
    const double b = 6356752.314245;

    // Vincenty Distance
    clock_t begin = clock();
    double distance = Vincenty_Distance( latitude_01, longitude_01,
                                         latitude_02, longitude_02,
                                         a, b);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin);

    std::cout << "Elapsed: " << (elapsed_secs) << std::endl;

    begin = clock();
    double distanceGC = greatCircle(latitude_01, longitude_01,
                                    latitude_02, longitude_02);
    end = clock();
    elapsed_secs = double(end - begin);
    std::cout << "Elapsed: " << (elapsed_secs) << std::endl;

    // Distance
    std::cout << std::fixed << "Distance: " << distance << std::endl;
    std::cout << std::fixed << "DistanceGC: " << distanceGC << std::endl;


    // Exit
    return 0;
}