// Constants
const double G  = 6.67430e-11; // Gravitational constant [m^3/kg*s^2]
const double M  = 5.972e24;    // Earth mass [kg]
const double R  = 6378137.0;   // Earth's radius [m]
const double f  = 1 / 298.257223563; // Earth's flattening
const double e2 = f * (2 - f); // Square of eccentricity
double      gst = 1.7528311; // greenwich sidereal time [rad]

int GRAV_MODEL_ORDER = 2;

enum gravityModel {
	Simple,
	JGM3,
    EGM2008
    };

// struct to hold coefficients
struct Coefficient{

    int n;
    int m;
    double Cnm;
    double Snm;

};
