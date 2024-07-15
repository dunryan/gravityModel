#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "../eigen/Eigen/Dense"
#include "./include/gravModel.h"

using namespace std;

// Select which gravity model to use 
// 1 = simple gravity model
// 2 = JGM-3 gravitational potential model
// 3 = EGM2008 gravitational undulation and field model 
gravityModel gravModelRun = EGM2008;

// Function to return gravitational force
// Input: 3d position vector 
// Output: 3d gravitational force vector
Eigen::Vector3d gravitationalForce(const Eigen::Vector3d &position) {
    double r = position.norm();
    return -G * M / (r * r * r) * position;
}

// Function to perform Euler integration step in 3D
void eulerStep(Eigen::Vector3d &position, Eigen::Vector3d &velocity, std::vector<double> &g, double dt) {
    // Update position components
    position += velocity * dt;

    // Update velocity components
    Eigen::Vector3d gEig = Eigen::Vector3d::Zero();

    gEig[0] = g[0]; // Assign x component
    gEig[1] = g[1]; // Assign y component
    gEig[2] = g[2]; // Assign z component

    velocity += gEig * dt;
}

// Main sim loop
int main() {

    if (gravModelRun == JGM3){
        provideHarmonicCoefficientsJGM3(C,S);
    }
    else if (gravModelRun == EGM2008){
        provideHarmonicCoefficientsEGM2008(C,S);
    }

    // Computation of ECI Initial Position Vector
    double altitude = 400e3; // 400 km above Earth's surface
    double initialRadius = R + altitude;

    // Inclination angle 
    double inclination = 28.5 * M_PI / 180.0;

    // Initial position
    Eigen::Vector3d position(initialRadius * cos(inclination), 0.0, initialRadius * sin(inclination));

    // Orbital velocity for a circular orbit: v = sqrt(G * M / r)
    double orbitalVelocity = sqrt(G * M / initialRadius);

    // Initial velocity
    Eigen::Vector3d velocity(0.0, orbitalVelocity * cos(inclination), orbitalVelocity * sin(inclination));

    // Simulation parameters
    double dt = 10.0;        // Time step (s)
    double simulationTime = 864000.0;  // Total simulation time [s] == 1 day

    cout << "Initial position: " << position[0] << " " << position[1] << " " << position[2] << endl;
    cout << "Initial velocity: " << velocity[0] << " " << velocity[1] << " " << velocity[2] << endl;
    cout << "Simulation time: " << simulationTime << endl;
    cout << "Simulation integration time: " << dt << endl;

    double r = position.norm(); // distance from Earth's center [m]
    double lambda = atan2(position[1], position[0]); // longitude 
    double theta = acos(position[2] / r); // latitude

    // Sanity check that radial component of gravity is -9.81... passed
    // double r = R + 1000.0;
    // double lambda = 0.0; 
    // double theta =  0.0;

    // Compute the gravitational acceleration in spherical coordinates using radius, latitude, longitudinal 
    std::vector<double> gSpherical = computeGravitationalAcceleration(r, theta, lambda, C, S);

    // Convert acceleration from spherical coordinates to cartesian coordinates 
    // Source: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/gravity-SphericalHarmonics.pdf
    gEcef[0] = cos(theta)*cos(lambda)*gSpherical[0] - sin(theta)*sin(lambda)*gSpherical[1] - sin(lambda)*gSpherical[2];
    gEcef[1] = cos(theta)*sin(lambda)*gSpherical[0] - sin(theta)*sin(lambda)*gSpherical[1] - cos(lambda)*gSpherical[2];
    gEcef[2] = cos(theta)*gSpherical[0] - sin(theta)*gSpherical[1];

    // Convert gravity in ECEF to gravity in ECI to pass into numerical integration
    gEci[0] = cos(gst)*gEcef[0] - sin(gst)*gEcef[1];
    gEci[1] = sin(gst)*gEcef[0] + cos(gst)*gEcef[1];
    gEci[2] = gEcef[2];

    // Open a file to store the simulation data
    ofstream outputFile;
    outputFile.open("satellite_orbit_3d.csv");
    outputFile << "Time ,X ,Y ,Z ,dX ,dY ,dZ ,aX ,aY ,aZ \n";

    // Simulation loop
    for (double t = 0.0; t <= simulationTime; t += dt) {

        double r = position.norm(); // distance from Earth's center [m]
        double lambda = atan2(position[1], position[0]); // longitude 
        double theta = acos(position[2] / r); // latitude

        // Compute the gravitational acceleration in spherical coordinates using radius, latitude, longitudinal 
        std::vector<double> gSpherical = computeGravitationalAcceleration(r, theta, lambda, C, S);

        // Convert acceleration from spherical coordinates to cartesian coordinates 
        // Source: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/gravity-SphericalHarmonics.pdf
        gEcef[0] = cos(theta)*cos(lambda)*gSpherical[0] - sin(theta)*sin(lambda)*gSpherical[1] - sin(lambda)*gSpherical[2];
        gEcef[1] = cos(theta)*sin(lambda)*gSpherical[0] - sin(theta)*sin(lambda)*gSpherical[1] - cos(lambda)*gSpherical[2];
        gEcef[2] = cos(theta)*gSpherical[0] - sin(theta)*gSpherical[1];

        // Convert gravity in ECEF to gravity in ECI to pass into numerical integration - same source 
        gEci[0] = cos(gst)*gEcef[0] - sin(gst)*gEcef[1];
        gEci[1] = sin(gst)*gEcef[0] + cos(gst)*gEcef[1];
        gEci[2] = gEcef[2];

        // Output current state to file
        outputFile << t << "," << position.x() << "," << position.y() << "," << position.z() << "," <<  velocity.x() << "," << velocity.y() << "," << velocity.z() << "," << gEci[0] << "," << gEci[1] << "," << gEci[2] << "\n";

        // Perform Euler integration step
        eulerStep(position, velocity, gEci, dt);
    }

    // Close the output file
    outputFile.close();

    cout << "Simulation complete. Output saved to satellite_orbit_3d.csv\n";

    return 0;
}