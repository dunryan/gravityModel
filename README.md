# gravityModel
Spherical Harmonic Implementation of Gravitational Potential 

Compile ---> g++ *.cpp -o satSim.exe
Run     ---> ./satSim.exe
Plot    ---> python3 /.plotter_sat.py

Script to clean build, compile, run, plot
./run.sh

Switch from JGM3 to EGM2008 via top gravityModel option 
need to uncomment gravModel.h line 81 for JGM3 to work... ran out of steam

#### Results ####
Spent too much time trying to get the spherical harmonic grav potential partial derivs to work, it seems the process of computing the 
partial derivs is correct but my results seem off. Gravity in the ECI seems right initially but then falls off a cliff. This could be due to an issue somewhere in 
the numerical integration or how I implemented the partial deriv updates as position changes. 

With more time, I would dig deeper into the steps of the sim and understand where and when the gravity term exponentially decreases. I would also interrogate the 
coefficients used as well as the routine for partial derivative computation.

I also had an RK4 method that I used to verify the satellite dynamics initially but did not include that here. My plan was to get good results with Euler and then use RK4.

#### Shortcuts ####
Used the eigen library for the legendre polynomial generation. Ellicited ChatGPT's assistance with coding the partial derivative code using the equations from the NASA document https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/gravity-SphericalHarmonics.pdf.

