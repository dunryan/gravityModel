import pandas as pd
import matplotlib.pyplot as plt

# import math as math

# Read the data from CSV
df = pd.read_csv('satellite_orbit_3d.csv')

# R = 6.371e6     # Earth Raduis

print(df.columns)

# Extracting columns from dataframe
time = df['Time ']

x_eci = df['X ']
y_eci = df['Y ']
z_eci = df['Z ']
vx_eci = df['dX ']
vy_eci = df['dY ']
vz_eci = df['dZ ']
ax_eci = df['aX ']
ay_eci = df['aY ']
az_eci = df['aZ ']


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x_eci, y_eci, z_eci)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()


# Satellite Radius vs Orbit Angle
plt.figure(figsize=(10, 6))  # Adjust the figure size if necessary
plt.plot(time, ax_eci, marker='o', linestyle='-', color='b', label='x acc [m/s^2]')
# plt.plot(time, ax_eci, marker='o', linestyle='-', color='r', label='y acc [m/s^2]')
# plt.plot(time, ax_eci, marker='o', linestyle='-', color='k', label='z acc [m/s^2]')
plt.title('Acceleration of Satellite Orbiting Earth')
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s^2]')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()