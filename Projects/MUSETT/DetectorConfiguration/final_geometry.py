import math
from scipy.integrate import dblquad
import numpy as np

def parse_detector_file(filename):
    with open(filename, 'r') as file:
        data = file.readlines()

    musett_data = []
    musett_points = []
    for line in data:
        if not line.strip():  # Skip empty lines
            continue
        parts = line.split("=")
        if parts[0].strip().startswith("X1_Y") or parts[0].strip().startswith("X128_Y"):  # Check if the line contains point coordinates
            x, y, z = map(float, parts[1].split()[:3])
            musett_points.append((x, y, z))
        elif line.startswith("MUSETT"):  # If line doesn't contain point coordinates, it's the start of a new MUSETT entry
            if musett_points:  # If there are points in musett_points, add them to musett_data
                musett_data.append(musett_points)
                musett_points = []  # Reset musett_points for the next entry

    # Add the last set of points
    if musett_points:
        musett_data.append(musett_points)

    return musett_data


def find_center(corner_points):
    # Extracting x, y, and z coordinates of each corner point
    x_coords = [point[0] for point in corner_points]
    y_coords = [point[1] for point in corner_points]
    z_coords = [point[2] for point in corner_points]

    # Computing the center point
    center_x = sum(x_coords) / len(corner_points)
    center_y = sum(y_coords) / len(corner_points)
    center_z = sum(z_coords) / len(corner_points)

    return (center_x, center_y, center_z)

def compute_distance_to_origin(center):
    # Computing the distance from the center to the origin (0, 0, 0)
    distance = math.sqrt(center[0]**2 + center[1]**2 + center[2]**2)
    return distance

def compute_solid_angle(corner_points):

    length = 97.220
    # Function to compute the differential solid angle
    def d_omega(x, y):
        # Extracting the z coordinate of the first corner point
        z = corner_points[0][2]
        z-= 0.150
        return z / math.sqrt(x**2 + y**2 + z**2)**3


    # Getting the minimum and maximum values of x and y coordinates
    min_x, max_x = -length/2,length/2
    min_y, max_y = -length/2,length/2

    # Compute the solid angle using double integration
    solid_angle, _ = dblquad(d_omega, min_x, max_x, lambda x: min_y, lambda x: max_y)

    return solid_angle

# Example usage:

#corner_points = [(18.16, 92.38, -39.51), (18.05, 92.32, 60.09), (88.66, 22.07, 60.12), (88.77, 22.13, -39.47)]
musett_data = parse_detector_file('MUSETT_TEST.detector')
corner_points = musett_data[0]
center = find_center(corner_points)
distance_to_origin = compute_distance_to_origin(center)
solid_angle = compute_solid_angle(corner_points)
print("Center of the detector:", center)
print("Distance from center to origin:", distance_to_origin)
print("Geometric efficiency = Ω/4π = ",solid_angle/(4*np.pi))
