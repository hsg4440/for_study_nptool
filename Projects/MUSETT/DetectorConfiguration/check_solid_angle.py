import math

def calculate_distance(point1, point2):
    return math.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2 + (point2[2] - point1[2])**2)

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

def compute_neighbor_distances(vertices):
    distances = []
    num_vertices = len(vertices)
    for i in range(num_vertices):
        v1 = vertices[i]
        v2 = vertices[(i + 1) % num_vertices]  # Next vertex (wraps around to the first vertex for the last one)
        dist = calculate_distance(v1, v2)
        distances.append(dist)
    return distances

def compute_center_distances(musett_data):
    distances = []
    for musett_points in musett_data:
        center_x = sum(point[0] for point in musett_points) / 4
        center_y = sum(point[1] for point in musett_points) / 4
        center_z = sum(point[2] for point in musett_points) / 4
        print(center_x,center_y,center_z)
        dist = calculate_distance((center_x, center_y, center_z), (0, 0, 0))
        distances.append(dist)
    return distances

if __name__ == "__main__":
    filename = "MUSETT_NPS.detector"
    musett_data = parse_detector_file(filename)

    # Compute distances between neighboring points for each MUSETT
    for i, musett_points in enumerate(musett_data, 1):
        distances = compute_neighbor_distances(musett_points)
        print(f"Distances between neighboring points for MUSETT {i}:")
        for j, dist in enumerate(distances, 1):
            print(f"  Distance {j}: {dist:.2f} mm")

    # Compute distances of MUSETT centers to the origin
    distances = compute_center_distances(musett_data)
    for i, dist in enumerate(distances, 0):
        print(f"Distance of MUSETT {i} center to origin: {dist:.2f} mm")
