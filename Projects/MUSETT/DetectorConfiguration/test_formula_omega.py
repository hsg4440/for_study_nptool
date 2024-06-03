import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})

def solid_angle_square_integration(z, width, length):
    # Function to compute the differential solid angle
    def d_omega(x, y):
        # z = 0  # Since we're integrating over x and y, we assume z = 0 (point of observation)
        return z / np.sqrt(x**2 + y**2 + z**2)**3

    # Perform the double integration
    solid_angle, _ = dblquad(d_omega, -width/2, width/2, lambda x: -length/2, lambda x: length/2)

    return solid_angle

# Example points of a square
width = 97.220
length = 97.220

# Compute the solid angle using integration over x and y
distances = np.array([78.95, 79.18, 78.41, 83.06])
distances = np.array([101.12])

solid_angles = []
total = 0
for i, z in enumerate(distances):
    angle_integration = solid_angle_square_integration(z, width, length)
    solid_angles.append(angle_integration / (4 * np.pi))
    print(f"Det {i}, Solid Angle / 4 π = {angle_integration / (4 * np.pi)}")
    total += angle_integration



print('total solid angle / 4π = ',total/(4*np.pi))
'''
# Normalize the values
max_solid_angle = max(solid_angles)
normalized_solid_angles = np.array([angle / max_solid_angle for angle in solid_angles])


solid_angles = np.array(solid_angles)
estimate = solid_angles * 10000000

#estimate *= 897158 / estimate[0]
#print(estimate)



#G4_output = np.array([9150, 8880, 9026, 8354])
#G4_output = np.array([9015, 8841, 9095, 8240])
G4_output = np.array([897158,890469,905700,827856])
G4_output_errors = np.sqrt(G4_output)
#print(G4_output_errors)
#G4_output_errors *=
# Plotting histogram with annotations
#plt.bar(range(len(normalized_solid_angles)), normalized_solid_angles, label='My Calculation')
plt.bar(range(len(G4_output)),G4_output, label = 'G4 sim')
plt.errorbar(range(len(G4_output)), G4_output, yerr=G4_output_errors, fmt='ro-', c = 'black', label = 'G4 +/ sigma')  # Plot the values from another calculation with error bars
plt.plot([0,1,2,3],estimate,'s', c = 'red', label = 'geom estimate')
plt.legend()
#plt.errorbar(range(len(G4_output)), G4_output, yerr=G4_output_errors, fmt='ro-', label='G4 Output')  # Plot the values from another calculation with error bars
plt.xlabel('Detector Number')
plt.ylabel(r'N')
plt.title('Geometric solid angles')
#plt.ylim(8200,9300)

# Annotate each bar with the corresponding solid angle value
#for i, angle in enumerate(normalized_solid_angles):
#    plt.text(i, solid_angles[i], r'$\bar{\Omega}_{%d} = %.3f$' % (i, angle), ha='center', va='bottom', fontsize=16, weight='bold', color = 'red')

#omega_ratio_2 = solid_angles[2] / solid_angles[2]  # Assuming Omega_2 != 0
#plt.text(3.1, solid_angles[2], r'$\bar{\Omega} = \Omega / \Omega_2$', ha='center', va='bottom', fontsize=18, weight='bold', color='red', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

plt.xticks([0,1,2,3],[0,1,2,3])

#plt.ylim((0.08, 0.0950))
plt.show()
'''
