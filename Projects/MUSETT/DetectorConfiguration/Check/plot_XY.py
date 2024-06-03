import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 18})  # Increase this number for larger text globally
# Load the data from the file
def load_data(filename):
    detectors = {}
    cnt = -1
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 5:
                # Increment the detector count when 'X1_Y1' is found
                if "X1_Y1=" in parts[0]:
                    cnt += 1
                # Create a tuple key consisting of detector number and point identifier
                key = (f"DET{cnt}", parts[0][:-1])
                x, y, z = map(float, parts[1:4])
                if key in detectors:
                    detectors[key].append((x, y, z))
                else:
                    detectors[key] = [(x, y, z)]
    return detectors


def plot_detectors(detectors):

    dist = 8
    theta = 45*np.pi/180
    x_target = -dist*np.cos(theta)
    y_target = dist*np.sin(theta)
    #print(x_target,y_target)
    pos_target = np.array([x_target,y_target])
    fig, ax = plt.subplots(figsize=(9,9))

    # Iterate through all possible detector keys assuming the maximum number of detectors
    for i in range(4):  # Dynamically find the highest DET number
        start_key = (f"DET{i}", 'X1_Y128')
        end_key = (f"DET{i}", 'X128_Y128')

        # Check if both start and end points are present in the dictionary
        if start_key in detectors and end_key in detectors:
            start = np.array(detectors[start_key][0][:2])  # X, Y coordinates of the start point
            end = np.array(detectors[end_key][0][:2])    # X, Y coordinates of the end point
            center = 0.5*(start+end)
            ax.arrow(start[0], start[1], end[0] - start[0], end[1] - start[1], head_width=2, head_length=3, fc='black', ec='black', linewidth = 2)
            #if center[0]>0:
            ax.text(1.25*center[0],1.25*center[1], f'DET{i}', fontsize=18, ha='center', va='bottom')


            center_vector = center
            arrow_vector = end - start
            center_norm = np.linalg.norm(center_vector)
            arrow_norm = np.linalg.norm(arrow_vector)

            #if center_norm != 0 and arrow_norm != 0:
            angle_with_center = np.arccos(np.clip(np.dot(center_vector, arrow_vector) / (center_norm * arrow_norm), -1.0, 1.0))
            angle_with_center = np.degrees(angle_with_center)
            ideal_pos = np.array([center[0]/np.abs(center[0]),center[1]/np.abs(center[1])])
            ideal_pos[0]*=np.cos(np.pi/4)
            ideal_pos[1]*=np.sin(np.pi/4)
            ideal_pos*=center_norm
            distance_to_source = np.linalg.norm(center-pos_target)
            print(rf"DET{i} : distance to source {distance_to_source} mm")
            #print(rf"Angle between (0, center) and (X1_Y1, X128_Y1) for DET{i}: %.2f degrees"%angle_with_center)
            #print()
            difference = np.linalg.norm(ideal_pos-center_vector)
                #print(ideal_pos-center_vector)

            line_length = 50

            if (center[0]*center[1]>0):
                angle = -np.pi/4
            else :
                angle = np.pi/4
            dx = line_length * np.cos(angle)
            dy = line_length * np.sin(angle)
            dist = np.sqrt(center[0]**2+center[1]**2)
            #r'$\bar{\Omega}_{%d} = %.3f$' % (i, angle)
            ax.plot([x_target,center[0]],[y_target,center[1]], linestyle = "--", c = 'red')
            ax.plot([x_target,100*center[0]/np.abs(center[0])],[y_target,100*center[1]/np.abs(center[1])], linestyle = ":", c = 'orange', zorder = 20)
            ax.text(110*center[0]/np.abs(center[0]),110*center[1]/np.abs(center[1]), r"Center of distrib", fontsize=14, ha='center', va='bottom',c ='orange')
            ax.text(1.9*center[0],center[1],f"∆ = %.2f mm"%difference, c = 'darkorange',ha='center', va='bottom' )

            xoffset_text = 20
            if center[0]>0:
                if center[1]>0:
                    ax.text(xoffset_text+1.1*center[0]/2,0.8*center[1]/2, r"%.2f mm"%(distance_to_source), fontsize=14, ha='center', va='bottom',c ='red')
                else :
                    ax.text(xoffset_text+1.1*center[0]/2,1.2*center[1]/2, r"%.2f mm"%(distance_to_source), fontsize=14, ha='center', va='bottom',c ='red')
            else :
                if center[1]>0:
                    ax.text(-xoffset_text+1.1*center[0]/2,0.8*center[1]/2, r"%.2f mm"%(distance_to_source), fontsize=14, ha='center', va='bottom',c ='red')
                else :
                    ax.text(-xoffset_text+1.1*center[0]/2,1.2*center[1]/2, r"%.2f mm"%(distance_to_source), fontsize=14, ha='center', va='bottom',c ='red')
    # Plot the origin and axes
    ax.arrow(0, -125, 0, 50, head_width=5, head_length=5, fc='gray', ec='gray', linewidth = 5)
    ax.text(0,-140, "BEAM", c = 'gray', ha='center', va='bottom')
    ax.scatter(x_target, y_target, color='red', s=100)  # big point at the center
    ax.axhline(0, color='black', linewidth=0.5)  # X axis
    ax.axvline(0, color='black', linewidth=0.5)  # Y axis

    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.xlabel('X(geometre) [mm]')
    plt.ylabel('Y(geometre) [mm]')
    plt.title('Detector Orientation, Y128')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    #plt.savefig("PlotGeometry_Y128_shift_target.pdf")


def plot_tilted_target(detectors):
    fig, ax = plt.subplots(figsize=(9,9))
    colors = ['black','red','blue','orange']
    # Iterate through all possible detector keys assuming the maximum number of detectors
    for i in range(4):  # Dynamically find the highest DET number
        start_key = (f"DET{i}", 'X1_Y1')
        end_key = (f"DET{i}", 'X128_Y1')
        col = colors[i]

        # Check if both start and end points are present in the dictionary
        if start_key in detectors and end_key in detectors:
            start = np.array(detectors[start_key][0][:2])  # X, Y coordinates of the start point
            end = np.array(detectors[end_key][0][:2])    # X, Y coordinates of the end point
            center = 0.5*(start+end)
            ax.arrow(start[0], start[1], end[0] - start[0], end[1] - start[1], head_width=2, head_length=3, fc=col, ec=col, linewidth = 2)
            #if center[0]>0:
            ax.text(1.25*center[0],1.25*center[1], f'DET{i}', fontsize=18, ha='center', va='bottom', c = col)
            plt.plot([0,start[0]],[0,start[1]],":", c = col, linewidth = 2, zorder = 50)
            plt.plot([0,end[0]],[0,end[1]],":", c = col, linewidth = 2, zorder = 50)
            ax.text(1.2*end[0],1.2*end[1], 'X128', fontsize=14, ha='center', va='bottom', c = col)
            ax.text(1.2*start[0],1.2*start[1], 'X1', fontsize=14, ha='center', va='bottom', c = col)


    # Plot the origin and axes
    ax.arrow(0, -125, 0, 50, head_width=5, head_length=5, fc='gray', ec='gray', linewidth = 5)
    ax.text(0,-140, "BEAM", c = 'gray', ha='center', va='bottom')
    #ax.text(1.25*center[0],1.25*center[1], f'DET{i}', fontsize=18,)

    ax.axhline(0, color='black', linewidth=0.5)  # X axis
    ax.axvline(0, color='black', linewidth=0.5)  # Y axis
    angle_target = 10
    size_target = 34
    angle_target*= np.pi/180
    x = size_target * np.cos(angle_target) /2
    y = size_target * np.sin(angle_target) /2
    ax.plot([-x,x],[-y,y], c= 'green', linewidth = 5, zorder =2 )
    #ax.plot([-x,x],[y,-y], c= 'gray')
    #ax.scatter(0, 0, color='red', s=15, marker = '+', linewidth=30, alpha = 0.7)  # big point at the center

    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.xlabel('X(geometre) [mm]')
    plt.ylabel('Y(geometre) [mm]')
    plt.title(r'Effect of Target. stripY = 1,$\theta$(target) = 10°')
    plt.grid(True)
    plt.tight_layout()
    #plt.savefig("PlotTarget_Y1_shift_target.pdf")
    plt.show()

# Replace 'MUSETT_NPS.detector' with your actual file path
detectors = load_data('../MUSETT_NPS.detector')
#for key,value in detectors.items() :
    #print(f"key = {key} -> value = {value}")
plot_detectors(detectors)
#plot_tilted_target(detectors)
