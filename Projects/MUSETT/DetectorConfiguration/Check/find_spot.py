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

def qualitative_hit_pattern(detectors, Xsource, Ysource):
    dist_strips = np.zeros((4,128))
    ignore_strips = [[10,26,38],[],[104],[22]]
    #ignore_strips = [[],[],[],[]]
    for i in range(4):  # Dynamically find the highest DET number
        start_key = (f"DET{i}", 'X1_Y128')
        end_key = (f"DET{i}", 'X128_Y128')
        pos_strips = []
        min_dist = 100000
        min_dist_strip = 0
        mean_strip = 0
        all_dist_square = 0
        # Check if both start and end points are present in the dictionary
        if start_key in detectors and end_key in detectors:
            start = np.array(detectors[start_key][0][:2])  # X, Y coordinates of the start point
            end = np.array(detectors[end_key][0][:2])    # X, Y coordinates of the end point
            center = 0.5*(start+end)
            for j in range(128):
                if j not in ignore_strips[i] :
                    xi = start[0] + j*(end[0]-start[0])/127
                    yi = start[1] + j*(end[1]-start[1])/127
                    dist_strips[i][j] = np.sqrt((xi-Xsource)**2 + (yi-Ysource)**2)
                    mean_strip += j / dist_strips[i][j]**2
                    all_dist_square += 1/dist_strips[i][j]**2
                    if(dist_strips[i][j] < min_dist):
                        min_dist = dist_strips[i][j]
                        min_dist_strip = j
            mean_strip *= 1/all_dist_square
            print("DET", i, ", max strip =", min_dist_strip, ", mean = ", mean_strip )

detectors = load_data('../MUSETT_NPS.detector')

def computeMean(Xsource,Ysource):
    ignore_strips = [[10,26,38],[],[104],[22]]
    l_meanstrip = np.zeros(4)
    l_maxstrip = np.zeros(4)
    for i in range(4):  # Dynamically find the highest DET number
        start_key = (f"DET{i}", 'X1_Y128')
        end_key = (f"DET{i}", 'X128_Y128')
        mean_strip = 0
        all_dist_square = 0
        min_dist = 100000
        min_dist_stirp = -2
        # Check if both start and end points are present in the dictionary
        if start_key in detectors and end_key in detectors:
            start = np.array(detectors[start_key][0][:2])  # X, Y coordinates of the start point
            end = np.array(detectors[end_key][0][:2])    # X, Y coordinates of the end point
            center = 0.5*(start+end)
            for j in range(128):
                if j not in ignore_strips[i] :
                    xj = start[0] + j*(end[0]-start[0])/127
                    yj = start[1] + j*(end[1]-start[1])/127
                    dist_j = np.sqrt((xj-Xsource)**2 + (yj-Ysource)**2)
                    mean_strip += j / dist_j**2
                    all_dist_square += 1/dist_j**2
                    if(dist_j < min_dist):
                        min_dist = dist_j
                        min_dist_strip = j
            mean_strip *= 1/all_dist_square
            l_meanstrip[i] = mean_strip
            l_maxstrip[i] = min_dist_strip
    return l_maxstrip

def plot2D():
    goal = np.array([61.14,64.18,61.84,63.91])

    x_values = np.linspace(-10, 10, 100)
    y_values = np.linspace(-10, 10, 100)
    refValue = goal[2]

    X, Y = np.meshgrid(x_values, y_values)

    # Initialize the result array
    result0 = np.zeros_like(X)
    result1 = np.zeros_like(X)
    result2 = np.zeros_like(X)
    result3 = np.zeros_like(X)
    l_resul = [result0,result1,result2,result3]


    x_size_target = 34
    y_size_target = 1
    for i in range(len(x_values)):
        for j in range(len(y_values)):
            res = computeMean(X[j, i], Y[j, i])
            result0[j, i] = res[0] - goal[0]  # Adjust refValue as needed
            result1[j, i] = res[1] - goal[1]  # Adjust refValue as needed
            result2[j, i] = res[2] - goal[2]  # Adjust refValue as needed
            result3[j, i] = res[3] - goal[3]  # Adjust refValue as needed


    # Plotting the results
    levels = np.arange(-20, 20.1, 1)
    for k in range(4):
        plt.figure(figsize=(10, 8))
        cp = plt.contourf(X, Y, l_resul[k], cmap='coolwarm',vmin=-20, vmax=20,levels=levels)
        #cp = plt.contourf(X, Y, l_resul[k], cmap='coolwarm')
        #ticks = np.arange(-4, 4.5, 0.5)  # Includes 4.0, steps by 0.5

        # Create colorbar with specified ticks
        colorbar = plt.colorbar(cp)
        colorbar.set_label('∆ strip')
        #plt.colorbar(cp)  # Add a colorbar to a plot
        #plt.plot(min_x, min_y, 'r*', markersize=15, label='Minimum Value')
        plt.title(f'Max (geo - Data), DET{k}')
        #plt.legend()
        plt.xlabel('X(source) [mm]')
        plt.ylabel('Y(source) [mm]')
        #plt.savefig(f"delta_max_Det{k}.pdf")
        plt.show()



def plot_config(posTarget,angleTarget, x0Beam):
    x0 = posTarget[0]
    y0 = posTarget[1]
    angleTarget*=-np.pi/180
    y0Beam = y0 - (x0-x0Beam)*np.tan(angleTarget)
    print('y0 Beam =', y0Beam )
    X_dets = 0
    Y_dets = 0
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
            start[0]+=X_dets
            start[1]+=Y_dets

            end[0]+=X_dets
            end[1]+=Y_dets
            center = 0.5*(start+end)
            ax.arrow(start[0], start[1], end[0] - start[0], end[1] - start[1], head_width=2, head_length=3, fc=col, ec=col, linewidth = 2)
            #if center[0]>0:
            ax.text(1.25*center[0],1.25*center[1], f'DET{i}', fontsize=18, ha='center', va='bottom', c = col)
            plt.plot([x0Beam,start[0]],[y0Beam,start[1]],":", c = col, linewidth = 2, zorder = 50)
            plt.plot([x0Beam,end[0]],[y0Beam,end[1]],":", c = col, linewidth = 2, zorder = 50)
            ax.text(1.2*end[0],1.2*end[1], 'X128', fontsize=14, ha='center', va='bottom', c = col)
            ax.text(1.2*start[0],1.2*start[1], 'X1', fontsize=14, ha='center', va='bottom', c = col)
            if(i==0):
                print('X128_Y1 = ', end, ',angle =', np.arctan((end[1]-y0)/(end[0]-x0))*180/np.pi)


    # Plot the origin and axes
    #ax.text(1.25*center[0],1.25*center[1], f'DET{i}', fontsize=18,)

    ax.axhline(0, color='black', linewidth=0.5)  # X axis
    ax.axvline(0, color='black', linewidth=0.5)  # Y axis


    size_target = 12
    angle_target = angleTarget
    #angle_target*= np.pi/180
    #for theta in np.linspace(-0.05*np.pi,0.05*np.pi,3):
    x = size_target * np.cos(angle_target) /2
    y = size_target * np.sin(angle_target) /2
    ax.plot([-x+x0,x+x0],[y0-y,y0+y], c= 'green', linewidth = 2, zorder =2 )


    size_holder = 34
    x_holder = size_holder * np.cos(angle_target) /2
    y_holder = size_holder * np.sin(angle_target) /2
    dx = x_holder-x
    dy = y_holder -y
    ax.plot([-x_holder+x0,-x_holder+x0+dx],[y0-y_holder,y0-y_holder+dy], '-',c= 'dimgray', linewidth = 2.5, zorder =20 )
    ax.plot([x_holder+x0,x_holder+x0-dx],[y0+y_holder,y0+y_holder-dy], '-',c= 'dimgray', linewidth = 2.5, zorder =20 )
    #ax.plot([x_holder+x0,x_holder+x0-dx],[y0-y,y0+y], '-',c= 'darkgray', linewidth = 2.5, zorder =20 )
    #ax.plot([-x+x0,x+x0],[y0-y,y0+y], '-',c= 'black', linewidth = 5, zorder =-20 )

    ax.arrow(x0Beam, -125, 0, 50, head_width=5, head_length=5, fc='gray', ec='gray', linewidth = 5)
    ax.text(x0Beam,-140, "BEAM", c = 'gray', ha='center', va='bottom')
    #ax.plot([-x,x],[y,-y], c= 'gray')
    #ax.scatter(0, 0, color='red', s=15, marker = '+', linewidth=30, alpha = 0.7)  # big point at the center

    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.xlabel('X(geometre) [mm]')
    plt.ylabel('Y(geometre) [mm]')
    plt.title(rf'Shift Target = ({x0},{y0}) mm, $\theta$(target) = {-angleTarget/np.pi*180} °')
    plt.grid(True)
    plt.tight_layout()
    #plt.savefig("test_no_shift.pdf")
    plt.show()


def calcul():
    N_ref = np.array([803420,801981,827156,880858])
    d2 = np.sqrt(1/N_ref)
    d2 *= 78.35/d2[0]
    for i in range(4):
        print(f'For Det{i}, dist = {d2[i]} mm ')
        print(f'For Det{i}, N(det3)/Ndet[{i}] = {N_ref[-1]/N_ref[i]}  ')
        print('\n')

    N_mes = np.array([392844,325642,298612,472568])
    d_mes = np.array([78.35,86.52,77.90,74.29])

    #print(N_mes*d_mes**2/(np.max(N_mes*d_mes**2)))
plot_config([-0.1,6],0,-6)
#calcul()
