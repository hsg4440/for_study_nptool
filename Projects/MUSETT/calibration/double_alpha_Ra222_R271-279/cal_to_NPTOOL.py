import matplotlib.pyplot
import re
MapX = [33, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 63,62, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 32,1,  2,  4,  6,  8,  10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 31,30, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11,  9,  7,  5,  3,  0,97, 98,100,102,104,106,108,110,112,114,116,118,120,122,124,127,126,125,123,121,119,117,115,113,111,109,107,105,103,101, 99, 96,65, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 95,94, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 64]


MapY =[30, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11,  9,  7,  5,  3,  0,1,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 31,62, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 32,33, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 63,94, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 64,65, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 95,126,125,123,121,119,117,115,113,111,109,107,105,103,101,99,  96,97, 98,100,102,104,106,108,110,112,114,116,118,120,122,124,127]

Ecal = 1
Tcal = 1
clean = 1
if Ecal :
    if clean == 1 :
        for mumu in range(4):
            for xy in ['X','Y']:

                input_file_path = "M" + str(mumu) + "_E/Cal_Str_" + xy + "_E_MM" + str(mumu) + ".cal"
                with open(input_file_path, 'w') as file:
                    pass
    f = open("calibE0_R271-279.cal1") #Energy calibration

    for line in f:
        temp = line.split()
        p0,p1,p2 = float(temp[0]),float(temp[1]),float(temp[2])

        input_string = temp[-1]
        pattern = r"Mu(\d+)Ch([XY])(\d+)"
        match = re.search(pattern, input_string)

        if match:
            # Extract the letter and number using group(2) and group(3)
            mu_number = match.group(1)
            X_Y = match.group(2)
            ch = match.group(3)
            if X_Y == "X":
                strip = MapX[int(ch)-1]
                #print("\n Det = ", mu_number)
                #print("ch = ", ch, ", strip =", strip)
                #print("p0 = ", p0, "p1 = ", p1)
            else :
                strip = MapY[int(ch)-1]
            strip = int(ch)-1
            #strip = int(ch)-1
        file_path = f"M{mu_number}_E/Cal_Str_{X_Y}_E_MM{mu_number}.cal"

        # Define the text to write to the file
        text_to_write = f"MUSETT_T{mu_number}_DSSD_{X_Y}{strip}_E" + " " + str(float(p0)/1000) + " " + str(float(p1)/1000) + "\n"
        # Use a try-except block to handle any potential errors
        try:
            # Open the file in write mode ('w') and write the text
            with open(file_path, 'a') as file:
                file.write(text_to_write)
        except Exception as e:
            print(f"An error occurred: {str(e)}")

    for mumu in range(4):
        for xy in ['X','Y']:

            input_file_path = "M" + str(mumu) + "_E/Cal_Str_" + xy + "_E_MM" + str(mumu) + ".cal"
            with open(input_file_path, 'r') as file:
                lines = file.readlines()

            def sort_key(line):
                # Use a regular expression to find the DSSD_Xnn number
                import re
                if xy == "X":
                    match = re.search(r'DSSD_X(\d+)_E', line)
                else :
                    match = re.search(r'DSSD_Y(\d+)_E', line)
                if match:
                    return int(match.group(1))
                return 0  # Return 0 for lines that don't match the pattern

            # Sort the lines based on the DSSD_Xnn number
            sorted_lines = sorted(lines, key=sort_key)

            # Define the path to your output file
            output_file_path = input_file_path

            # Write the sorted lines to the output file
            with open(output_file_path, 'w') as file:
                file.writelines(sorted_lines)


if Tcal :
    f = open("calibT0_R271-279.cal1") #Time calibration
    if clean == 1 :
        for mumu in range(4):
            for xy in ['X','Y']:

                input_file_path = "M" + str(mumu) + "_T/Cal_Str_" + xy + "_T_MM" + str(mumu) + ".cal"
                with open(input_file_path, 'w') as file:
                    pass
    for line in f:
        temp = line.split()
        p0,p1,p2 = float(temp[0]),float(temp[1]),float(temp[2])

        input_string = temp[-1]
        pattern = r"Mu(\d+)Ch([XY])(\d+)"
        match = re.search(pattern, input_string)

        if match:
            # Extract the letter and number using group(2) and group(3)
            mu_number = match.group(1)
            X_Y = match.group(2)
            ch = match.group(3)
            if X_Y == "X":
                strip = MapX[int(ch)-1]
            else :
                strip = MapY[int(ch)-1]
            strip = int(ch)-1
        file_path = f"M{mu_number}_T/Cal_Str_{X_Y}_T_MM{mu_number}.cal"

        # Define the text to write to the file
        text_to_write = f"MUSETT_T{mu_number}_DSSD_{X_Y}{strip}_T" + " " + str(float(p0)) + " " + str(float(p1)) + "\n"

        # Use a try-except block to handle any potential errors
        try:
            # Open the file in write mode ('w') and write the text
            with open(file_path, 'a') as file:
                file.write(text_to_write)
        except Exception as e:
            print(f"An error occurred: {str(e)}")

    for mumu in range(4):
        for xy in ['X','Y']:

            input_file_path = "M" + str(mumu) + "_T/Cal_Str_" + xy + "_T_MM" + str(mumu) + ".cal"
            with open(input_file_path, 'r') as file:
                lines = file.readlines()

            def sort_key(line):
                # Use a regular expression to find the DSSD_Xnn number
                import re
                if xy == "X":
                    match = re.search(r'DSSD_X(\d+)_T', line)
                else :
                    match = re.search(r'DSSD_Y(\d+)_T', line)
                if match:
                    return int(match.group(1))
                return 0  # Return 0 for lines that don't match the pattern

            # Sort the lines based on the DSSD_Xnn number
            sorted_lines = sorted(lines, key=sort_key)

            # Define the path to your output file
            output_file_path = input_file_path

            # Write the sorted lines to the output file
            with open(output_file_path, 'w') as file:
                file.writelines(sorted_lines)
