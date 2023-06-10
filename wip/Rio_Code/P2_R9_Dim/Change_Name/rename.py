import os;
#import numpy as np
def Creat_file(Scan_start,Scan_end):
    BasicPath =  os.path.expanduser(
        "~/EnvPBGEM_Linux/SimSave/P2_R9_Dim")
    Target = "/Latin_6para_200cases_narrow/"
    purpose = "Latin_6para_200cases_40oC_narrow"
    source_file =BasicPath+Target+ f"{purpose}_1_10.py"  # Name of the source .py file
    destination_file = BasicPath+Target+ f"{purpose}_{Scan_start}_{Scan_end}.py"  # Name of the destination .py file (copy)
    old_text = "Scan_start = 1;    Scan_end = 10;"  # Text to be replaced
    new_text = f"Scan_start = {Scan_start};    Scan_end = {Scan_end};"

    # Copy the source file to the destination file
    os.system(f"cp {source_file} {destination_file}")

    # Read the contents of the destination file
    with open(destination_file, "r") as file:
        content = file.read()

    # Modify the desired parts of the content
    modified_content = content.replace(old_text, new_text)

    # Write the updated contents to the destination file
    with open(destination_file, "w") as file:
        file.write(modified_content)




    ### change .pbs file
    source_file = BasicPath+Target+f"{purpose}_1_10.pbs"  # Name of the source .pbs file
    destination_file = BasicPath+Target+f"{purpose}_{Scan_start}_{Scan_end}.pbs"  # Name of the destination .pbs file (copy)
    old_content = f"python3 {purpose}_1_10.py"  # Text to be replaced
    new_content = f"python3 {purpose}_{Scan_start}_{Scan_end}.py"  # New text to replace the old_content

    # Read the contents of the source file
    with open(source_file, "r") as file:
        content = file.read()

    # Modify the desired parts of the content
    modified_content = content.replace(old_content, new_content)

    # Write the updated contents to the destination file
    with open(destination_file, "w") as file:
        file.write(modified_content)

    print(f"Finish creating {purpose}_{Scan_start}_{Scan_end}.py")


# create a bunch of files!
# Big_start = 101; Big_end = 200; case_no=10;
Scan_start_all = [
    1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 
    101, 111, 121, 131, 141, 151, 161, 171, 181, 191]
""" (
    np.arange(Big_start,Big_end+1,case_no)
    ).tolist() """
Scan_end_all = [
    10, 20, 30, 40, 50, 60, 70, 80, 90, 
    100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
""" (
    np.arange(Big_start+case_no-1,Big_end+case_no,case_no)
    ).tolist() """
print(Scan_start_all)
print(Scan_end_all)
for Scan_start,Scan_end in zip(Scan_start_all,Scan_end_all):
    Creat_file(Scan_start,Scan_end)