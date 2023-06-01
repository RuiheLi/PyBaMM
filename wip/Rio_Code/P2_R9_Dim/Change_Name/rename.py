import os;
#import numpy as np
def Creat_file(Scan_start,Scan_end):
    BasicPath =  os.path.expanduser(
        "~/EnvPBGEM_Linux/SimSave/P2_R9_Dim")
    Target = "/Exp2_AAT/"
    source_file =BasicPath+Target+ "Run_case_1_10.py"  # Name of the source .py file
    destination_file = BasicPath+Target+ f"Run_case_{Scan_start}_{Scan_end}.py"  # Name of the destination .py file (copy)
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
    source_file = BasicPath+Target+"Run_case_1_10.pbs"  # Name of the source .pbs file
    destination_file = BasicPath+Target+f"Run_case_{Scan_start}_{Scan_end}.pbs"  # Name of the destination .pbs file (copy)
    old_content = "python3 Run_case_1_10.py"  # Text to be replaced
    new_content = f"python3 Run_case_{Scan_start}_{Scan_end}.py"  # New text to replace the old_content

    # Read the contents of the source file
    with open(source_file, "r") as file:
        content = file.read()

    # Modify the desired parts of the content
    modified_content = content.replace(old_content, new_content)

    # Write the updated contents to the destination file
    with open(destination_file, "w") as file:
        file.write(modified_content)

    print(f"Finish creating Run_case_{Scan_start}_{Scan_end}.py")


# create a bunch of files!
# Big_start = 101; Big_end = 200; case_no=10;
Scan_start_all = [
    201, 211, 221, 231, 241, 251, 261, 271, 281, 291, 
    301, 311, 321, 331, 341, 351, 361, 371, 381, 391, 
    401, 411, 421, 431, 441, 451, 461, 471, 481, 491, 
    501, 511, 521, 531, 541, 551, 561, 571, 581, 591, 
    601, 611, 621, 631, 641, 651, 661, 671, 681, 691, 
    701, 711, 721, 731, 741, 751, 761, 771, 781, 791, 
    801, 811, 821, 831, 841, 851, 861, 871, 881, 891, 
    901, 911, 921, 931, 941, 951, 961, 971, 981, 991]
""" (
    np.arange(Big_start,Big_end+1,case_no)
    ).tolist() """
Scan_end_all = [
    210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 
    310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 
    410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 
    510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 
    610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 
    710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 
    810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 
    910, 920, 930, 940, 950, 960, 970, 980, 990, 1000]
""" (
    np.arange(Big_start+case_no-1,Big_end+case_no,case_no)
    ).tolist() """
print(Scan_start_all)
print(Scan_end_all)
for Scan_start,Scan_end in zip(Scan_start_all,Scan_end_all):
    Creat_file(Scan_start,Scan_end)