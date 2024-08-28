import os;
import numpy as np


def Create_file(purpose_i, original, Case_No):
    
    BasicPath =  os.path.expanduser(
        "~/EnvPBGEM_NC/SimSave")
    Target = "/Jobs/" 
    

    ########### First, change text for .py file
    source_file =BasicPath+Target+ f"{purpose_i}_use_Fun_NC_case_{original}.py"  # Name of the source .py file
    destination_file = BasicPath+Target+ f"{purpose_i}_use_Fun_NC_case_{Case_No}.py"  # Name of the destination .py file (copy)
    old_text = f"i_bundle = {original} # int(os.environ"  # Text to be replaced
    new_text = f"i_bundle = {Case_No} # int(os.environ"

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
    print(f"Finish creating {purpose_i}_use_Fun_NC_case_{Case_No}.py")





    ########### Second, change text for .py file change .pbs file
    source_file = BasicPath+Target+f"{purpose_i}_use_Fun_NC_case_{original}.pbs"  # Name of the source .pbs file
    destination_file = BasicPath+Target+f"{purpose_i}_use_Fun_NC_case_{Case_No}.pbs"  # Name of the destination .pbs file (copy)
    old_content = f"python3 {purpose_i}_use_Fun_NC_case_{original}.py"  # Text to be replaced
    new_content = f"python3 {purpose_i}_use_Fun_NC_case_{Case_No}.py"  # New text to replace the old_content

    # Read the contents of the source file
    with open(source_file, "r") as file:
        content = file.read()
    # Modify the desired parts of the content
    modified_content = content.replace(old_content, new_content)
    # Write the updated contents to the destination file
    with open(destination_file, "w") as file:
        file.write(modified_content)

    print(f"Finish creating {purpose_i}_use_Fun_NC_case_{Case_No}.pbs")

    return 



original = 10
purpose_i = "SEI_Dry_Exp1235_NC"
Case_List = [220, 264]
# create a bunch of files!
for Case_No in Case_List:
    Create_file(purpose_i, original, Case_No)

