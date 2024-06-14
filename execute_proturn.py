

import subprocess

def run_rscripts():
    filepath_common = 'proturn_output_'
    strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
    exp = ['ctrl', 'iso'] 
    hl = 'hl.out'
    dat = "dmi_np_"
    ids = "IDs_"
    mask = [1,2,3,4,5]
    commands = []
    for s in strains: 
        for e in exp: 
            for i in mask: 
                
                #example : 
                #Rscript proturn_dmi.R proturn_output_aj_ctrl hl.out dmi_np_1.csv IDs_1.csv 1
                
                commands.append(f"Rscript proturn_dmi.R {filepath_common}{s}_{e} {hl} {dat}{i}.csv {ids}{i}.csv {i}")
                          
    for cmd in commands:
        print(cmd)
        subprocess.run(cmd, shell=True)

# Call the function to execute the commands
run_rscripts()