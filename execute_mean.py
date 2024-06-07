import subprocess

def run_rscripts():
    filepath_common = 'proturn_output_'
    strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
    exp = ['ctrl', 'iso'] 
    dat = 'hl-data_masked'
    mask = [1,2,3,4,5]
    commands = []
    for s in strains: 
        for e in exp: 
            for i in mask: 
                commands.append(f"Rscript mean.R {i} {filepath_common}{s}_{e}/{dat}_{i}.csv")
                         
    for cmd in commands:
        subprocess.run(cmd, shell=True)

# Call the function to execute the commands
run_rscripts()