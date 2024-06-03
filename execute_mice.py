import subprocess

def run_rscripts():
    filepath_common = 'proturn_output_'
    strains = ['aj', 'balbc', 'c57', 'cej', 'dba', 'fvb']
    exp = ['ctrl', 'iso'] 
    hl = 'hl.out'
    dat = 'hl-data_masked'
    mask = [1,2,3,4]
    commands = []
    for s in strains: 
        for e in exp: 
            for i in mask: 
                commands.append(f'Rscript mice.R {i} {filepath_common}{s}_{e}/{hl} {filepath_common}{s}_{e}/{dat}_{i}.csv')
                
    # commands = [
    #     "Rscript mice.R 1 proturn_output_aj_ctrl/hl.out proturn_output_aj_ctrl/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_aj_ctrl/hl.out proturn_output_aj_ctrl/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_aj_ctrl/hl.out proturn_output_aj_ctrl/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_aj_ctrl/hl.out proturn_output_aj_ctrl/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_aj_iso/hl.out proturn_output_aj_iso/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_aj_iso/hl.out proturn_output_aj_iso/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_aj_iso/hl.out proturn_output_aj_iso/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_aj_iso/hl.out proturn_output_aj_iso/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_balbc_ctrl/hl.out proturn_output_balbc_ctrl/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_balbc_ctrl/hl.out proturn_output_balbc_ctrl/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_balbc_ctrl/hl.out proturn_output_balbc_ctrl/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_balbc_ctrl/hl.out proturn_output_balbc_ctrl/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_balbc_iso/hl.out proturn_output_balbc_iso/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_balbc_iso/hl.out proturn_output_balbc_iso/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_balbc_iso/hl.out proturn_output_balbc_iso/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_balbc_iso/hl.out proturn_output_balbc_iso/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_c57_ctrl/hl.out proturn_output_c57_ctrl/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_c57_ctrl/hl.out proturn_output_c57_ctrl/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_c57_ctrl/hl.out proturn_output_c57_ctrl/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_c57_ctrl/hl.out proturn_output_c57_ctrl/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_c57_iso/hl.out proturn_output_c57_iso/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_c57_iso/hl.out proturn_output_c57_iso/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_c57_iso/hl.out proturn_output_c57_iso/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_c57_iso/hl.out proturn_output_c57_iso/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_cej_ctrl/hl.out proturn_output_cej_ctrl/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_cej_ctrl/hl.out proturn_output_cej_ctrl/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_cej_ctrl/hl.out proturn_output_cej_ctrl/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_cej_ctrl/hl.out proturn_output_cej_ctrl/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_cej_iso/hl.out proturn_output_cej_iso/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_cej_iso/hl.out proturn_output_cej_iso/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_cej_iso/hl.out proturn_output_cej_iso/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_cej_iso/hl.out proturn_output_cej_iso/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_dba_ctrl/hl.out proturn_output_dba_ctrl/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_dba_ctrl/hl.out proturn_output_dba_ctrl/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_dba_ctrl/hl.out proturn_output_dba_ctrl/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_dba_ctrl/hl.out proturn_output_dba_ctrl/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_dba_iso/hl.out proturn_output_dba_iso/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_dba_iso/hl.out proturn_output_dba_iso/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_dba_iso/hl.out proturn_output_dba_iso/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_dba_iso/hl.out proturn_output_dba_iso/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_fvb_ctrl/hl.out proturn_output_fvb_ctrl/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_fvb_ctrl/hl.out proturn_output_fvb_ctrl/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_fvb_ctrl/hl.out proturn_output_fvb_ctrl/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_fvb_ctrl/hl.out proturn_output_fvb_ctrl/hl-data_masked_4.csv",
    #     "Rscript mice.R 1 proturn_output_fvb_iso/hl.out proturn_output_fvb_iso/hl-data_masked_1.csv",
    #     "Rscript mice.R 2 proturn_output_fvb_iso/hl.out proturn_output_fvb_iso/hl-data_masked_2.csv",
    #     "Rscript mice.R 3 proturn_output_fvb_iso/hl.out proturn_output_fvb_iso/hl-data_masked_3.csv",
    #     "Rscript mice.R 4 proturn_output_fvb_iso/hl.out proturn_output_fvb_iso/hl-data_masked_4.csv"
    # ]

    for cmd in commands:
        subprocess.run(cmd, shell=True)

# Call the function to execute the commands
run_rscripts()
