"""Add CHARGE=1 after every PEPMASS line of the merged MassQL MGFs (input for Falcon)."""

# Root of the local project folder; edit to run elsewhere.
PROJECT_DIR = '/media/dorrestein/Helena2/Bile_acids_update'


input_files = [
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Nonhydroxy/download_mgfs/final_mgf/C24_Nonhydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Monohydroxy/download_mgfs/final_mgf/C24_Monohydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Dihydroxy/download_mgfs/final_mgf/C24_dihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Trihydroxy/download_mgfs/final_mgf/C24_Trihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Tetrahydroxy/download_mgfs/final_mgf/C24_Tetrahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_final/Pentahydroxy/download_mgfs/final_mgf/C24_Pentahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_monohydroxy/download_mgfs/final_mgf/C23_Monohydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_dihydroxy/download_mgfs/final_mgf/C23_dihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_trihydroxy/download_mgfs/final_mgf/C23_trihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_tetrahydroxy/download_mgfs/final_mgf/C23_tetrahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C23/C23_pentahydroxy/download_mgfs/final_mgf/C23_pentahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_monohydroxy/download_mgfs/final_mgf/C27_monohydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_dihydroxy/download_mgfs/final_mgf/C27_dihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_trihydroxy/download_mgfs/final_mgf/C27_trihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_tetrahydroxy/download_mgfs/final_mgf/C27_tetrahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C27/C27_pentahydroxy/download_mgfs/final_mgf/C27_pentahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/Aceto/Dihydroxy/download_mgfs/final_mgf/Aceto_dihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_monohydroxy/download_mgfs/final_mgf/C24_alkylamine_monohydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_dihydroxy/download_mgfs/final_mgf/C24_alkylamine_dihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_trihydroxy/download_mgfs/final_mgf/C24_alkylamine_trihydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_tetrahydroxy/download_mgfs/final_mgf/C24_alkylamine_tetrahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C24_alkylamine/C24_alkylamine_pentahydroxy/download_mgfs/final_mgf/C24_alkylamine_pentahydroxy_input_falcon.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C28_alcohol/One_OH/download_mgfs/final_mgf/C28_alcohol_1OH_input_falcon_351_C28O1.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C28_alcohol/Two_OH/download_mgfs/final_mgf/C28_alcohol_2OH_input_falcon_367_C28O2.mgf",
    f"{PROJECT_DIR}/Queries_final_mgf/C28_alcohol/Three_OH/download_mgfs/final_mgf/C28_alcohol_3OH_input_falcon_383_C28O3.mgf"
]

for input_file in input_files:
    output_file = input_file.replace(".mgf", "_charge.mgf")
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            outfile.write(line)
            if line.startswith("PEPMASS="):
                outfile.write("CHARGE=1\n")
    
    print(f"Processed: {input_file}")