"""Renumber SCANS sequentially (1..n) in an MGF."""

# Root of the local project folder; edit to run elsewhere.
PROJECT_DIR = '/media/dorrestein/Helena2/Bile_acids_update'

in_file = f'{PROJECT_DIR}/Queries_final_mgf/Aceto/Dihydroxy/download_mgfs/final_mgf/Aceto_dihydroxy_input_falcon_charge.mgf'

with open(in_file, 'r') as in_f:
    with open(f'{PROJECT_DIR}/Queries_final_mgf/Aceto/Dihydroxy/download_mgfs/final_mgf/Aceto_dihydroxy_input_falcon_charge_renumbered.mgf',
              'a') as out_f:
        counter = 1
        for line in in_f:
            if line.startswith('SCANS'):
                print('SCANS=' + str(counter), file=out_f, end='\n')
                counter += 1
            else:
                print(line, file=out_f, end='')