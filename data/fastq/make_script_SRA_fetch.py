#!/users/9/kell3262/miniforge3/envs/ncbi/bin/python


sra_file = "Caballo_Moro_SRAdataset.txt"

print(f'Making script to fetch fastqs from {sra_file}')

list_sras = []
with open(sra_file, "r") as file:
    for line in file:
        line = line.rstrip()
        list_sras.append(line)

prefetch_cmd = f"prefetch --option-file {sra_file}"
dump_cmd_list = [f'fasterq-dump --split-files {sra}' for sra in list_sras]
dump_cmd = '\n'.join(dump_cmd_list)

with open("fetch_fastqs.sh", "w") as write_file:
    write_file.write(f"#!/usr/bin/bash \n {prefetch_cmd}\n{dump_cmd}")
