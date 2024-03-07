import pandas as pd
import subprocess
import sys

genes = pd.read_csv(f"{sys.argv[1]}/gene_table_of.csv")

for i,gene in enumerate(genes.genes):
    command= f'esearch -db protein -query "{gene}[GENE] AND MANE_select[KYWD]"| efetch -format fasta > {sys.argv[1]}/proteins/{genes.uniprot_accession.iloc[i]}.fa'

    subprocess.run(command, shell=True)
    print(f'{gene} downloaded')