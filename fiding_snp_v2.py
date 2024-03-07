from Bio import SeqIO
import pandas as pd
from pathlib import Path
import xmltodict
import subprocess
import argparse
import os
import re


#P4=0; P3=1; P2=2; P1=3;P1'=4

POS_DICTIONARY = {0: "P4", 1: "P3", 2: "P2", 3: "P1", 4: "P1'"}

AAS_CODE = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'Q': 'Gln',
    'E': 'Glu',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val'
}

def snp_finder(csv, main_dicretory):
    '''
    Main function
    '''

    main_df = pd.read_csv(csv)
    filtered_dfs = []

    for i, gene in enumerate(main_df.genes):
        fasta_file = SeqIO.parse(open(f"{main_dicretory}/input/proteins/{main_df.uniprot_accession.iloc[i]}.fa"), "fasta")
        motif = main_df.motifs.iloc[i]  
        #TODO: Não refazer a pesquisa quando for o mesmo gene, com motivos diferentes
        df1 = search(gene=gene)
        if df1.empty:
            continue
        df2 = xml_processing(f"snps_{gene}.xml", gene=gene)
        filtered_df = process_data(df1, df2, fasta_file, motif)
        filtered_dfs.append(filtered_df)

        print(f"{gene} SNPs have been analyzed")

    final_df = aggregator(filtered_dfs)

    return final_df

def aggregator(dataframes):
    '''
    This function concatenates all the SNPs dataframes from
    different genes
    '''
    
    aggregated_df = pd.concat(dataframes).reset_index(drop=True)

    return aggregated_df

def process_data(df1, df2, fasta_file, motif):

    '''
    Function to merge the dataframes that contain some biological
    information of SNPs and the dataframe that contains allele frequencies.
    It also processes biological annoations.
    '''
    # Setting offset, protein_id, position_range, three letter code for motif AAs 
    # and regex pattern. All of this is important for HGVS searching
    offset = []
    protein_id = []
    for f in fasta_file:
        offset.append((int(str(f.seq).find(motif)) + 1))
        protein_id.append(f.id)
    
    positions_range = [offset[0] + n for n in range(0, 5)]
    motif_code = [AAS_CODE.get(code, '') for code in motif]

    patterns = []
    for code, pos in zip(motif_code, positions_range):
        pattern = protein_id[0] + ":p." +  code + str(pos)
        patterns.append(pattern)

    pattern_regexes = [re.compile(re.escape(pattern) + r"(?!\d)") for pattern in patterns]
    print(pattern_regexes)
    
    df_to_filter = df1.merge(right=df2, 
                            left_on=["SNP_ID"], 
                            right_on=["snps"])   
    
    # Filtering only GnomAD studies
    df_to_filter = df_to_filter[df_to_filter['studies'].isin(['GnomAD', 'GnomAD_exomes'])]
    

    # Filtering by patterns (forma 1)
    matching_indices = [
    i
    for i in range(len(df_to_filter))
    for hgvs in df_to_filter.DOCSUM.iloc[i].split(",")
    if any(p_regex.search(hgvs) for p_regex in pattern_regexes) #We can have more than two matches for multiallelic
    ]
    

    filtered_df = df_to_filter.iloc[matching_indices]

    # Assigning the cleavage position
    cleavage_positions = []
    allele_class = []
    hgvs_notation = []
    bi_hgvs = []
    bi_cp = []
    bi_variants = pd.DataFrame(columns=filtered_df.columns)
    for i in range(len(filtered_df)):
        count = 0
        for hgvs in filtered_df["DOCSUM"].iloc[i].split(","):
            for p_regex in pattern_regexes:
                if p_regex.search(hgvs):
                    if count >= 1:
                    	#criar novo df com variantes bi ou multialélica
                    	bi_variants = bi_variants.append(filtered_df.iloc[i], ignore_index=True)
                    	bi_hgvs.append(hgvs)
                    	bi_cp.append(POS_DICTIONARY[pattern_regexes.index(p_regex)])
                    	print(f"One rsID multiallelic - {filtered_df.gene.iloc[i]}")
                    else:
                        cleavage_positions.append(POS_DICTIONARY[pattern_regexes.index(p_regex)])
                        hgvs_notation.append(hgvs)
                        count = count + 1
    #Completando o novo dataframe
    bi_variants.loc[:, "cleavage_position"] = bi_cp
    bi_variants.loc[:, "hgvs_notation"] = bi_hgvs
    
    #Adicionando as linhas no dataframe normal
    filtered_df.loc[:, "cleavage_position"] = cleavage_positions
    filtered_df.loc[:, "hgvs_notation"] = hgvs_notation
    #filtered_df.loc[:, "allele_classification"] = allele_class
    
    #concatenate
    if not bi_variants.empty:
    	filtered_df = pd.concat([filtered_df, bi_variants])
    if not filtered_df.empty:
    	filtered_df.loc[:, "motif"] = motif
    print(matching_indices)
    return filtered_df

def search(gene):
    '''
    This function search for all the missense
    mutations that exists at certain gene. The xml is processed and
    the function returns a dataframe containg snp_id, annotations
    and clinical significance (if exists). 
    '''
    ...

    error = False
    #TODO: usar pathlib ao invés de muitas strings
    if not os.path.exists(f"../output/{gene}/snps_{gene}.xml"):
        subprocess.run(["mkdir", f"../output/{gene}"])
        command = f'esearch -db snp -query \'({gene}[GENE]) AND "missense variant"[Function Class] AND "gnomad"[PROJECT]\' | efetch -format xml > ../output/{gene}/snps_{gene}.xml'

        subprocess.run(command, shell=True)
    try:
        data = pd.read_xml(f"../output/{gene}/snps_{gene}.xml")
    #Eception for when the XML file is truncated
    except Exception as e:
        error_message = str(e)
        print(e)
        line_number = int(error_message.split("line ")[1].split(",")[0])
        lines_to_exclude = [line_number - 2, line_number - 1]
        with open(f"../output/{gene}/snps_{gene}.xml", "r") as file:
            lines = file.readlines()
        filtered_lines = [line for idx, line in enumerate(lines) if idx not in lines_to_exclude]
        with open(f"../output/{gene}/snps_{gene}.xml", "w") as file:
            file.writelines(filtered_lines)
        try:
            data = pd.read_xml(f"../output/{gene}/snps_{gene}.xml")
        except Exception as e:
            error = True
            
    
    if error:
        main_data = pd.DataFrame()
        print(f"ERROR:{gene}")
    else:
        #Handle when there is no "clinical_significance" column
        if not "CLINICAL_SIGNIFICANCE" in data.columns:
            main_data = data[["SNP_ID", "DOCSUM"]]
        else:
            main_data = data[["SNP_ID", "DOCSUM", "CLINICAL_SIGNIFICANCE"]]
        main_data["SNP_ID"] = main_data["SNP_ID"].apply(str)
        main_data.loc[:, "gene"] = gene

        main_data.to_csv(f"../output/{gene}/snps_ids_{gene}.csv")

    return main_data

def xml_processing(snp_data, gene=""):
    
    '''
    This function processes the XML returned from NCBI and return a 
    dataframe with SNP_IDs, the studies where they were intentified
    and their respective allele frequencies.
    '''
    with open(f'../output/{gene}/{snp_data}', 'r', encoding='utf-8') as f:
        snp_xml = f.read()

    snp_dict = xmltodict.parse(snp_xml)

    studies = list()
    freqs = list()
    allele_counts = []
    snps = []
    
    for n_th_snp in range(len(snp_dict['ExchangeSet']['DocumentSummary'])):
        #TODO: recuperar o GENE também
        
        snp = snp_dict['ExchangeSet']['DocumentSummary'][n_th_snp]["SNP_ID"]

        if "GLOBAL_MAFS" in snp_dict['ExchangeSet']['DocumentSummary'][n_th_snp]:
            counter = 0
            if not type(snp_dict['ExchangeSet']['DocumentSummary'][n_th_snp]["GLOBAL_MAFS"]['MAF']) is dict:
                for n_th_study in range(len(snp_dict['ExchangeSet']['DocumentSummary'][n_th_snp]["GLOBAL_MAFS"]['MAF'])):
                    data = snp_dict['ExchangeSet']['DocumentSummary'][n_th_snp]["GLOBAL_MAFS"]['MAF'][n_th_study]

                    study = data["STUDY"]
                    allele, maf, ac = re.split('=|/', data["FREQ"])

                    studies.append(study)
                    freqs.append(float(maf))
                    allele_counts.append(int(ac))

                    counter = counter + 1
                snps = snps + [f"{snp}"] * counter
            #This condition is needed for SNPs that have only one study frequency
            else: 
                data = snp_dict['ExchangeSet']['DocumentSummary'][n_th_snp]["GLOBAL_MAFS"]['MAF']

                study = data["STUDY"]
                allele, maf, ac = re.split('=|/', data["FREQ"])

                studies.append(study)
                freqs.append(float(maf))
                allele_counts.append(int(ac))

                snps = snps + [f"{snp}"] * 1
        else:
            continue

    print(len(snps), len(studies), len(freqs))
    snp_df = pd.DataFrame({"snps": snps,
                  "studies": studies,
                  "freqs": freqs,
                  "allele_count": allele_counts})
    
    snp_df["snps"] = snp_df["snps"].apply(str)
    
    snp_df.to_csv(f"../output/{gene}/snps_freqs_{gene}.csv")
    
    return snp_df

def init_argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument("-main_dir", "--main_directory", nargs=1)
    parser.add_argument("-motifs", "--motifs_csv", nargs=1)


    return parser

def main():

    parser = init_argparser()
    args = parser.parse_args()
    #Fazer merge para saber o nome da proteína
    final_df = snp_finder(args.motifs_csv[0], args.main_directory[0])

    final_df.to_csv("../output/gnomad_report_v6.csv")


if __name__ == '__main__':
    main()

