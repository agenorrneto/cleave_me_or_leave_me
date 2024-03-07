import subprocess
import os
import pandas as pd
import urllib.request
import Bio.PDB


#PDB structures retriving from AlphaFold Protein Structure Database
#Setting the environment
#Downloading the needed files
motifs = pd.read_csv("/home/agenor/scoring_motifs_3dparsing/cleavage_sites_oficial.csv")
for i in range(len(motifs)):
    identifier = motifs.iloc[i]['uniprot accession']
    #There are some protein which have more than one motif, so I do not need to download them twice
    if os.path.isfile(f"/home/agenor/scoring_motifs_3dparsing_results/{identifier}.pdb"):
        print("Structure (PDB) already downloaded")
    else:
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{identifier}-F1-model_v2.pdb"
            request = urllib.request.Request(url)

            with urllib.request.urlopen(request) as response:
                res = response.readlines()
            with open(f"/home/agenor/scoring_motifs_3dparsing_results/{identifier}.pdb", mode="wb") as f:
                for line in res:
                    f.write(line)
                f.close()
        except:
            print(f"{identifier} (PDB) not found")
    if os.path.isfile(f"/home/agenor/scoring_motifs_3dparsing_results/{identifier}_result.rsa.csv"):
        print("DSSP analysis already done")
    else:
        #Perform DSSP assignment
        command = ["python3", "/home/agenor/scoring_motifs_3dparsing/3d_parsing.py", 
                    "-pdb", f"/home/agenor/scoring_motifs_3dparsing_results/{identifier}.pdb", 
                    "-o", f"{identifier}_result"]
        subprocess.call(command)
        

#Secondary structure and exposition profile assignment
identifiers = []
motifs_list = []
ss_list = []
rsa_profiles = []
model_conficence_values = {}
c = 0
repeat_count = 1
#Analyzing the files
for i in range(len(motifs)):
    try:
        identifier = motifs.iloc[i]['uniprot accession']
        motif = motifs.iloc[i]["motif"]
        #Retrieving the assignments
        results = pd.read_csv(f"/home/agenor/scoring_motifs_3dparsing_results/{identifier}_result.rsa.csv")
        string = ""
        for aa in list(results.pdb_aa):
            string = string + aa
        
        offset = string.index(motif)
        print(f'Motif found {string.count(motif)} times') #A flag that should always be equal to one
        #RSA
        profile = ""
        for state in list(results.iloc[offset:offset+5].rsa):
            if state >= 0.25:
                profile = profile + "E"
            elif state < 0.25:
                profile = profile + "B"
        
        #SS
        ss = ""
        for structure in list(results.iloc[offset:offset+5].structure):
            if str(structure) == "nan":
                ss += "C"
            else:
                ss += structure
        identifiers.append(identifier)
        motifs_list.append(motif)
        ss_list.append(ss)
        rsa_profiles.append(profile)

        #Parsing PDB for plldt value retrieving
        p = Bio.PDB.PDBParser()
        structure = p.get_structure('myStructureName', f"/home/agenor/scoring_motifs_3dparsing_results/{identifier}.pdb")
        plldt = [a.get_bfactor() for a in structure.get_atoms()]
        atoms = [a.get_name() for a in structure.get_atoms()]
        residues = [res.get_resname() for res in structure.get_residues()]
        if len (plldt) == len(atoms):
            print("Keep moving on")
        df = pd.DataFrame({"plldt": plldt, "atoms": atoms})
        
        plof = []
        for i in range(len(df)):
            if df.atoms.iloc[i] == "CA":
                plof.append(df.plldt.iloc[i])
        
        plldt_df = pd.DataFrame({"residue": residues, "ppldt": plof})

        if identifier in model_conficence_values:
            model_conficence_values[f'{identifier}_{repeat_count}'] = list(plldt_df.ppldt.iloc[offset:offset+5])
            repeat_count = repeat_count + 1
        else:
            model_conficence_values[identifier] = list(plldt_df.ppldt.iloc[offset:offset+5])
        c += 1
        print(f"{c} structures analyzed")
    except Exception as e:
    	print(e)



print(ss_list)
final_report = {
    "identifier": identifiers,
                    "motif": motifs_list,
                    "ss_8": ss_list,
                    "profile": rsa_profiles
                    }
pd.DataFrame(final_report).to_csv("final_report.csv")
pd.DataFrame({key: pd.Series(value) for key, value in model_conficence_values.items()}).T.to_csv("model_conficence_values.csv")





