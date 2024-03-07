from Bio import SeqIO
fasta_sequences = SeqIO.parse(open("../proteins.fasta"),'fasta') 

my_record = []

for fasta in fasta_sequences:
    fasta.id = fasta.description.split("|")[1]
    fasta.description = ""
    output_file=f"{fasta.id}.fa"
    my_record.append(fasta)
    SeqIO.write(my_record, output_file, "fasta")

    my_record.clear()
