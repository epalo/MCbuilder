# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
import argparse, os, sys, UserInteraction


#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""
    
    # obtaining fasta and pdb files 
    
    fasta_files, pdb_files = UserInteraction.getUserInput()

# TODO: insert case of empty fasta file 
    seq_record_list = []
    for seq in fasta_files:
        for seq_record in SeqIO.parse(seq, "fasta"):
            seq_record_list.append(seq_record)

# TODO: insert case of empty pdb-file 
    parser = PDB.PDBParser()
    interact_structure = []
    for pdb_struct in pdb_files:
        interact_structure.append(parser.get_structure(pdb_struct,pdb_struct))

    pdb_seq = []
    for i in range(len(pdb_files)):
        for record in SeqIO.parse(pdb_files[i], "pdb-seqres"):
            # saves the record together with the index of the pdb file
            pdb_seq.append([record,i])

# find the sequences that occur multiple times in pdb files and save all proteins for each structural aln in a separate list
similar_seq = []

for i in range(len(pdb_seq)):
    for m in range(i):
        print(pdb_seq[i][1])
        # just check sequence alignments if sequences are not in the same pair
        if not(pdb_seq[i][1] == pdb_seq[m][1]):
            # find the best alignment for two sequences (first element of list of alignments)
            alignment = pairwise2.align.globalxx(pdb_seq[i][0].seq, pdb_seq[m][0].seq)[0]
            aln_seq_1 = alignment[0]
            aln_seq_2 = alignment[1]
            al_length = len(alignment[0])
            ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))
    
            if ident/al_length >= 0.95:
                for elem in similar_seq:
                    if pdb_seq[i][1] in elem:
                        elem.append(pdb_seq[m][1])
                    elif pdb_seq[m][1] in elem:
                        elem.append(pdb_seq[i][1])
                    else: 
                        similar_seq.append([pdb_seq[i][1], pdb_seq[m][1]])
print(similar_seq)


# do a structural alignment for all pdb files that contain an identity higher than 95%
def calcStrucAln():
    for i in range(len(similar_seq)):
        f = open('{}.domains'.format(i),"w+")
        for elem in similar_seq[i]:
            #access to pdb_file with certain index
            f.write(pdb_files[elem])
        f.close
        # do structural superposition with all pdb files that the index refers to
        # for using STAMP we need like globin.domains file

        #./2hhb.pdb 2hhba {CHAIN A}
        #./2hhb.pdb 2hhbb {CHAIN B}
        #./1lh1.pdb 1lh1 {ALL}
        #./2lhb.pdb 2lhb {ALL}
        #./4mbn.pdb 4mbn {ALL}
        #./1ecd.pdb 1ecd {ALL}



        # create file that contain all of the pdb files with the indexes

        # f= open("guru99.txt","w+")
        # for i in range(10):
        #   f.write("This is line %d\r\n" % (i+1))
        # f.close() 

        # install STAMP in our program?? as a dependency
        # then run STAMP


        return
