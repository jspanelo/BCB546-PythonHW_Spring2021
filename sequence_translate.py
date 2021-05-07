#######################################
#####                             #####
#####   Python Translate Script   #####
#####                             #####
#######################################

# Libraries 

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

######################
##### Functions  #####
######################

# Get sequences

def get_sequences_from_file(fasta_fn):
    """Descrption of get_sequences_from_file() function
    
    Creates a dictionary.
    Populates the dictionary iterating trough a FASTA file.
        From each identifier, gets the species name as a key, and the nucleotide sequence as a value.
    
    Parameters:
        fasta_fn : a FASTA file, 
                        * a header line starting with a '>' followed by descriptions
                        * string of nucleotides starting in the second line
    Return:
        A dictionary of sequences.
                        * dict.key() as 'species_name' parses the second and the third element from the header line
                        * dict.values() as 'record.seq' is a Bio.Seq element which is the dna sequence for each identifier   
    """
    sequence_data_dict = {}
    for record in SeqIO.parse(fasta_fn, "fasta"):
        description = record.description.split()
        species_name = description[1] + " " + description[2]
        sequence_data_dict[species_name] = record.seq
    return(sequence_data_dict)


# 2. Translate function

def translate_function(dna_string):
    
    """Description of translate_function 
    
    Translates a DNA string into an aminoacid string.
    
    Parameters:
        dna_string : a string of nucleotides
    
    Return: aa sequence for the protein coded in dna_string, as a string
    
    Example of usage:
    
        >> translate_function('GCGGCTTCATAGGAG')
        
        Output: 
            'AASE'
   
    """
    
    aa_seq_string = "" # start aa string as an empty string
    
    for i in range(0, len(dna_string), 3):  # iterate through the string in steps of 3
        
        codon = dna_string[i: i + 3]        # subset the string from the element i to i+3, send it to `codon`
        
        codon_seq = Seq(codon)              # turns `codon` into a Seq() object
        
        # use df.translate() to translate codon into aminoacid (aa)  
        aa = codon_seq.translate(table = "Vertebrate Mitochondrial", to_stop = True)

        aa_seq_string += aa                 # add the last translated aa to the aa string first defined
    
    return str(aa_seq_string)               # return the aa sequence as a string

# 3. Translate function

from Bio.Seq import Seq
def translate_to_protein(dna_string):
    """Description of translate_to_protein function 
    
    Translates a DNA string into a aminoacid string when the appropriate datatype is provided, or give an error message.
    
    Parameters:
        string_nucleotides : a string of nucleotides
    
    Return: aa sequence for the protein coded in string_nucleotides
            Error messagge when an invalid parameter in used as input 
    
    Example of usage:
    
        >> translate_to_protein('GCGGCTTCATAGGAG')
        
        Output: 
            AAS*E
   
    """
    if type(dna_string) == str:
        
        aa_seq = Seq(dna_string).translate()
        
        return print(aa_seq)
    else:
        
        return print('Error: Invalid object')

# 4. Molecular weight 

def trans_func( dna_string ):
    """ Decription for translation formula
    
    Translate a DNA string into a aminoacid string, and returns its molecular weight 
    
    Parameters:
        dna_string : a string of nucleotides 
    
    Return:
        molecular weight for the encoded protein
    
    Example of usage:
        >> trans_func("ACGACCGGA")
        
        Output: 
            277.2744
        
    """
    coding_dna = Seq( dna_string )
    
    aa_string = coding_dna.translate( table = "Vertebrate Mitochondrial", to_stop = True) 
    
    analysed_seq = ProteinAnalysis( str( aa_string ) )
    
    molecular_weight = analysed_seq.molecular_weight()
    
    return molecular_weight

# 5. CG content calculator

def CG_content(string_nucleotides):
    """Description CG content calculator
    
    Calculates the CG content % for a given DNA string
    
    Parameters:
        string_nucleotides : a string of nucleotides or a Seq() object
    
    Return: CG content in the string as a percentage
    
    Example of usage:
        Using a string
        >> CG_content("ACGACCGGAC")
        
        Output: 
            70.0
   
        Using a Seq() object  
        >> ms = Seq("ACGGCTTACCGGAC")
        >> CG_content(ms)
        
        Output: 
            64.28571428571429   
                
    """
    nuc_string = str(string_nucleotides) # force string or Seq() object as a string

    seq_up = nuc_string.upper()          # force string to be uppercase
    
    my_seq = Seq(seq_up)                 # turn string into a seq() object
    
    # calculate CG content
    cg_count = 100*((my_seq.count("G") + my_seq.count("C")) /len(my_seq))
    
    return(cg_count) # return the percent of CG


#####################
####### Main  #######
#####################

cytb_seqs = get_sequences_from_file("bears_cytb.fasta") 

bears_df = pd.read_csv("bears_mass.csv") # Includes only data for body mass 

species_list = list(bears_df.species)    # get the list of species in the dataframe

# 6. Add columns for CG content and MW

bears_df['Mol weight'] = 'NaN' # Creates a new column setting a common value to each element in the column

bears_df['CG content'] = 'NaN'

# 7. For-loop to that translates each sequence, gets molecular weight and computes the GC content

bears_pd = bears_df.set_index('species') # set 'species' as index to facilitate iteration  

i = 'NULL'                               # set iterative variable as NULL 

for i in cytb_seqs.keys():               # iterate through each key (species) in the dictionary cytb_seqs  

    dna_string = cytb_seqs[i]               # get the value (dna sequence, string) for the ith key in the cytb_seqs, and call it 'dna_string' 
    
    MW = trans_func(str(dna_string))        # apply trans_func() to dna_string. Get molecular weight (float) for i and stored it as MW 
    
    CG = CG_content(dna_string)             # apply GC_content() to dna_string. Get CG content (float) for i and stored it as CG 
    
    bears_pd.at[i, 'CG content'] = CG       # replace 'CG content' for index (specie) i, by the ith float stored as CG
    
    bears_pd.at[i, 'Mol weight'] = MW       # replace 'Mol weight' for index (specie) i, by the ith float stored as MW
    
bears_pd = bears_pd.reset_index()

# 8. Bar-plot: mass by species.

weights = sns.barplot(
            data = bears_pd,
            x = 'species', 
            y = "mass", 
            order = bears_pd.sort_values('mass').species, # sort species by weight 
            palette = "colorblind")

plt.xticks(rotation = 45,                   # rotate labels for bars
           horizontalalignment = 'right')   # align labels properly below each bar
plt.xlabel('Species')
plt.ylabel("Body mass [Kg]")
plt.title("Bears body mass by species")
weights

# 9. Molecular Weight (y-axis) as a function of GC-content (x-axis)

CG_to_MW = sns.scatterplot(data = bears_pd,
                           x = "CG content",
                           y = "Mol weight",
                           hue = 'species',
                           palette = "colorblind")

CG_to_MW.legend(title = 'Species', loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('CG content [%]')
plt.ylabel("Molecular weight")
plt.title("Cytochrome b: Molecular weight vs CG content")

# 10. Save the new DataFrame as "bears_mass_cytb.csv"

bears_pd.to_csv("bears_mass_cytb.csv", index = False) 