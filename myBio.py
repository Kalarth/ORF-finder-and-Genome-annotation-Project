
def getStandardCode():
"""This function gives the standard genetic code in the form of a dictionary. Author : Thomas Blanc
Args : None
Returns : The function returns a dictionary containing the amino acids coded by each codons (codons are used as keys.)
"""

    table={}
    base1="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    base2="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    base3="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
    AAs="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

    for i in range(0,len(AAs)):
        codon=base1[i]+base2[i]+base3[i]
        aa=AAs[i]
        table[codon]=aa
    return table

def getGeneticCode(transl_table=4):
"""Returns a dictionary with the coding table corresponding to the number"""
    table={}





#Functions used to invert the dna strand.
def brincomp(seq):
    comp = ""
    for i in range(len(seq)):
        if seq[i] == "A":
            comp = comp + "T"
        if seq[i] == "C":
            comp = comp + "G"
        if seq[i] == "G":
            comp = comp + "C"
        if seq[i] == "T":
            comp = comp + "A"
    return comp


def invert(seq):
    comp = brincomp (seq)
    print comp
    revers=comp[::-1]
    #revers = ""
    """
    for i in range(len(seq)):
        if i != 0 :
            revers = revers + comp[-i]
        if i == 0 :
    """
    print revers

    #end of reverse Functions.
