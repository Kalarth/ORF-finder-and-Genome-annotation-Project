
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
