import myBio as bio

def translate(seq,codeTable) :
    """
    This function translate a Dna sequence into a proteic sequence.
    Arg : seq - the gene sequence
          codeTable - the genetic code associated with the species.
    Returns : The proteic sequence as a string.
    """
    seqprot=""
    flag=bio.isDNA(seq)
    if flag==True:
        for i in range(len(seq)):
            codon=bio.oneWord(seq,i,3)
            for tablecodon in table.keys:
                if codon==tablecodon:
                    seqprot=seqprot+table[tablecodon]
        return seqprot
    else:
        return "error, sequence is not dna."

def findORF (seq,threshold,codeTable):
    """
    This function return a list of ORFs in the form of a dictionary.
    Author : Thomas Blanc
    Arg : seq - the gene sequence
          threshold - the minimum ORF lsngth expressed in base pairs.
          codeTable - the genetic code associated with the species.
    Returns : A dictionary containing a list of information for each ORFs,
    namely their start and stop positions, their reading frames, their length ad the translated protein sequence.
    """
    startdic={}
    stopdic={}

#Research of all the start and stop codons in the 3 frames.
    for i in len(seq):
        strt=bio.isCodonStart(seq,i)
        if strt==True:
            frame=i%3
            pos=i
            startdic[i]=frame

        stp=bio.isCodonStop(seq,i)
        if stp==True:
            frame=i%3
            pos=i
            stopdic[pos]=frame

#dico[clé][3]
startframelist=[]
startposlist=[]
stopframelist=[]
stopposlist=[]

for i in startdic.keys:
    startframelist.append(stopdic[i])
    startposlist.append(i)



for i in stopdic.keys:
    stopframelist.append(stopdic[i])
    stopposlist.append(i)

 position (in bp)
stop po
ORFs={}
finalstartpos=[]
finalstoppos=[]
finalframe=[]
finallength=[]
finaltranslation=[]
cpt=0
for i in range(len(stopposlist):
    for j in range (i+1,len(poslist)):
        if stopframelist[i]==stopframelist[j]:

            for k in range (len(startposlist)):
                if startposlist[k]>stopposlist[i] and startposlist[k]<stopposlist[j] and startframelist[k]==stopframelist[i]:
                    #extract sequence
                    extractseq=seq[startposlist[k]:stopposlist[j]]

                    if len(extractseq)>threshold:
                        #translate
                        protein=translate(extractseq,codeTable)


                        #store all the useful data in lists
                        finalstartpos.append(startposlist[k])
                        finalstoppos.append(stopposlist[j])
                        finalframe.append(stopframelist[i])
                        finallength.append(len(extractseq))
                        finaltranslation.append(protein)
                        cpt=cpt+1

#store data from the lists into a dictionary
    for orf in range (0,cpt):
        ORFs[orf]=[finalstartpos[orf],finalstoppos[orf],finalframe[orf],finallength[orf],finaltranslation[orf]]
    return ORFs


def getLengths(ORFs):
    """Short description : From an ORF list return a list of ORF lengths

    this function is written by Michele RAFFAELLI .

    Args:
        ORFs - dictionary retuned by the findORF() function.

    Returns:
        The function getLengths returns a list of ORF lengths.
    """

    orf_list=[]
    for i in ORFs.keys:
        orf_list.extend(ORFs[i][3])

    return orf_list

def getLongestORF(orf_list):
     """Short description : From an list of  ORF length return the longest ORF

    this function is written by Michele RAFFAELLI .

    Args:
        orf_list: list of orf lengths returned by the function getLengths()

    Returns:
        The function getLongestORF returns the longest ORF
    """
    max=0
    for i in range (len(orf_list)):
        if(orf_list[i]>max):topValue = getTopLongestORF(orflist,value)
print topValue

            max=orf_list[i]

    return max



def getTopLongestORF(orf_list,percentage):
    """
    Short description:  auteur: DUVAL Ludwig
    Args: orflist (liste des ORF) et value: represente le pourcentage (exprimer entre 0.0 et 1.0)
    Returns: la liste des valeur
    """
    topValue = []
    nbmax = len(orf_list)*percentage
    tmp = 0
    for i in range(len(orf_list)):
        if len(topValue)<nbmax :
            topValue.append(orf_list[i])
        elif len(topValue) >= nbmax :
            for j in range(len(topValue)):
                if tmp > topValue[j] :
                    tmp,topValue[j] = topValue[j],tmp
                if orf_list[i] > topValue[j]:
                    tmp = topValue[j]
                    topValue[j] = orf_list[i]
                    break

    return topValue
