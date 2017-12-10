#-*- coding: utf-8 -*-
"""
THIS PROGRAM WAS MADE BY THOMAS BLANC, LUDVIG DUVAL, MICHELLE RAFFAELLI AND MARC MONGY.
"""
#import myBio as bio
"""
Beginning of myBio library
"""
def getStandardCode():
    """
    This function gives the standard genetic code in the form of a dictionary. Author : Thomas Blanc
    Args : None
    Returns : The function returns a dictionary containing the amino acids coded by each codons (codons are used as keys.)
    """

    codetable={}
    base1="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    base2="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    base3="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
    AAs="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    Start="---M------**--*----M---------------M----------------------------"

    for i in range(0,len(AAs)):
        codon=base1[i]+base2[i]+base3[i]
        aa=AAs[i]
        starter=Start[i]
        codetable[codon][0]=aa
        codetable[codon][1]=starter

    return codetable

def getGeneticCode(transl_table):
    """Returns a dictionary with the coding table corresponding to the number"""
    codetable={}

    if transl_table==4:
        base1="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
        base2="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
        base3="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
        AAs="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
        Start="--MM------**-------M------------MMMM---------------M------------"

    for i in range(0,len(AAs)):
        codon=base1[i]+base2[i]+base3[i]
        aa=AAs[i]
        starter=Start[i]
        codetable[codon]=[aa, starter]


    return codetable

def isDNA(seq):
    flag=""
    for nuc in seq:
        if nuc in ["A","T","G","C","a","t","g","c"]:
            flag=True
        else:
            flag=False
            return flag
    return flag

def oneWord(seq,start,wlen):
    w=seq[start:start+wlen]
    return w

def countWord(seq,word):
    cpt=0
    wlen=len(word)
    for i in range (len(seq)):
        v=oneWord(seq,i,wlen)
        if v==word:
            cpt=cpt+1
    return cpt

def isCodonStart (seq,pos,codetable):
    flag=False
    w=oneWord(seq,pos,3)
    for codon in codetable.keys():
        if w==codon and codetable[codon][1]=="M":
            flag=True
    return flag

def isCodonStop (seq,pos,codetable):
    flag=False
    w=oneWord(seq,pos,3)
    for codon in codetable.keys():
        if w==codon and codetable[codon][1]=="*":
            flag=True
    return flag

def isGene(seq):
    i=0
    ww=False
    while i<len(seq):
        w=isCodonStart(seq,i)
        if w==True:
            ww=w
        x=isCodonStop(seq,i)
        if ww==True and x==True:
            return True
        i=i+3
    return False

"""
def isGene3(seq):
    frame=[]
    for i in range (len(seq)):
        w=isCodonStart(seq,i)
        if w==True:
            frame.append(i%3)
            for j in range (i,len(seq)):
                x=isCodonStop(seq,j)
                if x==True:
                    if i%3==frame[]:
                        return
    return
"""
#seq='TGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'
#flag = isGene3(seq) # True
#print flag





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
    #comp = brincomp (seq)
    #print comp
    reverse=seq[::-1]
    return reverse

def brinAntiSens(seq):
    brinSens=seq
    inverseBrinSens=invert(brinSens)
    complementaire=brincomp(inverseBrinSens)
    return complementaire

    #end of reverse Functions.

    """
    END of myBio library
    """


def loadFASTA(nomFichier):
    rawFASTA=""
    fichier=open(nomFichier, "r")
    for ligne in fichier.readlines():
        if not ligne:
            break
        else:
            rawFASTA=rawFASTA+ligne
    fichier.close()
    return rawFASTA


def dictionaryFromFASTA(txt, moltype='auto'):
    '''Return a dictionary separating the header from the data
       Author : Marc Mongy
    '''
    #Init dictionary
    seq = {'description':'myBio sequence','data':'None','type':'auto','DB':'db','AC':'ac','ID':'id'}
    #TODO - Extract the header from the 2nd character to the end of line
    #TODO - The end of line in text corresponds to the character newline '\n'
    #TODO - Extract the data after the first newline up to the end of txt
    #TODO - Remove all the newline characters

    seq['description']=txt.split("\n")[1]
    seq['description']=seq['description'].replace('>', '')
    #seq['DB']=
    #seq['AC']=
    #seq['ID']=
    #print (seq['description'])

    seq['data'] = ''.join(txt.split("\n")[2:])
    seq['data']=seq['data'].upper()
    #print (seq['data'])

    return seq


def readFASTA(txt):
    """
    Get the DNA sequence from the dictionary into a string
    Author : Marc Mongy
    """
    sequence=dictionaryFromFASTA(txt, moltype='auto')
    seqString=sequence['data']
    print "FASTA OK"
    return seqString


def translate(seq,codetable) :
    """
    This function translate a Dna sequence into a proteic sequence.
    Author : Thomas Blanc
    Arg : seq - the gene sequence
          codetable - the genetic code associated with the species.
    Returns : The proteic sequence as a string.
    """
    seqprot=""
    flag=isDNA(seq)
    if flag==True:
        for i in range(0,len(seq),3):
            codon=oneWord(seq,i,3)
            for tablecodon in codetable.keys():
                if codon==tablecodon:
                    seqprot=seqprot+codetable[tablecodon][0]
        return seqprot
    else:
        return "error, sequence is not dna."


def getORF(seq,threshold,codetable,startposlist,startframelist,stopposlist,stopframelist,position,finalstartpos,finalstoppos,finalframe,finallength,finaltranslation):
    """
    This function looks for an ORF for one particular Stop Codon
    ARgs:seq - the gene sequence
         threshold - the minimum ORF lsngth expressed in base pairs.
         codeTable - the genetic code associated with the species.
         startposlist - the list of positions of each start codon.
         startframelist - the list of reading frames for each start codon
         stopposlist - the list of positions of each stop codon.
         stopframelist - the list of reading frames for each stop codon
         position - the index of the stop codon used to find the ORF
         finalstartpos - list containing the start position of each ORF
         finalstoppos - list containing the stop position of each ORF
         finalframe - list containing the reading frame of each ORF
         finallength - list containing the sequence length of each ORF (in nucleotides)
         finaltranslation - list containing the protein sequence of each ORF

    Returns : finalstartpos - list containing the start position of each ORF
    finalstoppos - list containing the stop position of each ORF
    finalframe - list containing the reading frame of each ORF
    finallength - list containing the sequence length of each ORF (in nucleotides)
    finaltranslation - list containing the protein sequence of each ORF
    """


    cpt=0

    for j in range (len(stopposlist)):
        if stopframelist[position]==stopframelist[j] and position < j and stopposlist[position] < stopposlist[j]:

            for k in range (len(startposlist)):
                if startposlist[k]>stopposlist[position] and startposlist[k]<stopposlist[j]:
                    if startframelist[k]==stopframelist[j]:

                        #extract sequence
                        extractseq=seq[startposlist[k]:stopposlist[j]+2]
                        #print "Sequence extraite: ",extractseq
                        if len(extractseq)>=threshold:
                            #translate
                            protein=translate(extractseq,codetable)
                            #print "Traduction: ",protein
                            #store all the useful data in lists
                            finalstartpos.append(startposlist[k])
                            finalstoppos.append(stopposlist[j])
                            finalframe.append(startframelist[k])
                            finallength.append(len(extractseq))
                            finaltranslation.append(protein)
                            cpt=cpt+1

                            print "i", position, "pos", stopposlist[position]
                            print "j", j, "pos", stopposlist[j]
                            print "ORF +1", "frame", stopframelist[position], "start", startposlist[k], "stop", stopposlist[j], "length", len(extractseq)
                            #print cpt


                            #print orflist
                            return finalstartpos, finalstoppos, finalframe, finallength, finaltranslation
                        else:
                            print "trop petit"
                            return finalstartpos, finalstoppos, finalframe, finallength, finaltranslation

            print "Pas de Start"
            return finalstartpos, finalstoppos, finalframe, finallength, finaltranslation

def findORF (seq,threshold,codetable,orflist,sens):
    """
    This function return a list of ORFs in the form of a dictionary.
    Author : Thomas Blanc and Marc Mongy
    Arg : seq - the gene sequence
          threshold - the minimum ORF lsngth expressed in base pairs.
          codeTable - the genetic code associated with the species.
    Returns : A dictionary containing a list of information for each ORFs,
    namely their start and stop positions, their reading frames, their length ad the translated protein sequence.
    """


#Research of all the start and stop codons in the 3 frames.
    #print seq
    #StartTable={}
    #StopTable={}
    startframelist=[]
    startposlist=[]
    stopframelist=[]
    stopposlist=[]
    #print len(seq)
    for i in range(len(seq)):
        print i
        strt=isCodonStart(seq,i,codetable)
        if strt==True:
            frame=(i%3)+1
            if sens == 1:
                frame=-frame

            startframelist.append(frame)
            startposlist.append(i)

            #StartTable[i]=frame
            print "Cadre: ",frame," , "," Position codon start: ",i

        stp=isCodonStop(seq,i,codetable)
        if stp==True:
            frame=(i%3)+1
            if sens == 1:
                frame=-frame
            pos=i

            stopframelist.append(frame)
            stopposlist.append(i)

            #StopTable[i]=frame
            print "Cadre: ",frame," , "," Position codon stop: ",i
    print "CODONS OK"





    #print stopposlist
    ORFs={}
    finalstartpos=[]
    finalstoppos=[]
    finalframe=[]
    finallength=[]
    finaltranslation=[]
    cpt=0

    for i in range(len(stopposlist)):
        #print i
        getORF(seq,threshold,codetable,startposlist,startframelist,stopposlist,stopframelist,i,finalstartpos,finalstoppos,finalframe,finallength,finaltranslation)

    print "LISTE BEGIN"
    orflist.append(finalstartpos)
    print "LISTE 1"
    orflist.append(finalstoppos)
    print "LISTE 2"
    orflist.append(finalframe)
    print "LISTE 3"
    orflist.append(finallength)
    print "LISTE 4"
    orflist.append(finaltranslation)
    print "LISTE OK"

    return orflist

def ORFtableToDict(orflist):
    ORFs={}
    startposition=orflist[0]
    startposition.extend(orflist[5])

    stopposition=orflist[1]
    stopposition.extend(orflist[6])

    cadre=orflist[2]
    cadre.extend(orflist[7])

    longueur=orflist[3]
    longueur.extend(orflist[8])

    traduction=orflist[4]
    traduction.extend(orflist[9])

    """
    for n in range(len(startposition)):
        ORFs[n]=[startposition[n], stopposition[n], cadre[n], longueur[n], traduction[n]]
    """
    listORFs=[]
    for n in range(len(startposition)):
        dictORF={'start':startposition[n], 'stop':stopposition[n], 'frame':cadre[n], 'length':longueur[n], 'sequence_prot':traduction[n]}
        listORFs.append(dictORF)
        print "DICO"
    return listORFs



def getLengths(listORFs):
    """Short description : From an ORF list return a list of ORF lengths

    this function is written by Michele RAFFAELLI .

    Args:
        ORFs - dictionary retuned by the findORF() function.

    Returns:
        The function getLengths returns a list of ORF lengths.
    """

    orf_list=[]
    for i in range(len(listORFs)):
        orf_list.append(listORFs[i]['length'])

    return orf_list

def getLongestORF(orf_list):
    """
    Short description : From an list of  ORF length return the longest ORF
    this function is written by Michele RAFFAELLI .
    Args:
       orf_list: list of orf lengths returned by the function getLengths()
    Returns:
        The function getLongestORF returns the longest ORF
   """
    maxi=0
    for i in range (len(orf_list)):
        if(orf_list[i]>maxi):
            maxi=orf_list[i]
    return maxi



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



def compare(listORFs1,listORFs2):
    """Short description

        is written by Ludwig DUVAL

    Args:
        orfliste1:
        orfliste2:
    Returns:
        la fonction compare retourne une liste contenant les orf produisant les meme
        proteine
    """
    listgeneidentique = []

    for i in range(len(listORFs1)):
        for j in range(len(listORFs_2)):
            if ORFs1[i]['sequence_prot'] == orf2[j]['sequence_prot'] :
                listgeneidentique.append([i, ORFs1[i]['sequence_prot'], j, orf2[j]['sequence_prot']])
    return listgeneidentique



#Begin

orflist=[]

rawFASTA=loadFASTA("my_genome.fasta")
seq=readFASTA(rawFASTA)



print len(seq)

#seq='CTGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'

invert_seq=brinAntiSens(seq)

CodeTable=getGeneticCode(4)

threshold=1000

findORF(seq, threshold, CodeTable, orflist, 0)

findORF(invert_seq, threshold, CodeTable, orflist, 1)



ORFs_FINAL_List=ORFtableToDict(orflist)

for i in range (len(ORFs_FINAL_List)):
    print ORFs_FINAL_List[i]



print "Nombre d'ORFs", len(ORFs_FINAL_List)

orf_length=getLengths(ORFs_FINAL_List)
#print orf_length

print "ORF le plus long", getLongestORF(orf_length)



print "Les longueurs des ", (len(ORFs_FINAL_List)*0.1), "ORFs les plus longs :", getTopLongestORF(orf_length,0.1)
