#-*- coding: utf-8 -*-
import os, sys
#import myBio_2017 as bio


'''---------------Partie 0 - Charger le FASTA-------------'''

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
    '''Return a dictionary separating the header from the data'''
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
    """
    #Tout ce merdier est un reliquat de l'année dernière qui identifiait automatiquement le type de séquence. Ici, il gène plus qu'autre chose.
    if moltype=='auto':
        if isRNA(seq):
            moltype='rna'
        elif isDNA(seq):
            moltype='dna'
        else:
            moltype='protein'

    seq['type']=moltype
    """
    return seq


def readFASTA(txt): #récupère la séquence en String dans le dico
    sequence=dictionaryFromFASTA(txt, moltype='auto')
    seqString=sequence['data']
    return seqString


'''---------------Partie 1 - ORF-------------'''

def getGeneticCode(codeTable):
    """Returns a dictionary with the coding table corresponding to the number"""
    table={}

    if codeTable==4: #code des Mycoplasmes
        base1="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
        base2="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
        base3="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
        AcAms="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
        Start="--MM------**-------M------------MMMM---------------M------------"

    for i in range(0,len(AcAms)):
        codon=base1[i]+base2[i]+base3[i]
        matchingAA=AcAms[i]
        starter=Start[i]
        aa=[]
        aa.append(matchingAA)
        aa.append(starter)
        table[codon]=aa
    #print table
    return table



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


def getCodonsStart(codeTable):
    codonsStart=[]
    table=getGeneticCode(codeTable)
    """ CE TRUC NE SERT A RIEN
    for n in table.keys():
        acideAmine=table[n]
    matchingAA=acideAmine[0]
    starter=acideAmine[1]
"""
    for n in table.keys():
        if table[n][1]=="M":
            codonsStart.append(n)
    #print codonsStart
    return codonsStart


def getCodonsStop(codeTable):
    codonsStop=[]
    table=getGeneticCode(codeTable)
    """ CE TRUC NE SERT A RIEN
    for n in table.keys():
        acideAmine=table[n]
    matchingAA=acideAmine[0]
    starter=acideAmine[1]
"""
    for n in table.keys():
        if table[n][1]=="*":
            codonsStop.append(n)
    #print codonsStop
    return codonsStop

""" INUTILE, MA FONCTION FAIS LA MEME EN MIEUX

def oneWord(seq, start, wlen):
    string=""
    fin = start + wlen
    for i in range(len(seq)):
        if i>=start and i<fin:
            #for start in range(start,start+wlen,1):
            string=string+seq[i]
    return string

"""
def isCodonStart(seq, pos, codeTable):
    codonsStart=getCodonsStart(codeTable)
    flag=False
    codon=oneWord(seq,pos,3)
    #print codon
    for j in range(len(codonsStart)):
        if codon==codonsStart[j]:
            flag=True
    return flag


def isCodonStop(seq, pos, codeTable):
    codonsStop=getCodonsStop(codeTable)
    flag=False
    codon=oneWord(seq,pos,3)
    #print codon
    for j in range(len(codonsStop)):
        if codon==codonsStop[j]:
            flag=True
    return flag


def parcoursORF(seq, threshold, codeTable, cadreDeLecture):

    complementaireSeq=brinAntiSens(seq)

    if cadreDeLecture=="1":
        sequence=seq
    elif cadreDeLecture=="2":
        sequence=seq[1:]
    elif cadreDeLecture=="3":
        sequence=seq[2:]
    elif cadreDeLecture=="reverse1":
        sequence=complementaireSeq
    elif cadreDeLecture=="reverse2":
        sequence=complementaireSeq[1:]
    elif cadreDeLecture=="reverse3":
        sequence=complementaireSeq[2:]

    cadrePositionsCodonsStart=[]
    cadrePositionsCodonsStop=[]
    for positionCodonStart in range(0, len(sequence), 3):
        if isCodonStart(sequence, positionCodonStart, codeTable)==True:
            cadrePositionsCodonsStart.append(positionCodonStart)
    for positionCodonStop in range(0, len(sequence), 3):
        if isCodonStop(sequence, positionCodonStop, codeTable)==True:
            cadrePositionsCodonsStop.append(positionCodonStop)
    #print "Position des codons Start pour l'ORF", cadrePositionsCodonsStart
    #print "Position des codons Stop pour l'ORF", cadrePositionsCodonsStop

    ORF=[]
    numeroSeqORF=[]
    positionCodonStartORF=[]
    positionCodonStopORF=[]
    tailleSeqORF=[]
    sequenceORF=[]

    stringORF=""
    numeroCodonStartCadre=0
    numeroCodonStopCadre=0
    numero=0
    while numeroCodonStopCadre < len(cadrePositionsCodonsStop):
        while numeroCodonStartCadre < len(cadrePositionsCodonsStart):
            if cadrePositionsCodonsStart[numeroCodonStartCadre] < cadrePositionsCodonsStop[numeroCodonStopCadre]:
                #print cadrePositionsCodonsStart[numeroCodonStartCadre]
                #print cadrePositionsCodonsStop[numeroCodonStopCadre]+3
                tailleSequence=(cadrePositionsCodonsStop[numeroCodonStopCadre]+3)-(cadrePositionsCodonsStart[numeroCodonStartCadre])
                if tailleSequence >= threshold:
                    n=cadrePositionsCodonsStart[numeroCodonStartCadre]
                    while n < cadrePositionsCodonsStop[numeroCodonStopCadre]+3:
                        stringORF=stringORF+sequence[n]
                        n=n+1
                    numero=numero+1
                    ORF.append(cadreDeLecture)
                    numeroSeqORF.append(numero)
                    positionCodonStartORF.append(cadrePositionsCodonsStart[numeroCodonStartCadre])
                    positionCodonStopORF.append(cadrePositionsCodonsStop[numeroCodonStopCadre])
                    tailleSeqORF.append(tailleSequence)
                    sequenceORF.append(stringORF)
            numeroCodonStartCadre=numeroCodonStartCadre+1
        numeroCodonStopCadre=numeroCodonStopCadre+1

    tableauGlobalORF=[]
    tableauGlobalORF.append(ORF)
    tableauGlobalORF.append(numeroSeqORF)
    tableauGlobalORF.append(positionCodonStartORF)
    tableauGlobalORF.append(positionCodonStopORF)
    tableauGlobalORF.append(tailleSeqORF)
    tableauGlobalORF.append(sequenceORF)

    return tableauGlobalORF


def ORFtableToDict(tableauGlobalORF):
    ORF=tableauGlobalORF[0]
    numeroSeqORF=tableauGlobalORF[1]
    positionCodonStartORF=tableauGlobalORF[2]
    positionCodonStopORF=tableauGlobalORF[3]
    tailleSeqORF=tableauGlobalORF[4]
    sequenceORF=tableauGlobalORF[5]
    listORFs=[]
    for n in range(len(sequenceORF)):
        dictORF={'numero':numeroSeqORF[n] , 'start':positionCodonStartORF[n], 'stop':positionCodonStopORF[n], 'frame':ORF[n], 'length':tailleSeqORF[n], 'sequence_ADN':sequenceORF[n]}
        listORFs.append(dictORF)
    return listORFs


def findORF(seq,threshold,codeTable):

    #-----------------------------------
    #-----Cadre de lecture 1 (ORF1)-----
    #-----------------------------------
    print "Finding genes in ORF1"
    tableauORF1=parcoursORF(seq, threshold, codeTable, "1")
    listORFs_1=ORFtableToDict(tableauORF1)
    print listORFs_1
    #-----------------------------------
    #-----Cadre de lecture 2 (ORF2)-----
    #-----------------------------------
    print "Finding genes in ORF2"
    tableauORF2=parcoursORF(seq, threshold, codeTable, "2")
    listORFs_2=ORFtableToDict(tableauORF2)
    print listORFs_2
    #-----------------------------------
    #-----Cadre de lecture 3 (ORF3)-----
    #-----------------------------------
    print "Finding genes in ORF3"
    tableauORF3=parcoursORF(seq, threshold, codeTable, "3")
    listORFs_3=ORFtableToDict(tableauORF3)
    print listORFs_3
    #--------------------------------------------
    #-----Cadre de lecture Reverse 1 (ORF-1)-----
    #--------------------------------------------
    print "Finding genes in ORF reverse1"
    tableauORFReverse1=parcoursORF(seq, threshold, codeTable, "reverse1")
    listORFs_reverse1=ORFtableToDict(tableauORFReverse1)
    print listORFs_reverse1
    #--------------------------------------------
    #-----Cadre de lecture Reverse 2 (ORF-2)-----
    #--------------------------------------------
    print "Finding genes in ORF reverse2"
    tableauORFReverse2=parcoursORF(seq, threshold, codeTable, "reverse2")
    listORFs_reverse2=ORFtableToDict(tableauORFReverse2)
    print listORFs_reverse2
    #--------------------------------------------
    #-----Cadre de lecture Reverse 3 (ORF-3)-----
    #--------------------------------------------
    print "Finding genes in ORF reverse3"
    tableauORFReverse3=parcoursORF(seq, threshold, codeTable, "reverse3")
    listORFs_reverse3=ORFtableToDict(tableauORFReverse3)
    print listORFs_reverse3


'''---------------Partie 2 - Stats and utilitary-------------'''


def writeCSV(filename, separator, data):
    '''Marc MONGY'''
    '''Cette fonction utilise une série de fonctions de manipulation de listes, de chaînes de caractères et de dictionnaires pour convertir l'ORF (sous forme de dictionnaire Python) et l'enregistrer sous un format de fichier CSV. Le paramètre "separator" représente le séparateur utilisé (ici, le point-virgule)'''
    keys=list(data.keys())
    values=list(data.values())
    result=readCSV(filename, separator)
    with open ("my_data.csv", "w") as my_data:
        separator=str(separator)
        keys=str(keys)
        keys=keys.replace(",", separator)
        keys=keys.replace("'", "")
        keys=keys.replace(" ","")
        keys=keys.replace("[","")
        keys=keys.replace("]","")
        values=str(values)
        values=values.replace(',', separator)
        values=values.replace(" ","")
        values=values.replace("[","")
        values=values.replace("]","")
        error=my_data.write(keys)
        error=my_data.write(str('\n'))
        error=my_data.write(values)
        my_data.close()
    return error



def readCSV(filename, separator):
    '''Marc MONGY'''
    '''Cette fonction effectue la conversion du format CSV vers un dictionnaire (opération inverse de writeCSV) afin de permettre l'utilisation de l'ORF sous forme de dictionnaire Python par d'autres scripts. '''
    with open ("my_data.csv", "r") as my_data:
        result=my_data.read()
        separator=str(separator)
        result=result.replace("'","")
        result=result.replace(separator, "', '")
        result="['" + result + "']"
        result=result.replace("\n", "']\n['")
        result=result.split('\n')
        result[0]=result[0].replace("'", "")
        result[0]=result[0].replace("[", "")
        result[0]=result[0].replace("]", "")
        result[0]=result[0].replace(" ","")
        result[0]=result[0].split(',')
        keys=result[0]
        result[1]=result[1].replace("'", "")
        result[1]=result[1].replace("[", "")
        result[1]=result[1].replace("]", "")
        result[1]=result[1].replace(" ","")
        result[1]=result[1].split(',')
        values=result[1]
        result={key:value for key, value in zip(keys,values)}
        result['id']=int(result['id'])
        result['start']=int(result['start'])
        result['stop']=int(result['stop'])
        result['name']=str(result['name'])
        my_data.close()
    return result




rawFASTA=loadFASTA("my_genome.fasta")
seq='CTGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'
#threshold=3*90
threshold=1
codeTable=4
findORF(seq,threshold,codeTable)
