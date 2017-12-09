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


def oneWord(seq, start, wlen):
    string=""
    fin = start + wlen
    for i in range(len(seq)):
        if i>=start and i<fin:
            #for start in range(start,start+wlen,1):
            string=string+seq[i]
    return string


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


def defineORF(sequence, numeroCadreORF, listePositionsCodonsStart, listePositionsCodonsStop):
    listORFs=[]
    #stringORF=""
    numeroCodonStartCadre=0
    numeroCodonStopCadre=0
    numero=0
    while numeroCodonStopCadre < len(listePositionsCodonsStop):
        while numeroCodonStartCadre < len(listePositionsCodonsStart):
            if listePositionsCodonsStart[numeroCodonStartCadre] < listePositionsCodonsStop[numeroCodonStopCadre]:
                #print cadrePositionsCodonsStart[numeroCodonStartCadre]
                #print cadrePositionsCodonsStop[numeroCodonStopCadre]+3
                tailleSequence=(listePositionsCodonsStop[numeroCodonStopCadre]+3)-(listePositionsCodonsStart[numeroCodonStartCadre])
                if tailleSequence >= threshold:
                    stringORF=""
                    n=listePositionsCodonsStart[numeroCodonStartCadre]
                    while n < listePositionsCodonsStop[numeroCodonStopCadre]+3:
                        stringORF=stringORF+sequence[n]
                        n=n+1
                    print stringORF
                    dictORF={'numero':numero, 'start':listePositionsCodonsStart[numeroCodonStartCadre], 'stop':listePositionsCodonsStop[numeroCodonStopCadre], 'frame':numeroCadreORF, 'length':tailleSequence, 'sequence_ADN':stringORF}
                    listORFs.append(dictORF)
                    print "Appending sequence ",numero," , ORF ",numeroCadreORF
                    numero=numero+1
            numeroCodonStartCadre=numeroCodonStartCadre+1
        numeroCodonStopCadre=numeroCodonStopCadre+1


    return listORFs


def parcoursORF(seq, threshold, codeTable):
    sequence=seq
    complementaireSeq=brinAntiSens(seq)

    cadreORF1PositionsCodonsStart=[]
    cadreORF1PositionsCodonsStop=[]
    cadreORF2PositionsCodonsStart=[]
    cadreORF2PositionsCodonsStop=[]
    cadreORF3PositionsCodonsStart=[]
    cadreORF3PositionsCodonsStop=[]
    cadreORFreverse1PositionsCodonsStart=[]
    cadreORFreverse1PositionsCodonsStop=[]
    cadreORFreverse2PositionsCodonsStart=[]
    cadreORFreverse2PositionsCodonsStop=[]
    cadreORFreverse3PositionsCodonsStart=[]
    cadreORFreverse3PositionsCodonsStop=[]

    for positionCodon in range(0, len(sequence)):
        if positionCodon%3==0:
            if isCodonStart(sequence, positionCodon, codeTable)==True:
                cadreORF1PositionsCodonsStart.append(positionCodon)
                print "ORF 1, Codon Start, position ", positionCodon
            elif isCodonStop(sequence, positionCodon, codeTable)==True:
                cadreORF1PositionsCodonsStop.append(positionCodon)
                print "ORF 1, Codon Stop, position ", positionCodon
        elif positionCodon%3==1:
            if isCodonStart(sequence, positionCodon, codeTable)==True:
                cadreORF2PositionsCodonsStart.append(positionCodon)
                print "ORF 2, Codon Start, position ", positionCodon
            elif isCodonStop(sequence, positionCodon, codeTable)==True:
                cadreORF2PositionsCodonsStop.append(positionCodon)
                print "ORF 2, Codon Stop, position ", positionCodon
        elif positionCodon%3==2:
            if isCodonStart(sequence, positionCodon, codeTable)==True:
                cadreORF3PositionsCodonsStart.append(positionCodon)
                print "ORF 3, Codon Start, position ", positionCodon
            elif isCodonStop(sequence, positionCodon, codeTable)==True:
                cadreORF3PositionsCodonsStop.append(positionCodon)
                print "ORF 3, Codon Stop, position ", positionCodon

    for positionCodon in range(0, len(complementaireSeq)):
        if positionCodon%3==0:
            if isCodonStart(complementaireSeq, positionCodon, codeTable)==True:
                cadreORFreverse1PositionsCodonsStart.append(positionCodon)
                print "ORF -1, Codon Start, position ", positionCodon
            elif isCodonStop(complementaireSeq, positionCodon, codeTable)==True:
                cadreORFreverse1PositionsCodonsStop.append(positionCodon)
                print "ORF -1, Codon Stop, position ", positionCodon
        elif positionCodon%3==1:
            if isCodonStart(complementaireSeq, positionCodon, codeTable)==True:
                cadreORFreverse2PositionsCodonsStart.append(positionCodon)
                print "ORF -2, Codon Start, position ", positionCodon
            elif isCodonStop(complementaireSeq, positionCodon, codeTable)==True:
                cadreORFreverse2PositionsCodonsStop.append(positionCodon)
                print "ORF -2, Codon Stop, position ", positionCodon
        elif positionCodon%3==2:
            if isCodonStart(complementaireSeq, positionCodon, codeTable)==True:
                cadreORFreverse3PositionsCodonsStart.append(positionCodon)
                print "ORF -3, Codon Start, position ", positionCodon
            elif isCodonStop(complementaireSeq, positionCodon, codeTable)==True:
                cadreORFreverse3PositionsCodonsStop.append(positionCodon)
                print "ORF -3, Codon Stop, position ", positionCodon
    #print "Position des codons Start pour l'ORF", cadrePositionsCodonsStart
    #print "Position des codons Stop pour l'ORF", cadrePositionsCodonsStop

    listORFs_1=defineORF(sequence, "1", cadreORF1PositionsCodonsStart, cadreORF1PositionsCodonsStop)
    listORFs_2=defineORF(sequence, "2", cadreORF2PositionsCodonsStart, cadreORF2PositionsCodonsStop)
    listORFs_3=defineORF(sequence, "3", cadreORF3PositionsCodonsStart, cadreORF3PositionsCodonsStop)
    listORFs_reverse1=defineORF(complementaireSeq, "-1", cadreORFreverse1PositionsCodonsStart, cadreORFreverse1PositionsCodonsStop)
    listORFs_reverse2=defineORF(complementaireSeq, "-2", cadreORFreverse2PositionsCodonsStart, cadreORFreverse2PositionsCodonsStop)
    listORFs_reverse3=defineORF(complementaireSeq, "-3", cadreORFreverse3PositionsCodonsStart, cadreORFreverse3PositionsCodonsStop)

    fullListORF=listORFs_1+listORFs_2+listORFs_3+listORFs_reverse1+listORFs_reverse2+listORFs_reverse3

    return fullListORF


def findORF(seq,threshold,codeTable):
    fullListORF=parcoursORF(seq, threshold, codeTable)
    for i in range(len(fullListORF)):
        print fullListORF[i]
    print len(fullListORF)


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
#seq='CTGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'
seq=readFASTA(rawFASTA)

#seq=seq2[0:50000]
#threshold=3*90
#threshold=1
threshold=600
codeTable=4
findORF(seq,threshold,codeTable)
