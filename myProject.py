
import myBio as bio

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

    return seq


def readFASTA(txt): #Get the DNA sequence from the dictionary into a string
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
    flag=bio.isDNA(seq)
    if flag==True:
        for i in range(len(seq)):
            codon=bio.oneWord(seq,i,3)
            for tablecodon in codetable.keys():
                if codon==tablecodon:
                    seqprot=seqprot+codetable[tablecodon][0]
        return seqprot
    else:
        return "error, sequence is not dna."

def findORF (seq,threshold,codetable,orflist,sens):
    """
    This function return a list of ORFs in the form of a dictionary.
    Author : Thomas Blanc
    Arg : seq - the gene sequence
          threshold - the minimum ORF lsngth expressed in base pairs.
          codeTable - the genetic code associated with the species.
    Returns : A dictionary containing a list of information for each ORFs,
    namely their start and stop positions, their reading frames, their length ad the translated protein sequence.
    """


#Research of all the start and stop codons in the 3 frames.
    #print seq
    startframelist=[]
    startposlist=[]
    stopframelist=[]
    stopposlist=[]
    #print len(seq)
    for i in range(len(seq)):
        strt=bio.isCodonStart(seq,i,codetable)
        if strt==True:
            frame=(i%3)+1
            if sens == 1:
                frame=-frame
            startframelist.append(frame)
            startposlist.append(i)
            print "Cadre: ",frame," , "," Position codon start: ",i

        stp=bio.isCodonStop(seq,i,codetable)
        if stp==True:
            frame=(i%3)+1
            if sens == 1:
                frame=-frame
            pos=i
            stopframelist.append(frame)
            stopposlist.append(i)
            print "Cadre: ",frame," , "," Position codon stop: ",i
    print "CODONS OK"
    #dico[cle][3]




    print stopposlist
    ORFs={}
    finalstartpos=[]
    finalstoppos=[]
    finalframe=[]
    finallength=[]
    finaltranslation=[]
    cpt=0
    i=0
    x=0
    """
    for i in range (len(stopposlist)):

        for j in range (i+1,len(stopposlist)):
            """
    while i <= (len(stopposlist)-1):
        j=i+1
        print i
        print j
        while x == 0:
            if j > len(stopframelist):
                x=1

            if stopframelist[i]==stopframelist[j]:


                for k in range (len(startposlist)):

                 if startposlist[k]>stopposlist[i] and startposlist[k]<stopposlist[j] and startframelist[k]==stopframelist[i]:
                        #extract sequence
                        extractseq=seq[startposlist[k]:stopposlist[j]]
                        print "Sequence extraite: ",extractseq

                        if len(extractseq)>threshold:
                            #translate
                            protein=translate(extractseq,codetable)
                            print "Traduction: ",protein


                            #store all the useful data in lists
                            finalstartpos.append(startposlist[k])
                            finalstoppos.append(stopposlist[j])
                            finalframe.append(stopframelist[i])
                            finallength.append(len(extractseq))
                            finaltranslation.append(protein)
                            cpt=cpt+1

                            print cpt
                            x=1

            j=j+1
        i=i+1

    orflist.append(finalstartpos)
    orflist.append(finalstoppos)
    orflist.append(finalframe)
    orflist.append(finallength)
    orflist.append(finaltranslation)
    #print orflist
    return orflist

def ORFtableToDict(tableauGlobalORF):
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
        orf_list.extend(ORFs[i]['length'])

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



def writeCSV(filename, separator, data):
    '''Marc MONGY'''
    '''Cette fonction utilise une serie de fonctions de manipulation de listes,
    de chaines de caracteres et de dictionnaires pour convertir l'ORF (sous forme de dictionnaire Python)
    et l'enregistrer sous un format de fichier CSV.
    Le parametre "separator" represente le separateur utilise (ici, le point-virgule)
    '''
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
    """Marc MONGY"""
    """Cette fonction effectue la conversion du format CSV vers un dictionnaire (operation inverse de writeCSV)
    afin de permettre l'utilisation de l'ORF sous forme de dictionnaire Python par d'autres scripts.
    """
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

#Begin
orflist=[]
"""
rawFASTA=loadFASTA("my_genome.fasta")
seq=readFASTA(rawFASTA)
"""
seq='CTGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'

invert_seq=bio.brinAntiSens(seq)

CodeTable=bio.getGeneticCode(4)

threshold=0

findORF(seq, threshold, CodeTable, orflist, 0)

findORF(invert_seq, threshold, CodeTable, orflist, 1)



ORFs_FINAL_List=ORFtableToDict(orflist)



print len(ORFs_FINAL_List)
