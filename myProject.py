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


def translate(seq,codetable) :
    """
    This function translate a Dna sequence into a proteic sequence.
    Arg : seq - the gene sequence
          codetable - the genetic code associated with the species.
    Returns : The proteic sequence as a string.
    """
    seqprot=""
    flag=bio.isDNA(seq)
    if flag==True:
        for i in range(len(seq)):
            codon=bio.oneWord(seq,i,3)
            for tablecodon in codetable.keys:
                if codon==tablecodon:
                    seqprot=seqprot+codetable[tablecodon][0]
        return seqprot
    else:
        return "error, sequence is not dna."

def findORF (seq,threshold,codetable):
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
        strt=bio.isCodonStart(seq,i,codetable)
        if strt==True:
            frame=i%3
            pos=i
            startdic[i]=frame

        stp=bio.isCodonStop(seq,i,codetable)
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
                        protein=translate(extractseq,codetable)


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
        if(orf_list[i]>max):
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



def compare(ORFs1,ORFs2):
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

    for i in ORFs1.keys:
        for j in ORFs2.keys:
            if ORFs1[i][4] == orf2[j][4] :
                listgeneidentique.append([i, ORFs1[i][4], j, orf2[j][4])
    return listgeneidentique



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
