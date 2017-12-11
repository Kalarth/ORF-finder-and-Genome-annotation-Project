#-*- coding: utf-8 -*-
import os, sys


def readFlatFile(filename):
    """This function is used to extract the text from a Flat File containing a GenBank entry.

    This function is written by MONGY MARC.
    Args: filename: path of the Flat File you want to read.
    Returns: The function return the content of the file in the form of a string
    """    

    fichier=open(filename,'r')
    rawFlatFile=fichier.read()
    fichier.close()
    return rawFlatFile


def getFeatures(txt):
    """This function is used to extract the FEATURES entry from a GenBank entry.
    This function is rewritten by MONGY MARC.
    Args: txt: text from which you want to extract the FEATURES table in string format.
    Returns: The function return the FEATURES table in string format.

    """    
    start_index=txt.find('FEATURES')
    end_index=txt.find('ORIGIN')
    features=txt[start_index:end_index]
    return features



def findEntry(txt, entry_name, entry_end):
    """
    This function is used to extract an entry from a text.
    This function is written by MARC MONGY.
    Args:txt: the text from which you want to extract an entry in string format.
        entry_name: the start of the entry.
        entry_end: the end of the entry.
    Returns: The function return the content of the file in the form of a string.
    """
    valueEntryStart=txt.index(entry_name)
    #print "valueEntryStart ",valueEntryStart
    valueEntryEnd=txt.index(entry_end)
    #print "valueEntryEnd ", valueEntryEnd
    #print "valeur: ", txt.index(entry_end,valueEntryStart)
    entry_line=txt[valueEntryStart:txt.index(entry_end,valueEntryStart)] #En fait, renvoie la plus forte valeur de l'index de entry_end avant le prochain valueEntryStart. J'ai dû m'accrocher avant de comprendre la syntaxe de ce truc.
    tabledEntryLine=entry_line.split(entry_name)
    selectedEntryOnly=tabledEntryLine[1].lstrip(' ').rstrip(' ') #les lstrip et rstrip retirent les espaces qui traînent
    return selectedEntryOnly



def getGenes(txt):
    """
    This function is used to extract the information concerning the coding sequences from the FEATURES table.
    This function is rewritten by MONGY MARC.
    Args: txt: text containing the gene and CDS subsections in string format.
    Returns: The function a list of dictionaries. Each dictionary contain a start, stop, frame, length, name, protein and product key with an associated value.
    """    
    nb_genes=txt.count('     gene            ')
    features_split=txt.split('     gene            ')
    #print features_split

    list_split=[]
    for i in range(0,len(features_split)):
        if 'CDS' in features_split[i]: 
            list_split.append(features_split[i])
        if 'rRNA' in features_split[i]:
            list_split.append(features_split[i])
        if 'tRNA' in features_split[i]:        
            list_split.append(features_split[i])
    #print list_split[i]
    list_genes=[]
    for i in range(0,len(list_split)):
        geneLine='     gene            '+list_split[i]
        #print geneLine
        dict_gene={'start':'','stop':'','frame':'','length':'','name':'','protein':'','product':''}

        if '/gene=' in geneLine:
            dict_gene['name']=findEntry(geneLine,entry_name='/gene="',entry_end='"\n')
        else:
            dict_gene['name']='unknown'

        if  '/translation=' in geneLine:
            dict_gene['protein']=findEntry(geneLine,entry_name='/translation="',entry_end='"\n').replace(' ','').replace('\n','')
        else:
            dict_gene['protein']='xxx'

        if  '/product=' in geneLine:            
            dict_gene['product']=findEntry(geneLine,entry_name='/product="',entry_end='"\n')
        else:
            dict_gene['product']='unknown'

        geneEntryStart='     gene            '
        geneEntryEnd='\n'
        startstop=findEntry(geneLine, geneEntryStart, geneEntryEnd)
        #print startstop

        comp=False
        if 'join' in startstop: #on trouve qqfois des "join" dans les entrées de genbank
            startstop=startstop.replace('join(','').replace(')','').replace('<','').replace('>','').split(',')
            startstop1=startstop[0].split('..')
            #print startstop1
            startstop2=startstop[1].split('..')
            #print startstop2

            dict_gene['start']=int(startstop1[0])
            dict_gene['stop']=int(startstop2[1])
            gene_frame=dict_gene['start']%3
            if gene_frame==0:
                gene_frame=3
            dict_gene['frame']=gene_frame   
            dict_gene['length']=int(startstop1[1])-dict_gene['start']+dict_gene['stop']+2
        else:
            if 'complement' in startstop:
                comp=True
            startstop=startstop.replace('complement(','').replace(')','').replace('<','').replace('>','').split('..') #de cette façon on obtient une liste de positions
            dict_gene['start']=int(startstop[0])
            dict_gene['stop']=int(startstop[1])
            gene_frame=dict_gene['start']%3
            if gene_frame==0:
                gene_frame=3
            if comp==True:
                if gene_frame==1:
                    gene_frame=-1
                if gene_frame==2:
                    gene_frame=-2
                if gene_frame==3:
                    gene_frame=-3
            dict_gene['frame']=gene_frame
            dict_gene['length']=dict_gene['stop']-dict_gene['start']+1
        list_genes.append(dict_gene)
    #list_split.remove(list_split[0])
    return list_genes


def readGenBank(filename):
    """This function is used to extract all the general and genes informations from a GenBank entry.

    This function is written by MARC MONGY.

    Args:
        filename: path of the Flat File you want to read..

    Returns:
        The function return a dictionary containing a description, type, data, ID, length, gbtype, organism, codeTableID and genes key with associated value.
        The genes key value is a list of dictionaries, each containing a start, stop, frame, length, name, protein and product key with an associated value.

    """   
    txt=readFlatFile(filename)
    dict_seq={'description':'','type':'','data':'','ID':'','length':'','gbtype':'','organism':'','codeTableID':'','genes':''}

    descriptionEntry='DEFINITION'
    descriptionEnd='\n'
    description=findEntry(txt, descriptionEntry, descriptionEnd)#.replace(' ','').replace('\n','')
    if description != '':
        dict_seq['description']=findEntry(txt, descriptionEntry, descriptionEnd)
    else:
        dict_seq['description']='No Definition'

    originEntry='ORIGIN'
    originEnd='//'
    data=findEntry(txt, originEntry ,originEnd).replace(' ','').replace('\n','')
    if data != '':
        dict_seq['data']=findEntry(txt, originEntry, originEnd)
    else:
        dict_seq['data']='No Data'

    locusEntry='LOCUS'
    locusEnd='\n'
    locus=findEntry(txt, locusEntry, locusEnd).split('  ')
    while '' in locus:
        locus.remove('')
    dict_seq['ID']=locus[0]
    dict_seq['length']=locus[1].strip(' ')
    dict_seq['gbtype']=locus[2].strip(' ')
    if dict_seq['gbtype']=='DNA':
        dict_seq['type']='DNA'
    elif dict_seq['gbtype']=='RNA':
        dict_seq['type']='RNA'
    else:
        dict_seq['type']='protein'

    organismEntry='ORGANISM'
    organismEnd='\n'
    dict_seq['organism']=findEntry(txt, organismEntry, organismEnd)
    if '/transl_table=' in txt:
        translationTableEntry='/transl_table='
        translationTableEnd='\n'
        dict_seq['codeTableID']=findEntry(txt, translationTableEntry, translationTableEnd)    
    features=getFeatures(txt)
    dict_seq['genes']=getGenes(features)
    return dict_seq



dictionary=readGenBank("genbank_entry.gb")
print dictionary['description']
print dictionary['type']
print dictionary['ID']
print dictionary['length']
print dictionary['gbtype']
print dictionary['organism']
print dictionary['codeTableID']
print dictionary['data']
genes=dictionary['genes']
for i in range(len(genes)):
    print "Gène ",i, ": ",genes[i]

