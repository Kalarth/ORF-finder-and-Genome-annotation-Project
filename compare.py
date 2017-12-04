


def compare(orfliste1,orfliste2):
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
    listProteine1 = translate(orfliste1)
    listProteine2 = translate(orfliste2)
    for i in range(len(listProteine1)):
        for j in range(len(listProteine2)):
            orf1 = oneWord(listProteine1[i],0,len(listProteine1[i]))
            orf2 = oneWord(listProteine2[j],0,len(listProteine2[j]))
            if orf1 == orf2 :
                listgeneidentique.append([i,listProteine1[i],j,listProteine2[j]])
    return listgeneidentique
