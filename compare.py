


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
