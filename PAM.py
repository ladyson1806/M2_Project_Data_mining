"""
def distance_map():
    distance_map = {}
    list_prot = ["a", "b", "c", "d", "e"]
    for prot1 in list_prot:
        distance_map_prot = {}
        for prot2 in list_prot:
            dist = ???
            distance_map_prot[prot2] = dist
        distance_map[prot1] = distance_map_prot
"""
        
def PAM(liste_prot, distance_map, k_medoid):
    for prot in liste_prot:
        k_medoid = PAM2(prot, distance_map, k_medoid)
    i = 0
    while i < len(k_medoid):
        result = swap(k_medoid[i][0], k_medoid[i])
        if result[0] == False:
            k_medoid = k_medoid[:i]
            for cluster_definitif in k_medoid:
                for prot in cluster_definitif:
                    liste_prot.remove(prot)
            k_medoid += PAM(liste_prot, distance_map, [result[1]])
        else:
            i += 1
    return k_medoid

def PAM2(prot, distance_map,k_medoid):
        for i in range (len(k_medoid)):
            if prot not in k_medoid[i]:
                if distance_map[k_medoid[i][0]][prot] <= 5.0:
                    k_medoid[i].append(prot)
                    return k_medoid
            else:
                return k_medoid
        new_cluster = [prot]
        k_medoid.append(new_cluster)
        return k_medoid

def swap(k,prots):
    distance_totale = calcul_distance(k, prots)
    new_cluster = []
    for prot in prots:
        distance_totale_differee = calcul_distance(prot, prots)
        if distance_totale_differee < distance_totale:
            new_cluster = [prot]
            distance_totale = distance_totale_differee
    if new_cluster != []:
        return False, new_cluster
    else:
        return True, prots

def calcul_distance(k, prots):
    distance_totale = 0
    for prot in prots:
        distance_totale += distance_map[k][prot]
    return distance_totale



"""exemple de jeu de donnees"""
distance_map_g1 = {"g1":0.0, "g2":8.1, "g3":9.2, "g4":7.7, "g5":9.3, "g6":2.3, "g7":5.1, "g8":10.2, "g9":6.1, "g10":7.0}
distance_map_g2 = {"g1":8.1, "g2":0.0, "g3":12.0, "g4":0.9, "g5":12.0, "g6":9.5, "g7":10.1, "g8":12.8, "g9":2.0, "g10":1.0}
distance_map_g3 = {"g1":9.2, "g2":12.0, "g3":0.0, "g4":11.2, "g5":0.7, "g6":11.1, "g7":8.1, "g8":1.1, "g9":10.5, "g10":11.5}
distance_map_g4 = {"g1":7.7, "g2":0.9, "g3":11.2, "g4":0.0, "g5":11.2, "g6":9.2, "g7":9.5, "g8":12.0, "g9":1.6, "g10":1.1}
distance_map_g5 = {"g1":9.3, "g2":12.0, "g3":0.7, "g4":11.2, "g5":0.0, "g6":11.2, "g7":8.5, "g8":1.0, "g9":10.6, "g10": 11.6}
distance_map_g6 = {"g1":2.3, "g2":9.5, "g3":11.1, "g4":9.2, "g5":11.2, "g6":0.0, "g7":5.6, "g8":12.1, "g9":7.7, "g10":8.5}
distance_map_g7 = {"g1":5.1, "g2":10.1, "g3":8.1, "g4":9.5, "g5":8.5, "g6":5.6, "g7":0.0, "g8":9.1, "g9":8.3, "g10":9.3}
distance_map_g8 = {"g1":10.2, "g2":12.8, "g3":1.1, "g4":12.0, "g5":1.0, "g6":12.1, "g7":9.1, "g8":0.0, "g9":11.4, "g10":12.4}
distance_map_g9 = {"g1":6.1, "g2":2.0, "g3":10.5, "g4":1.6, "g5":10.6, "g6":7.7, "g7":8.3, "g8":11.4, "g9":0.0, "g10":1.1}
distance_map_g10 = {"g1":7.0, "g2":1.0, "g3":11.5, "g4":1.1, "g5":11.6, "g6":8.5, "g7":9.3, "g8":12.4, "g9":1.1, "g10":0.0}
distance_map = {}
distance_map["g1"] = distance_map_g1
distance_map["g2"] = distance_map_g2
distance_map["g3"] = distance_map_g3
distance_map["g4"] = distance_map_g4
distance_map["g5"] = distance_map_g5
distance_map["g6"] = distance_map_g6
distance_map["g7"] = distance_map_g7
distance_map["g8"] = distance_map_g8
distance_map["g9"] = distance_map_g9
distance_map["g10"] = distance_map_g10
liste_prot = ["g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10"]

k_medoid = PAM(liste_prot, distance_map, [])
print k_medoid
