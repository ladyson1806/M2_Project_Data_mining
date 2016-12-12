#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import division
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
import numpy as np
import sys
import math
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import sklearn
import pprint

def loadfile(protein_list):
    """
    Load a tab separated file in a list
    """
    try:
        f = open(sys.argv[1], 'r')
    except:
        print ("Le fichier ", sys.argv[1], " est introuvable")
        sys.exit()
        
    allLines = [x.strip() for x in f.readlines()]
    f.close()
    key_list = allLines[0].split("\t") #List of keywords
    for line in allLines:
        protein = {}
        line = line.split("\t")
        for x in range(len(line)):
            if x==range(len(line)):
                protein[key_list[x]] = line[x].rstrip("\n")
            else:
                protein[key_list[x]] = line[x]
        protein_list.append(protein)
    return protein_list[1:len(protein_list)] #Doesn't return the first line with keywords

def formatOntoData(protein_list):
    """
    A function to add Gene ontology if it's missing from data (maybe it's a tab format problem)
    then format every GO for every protein as a list (Gene ontology is None if it's empty).
    It also makes sure that mass and length are in the correct format.
    """
    onto = "Gene ontology (GO)"
    mass = "Mass"
    length = "Length"
    for protein in protein_list:
        if not protein.has_key(onto):
            protein[onto] = None
        if not protein.has_key(mass):
            protein[mass] = None
        if not protein.has_key(length):
            protein[length] = None
        if not protein[onto] is None:
            protein[onto] = protein[onto].split(";")
        if not protein[mass] is None:
            protein["Mass"] = float(protein["Mass"].replace(",", ""))
        if not protein[length] is None:
            protein["Length"] = float(protein["Length"].replace(",", ""))

def geneOntoAsList(protein_list):
    """ 
    For every protein, get the list of a particular GO, and its entry ID.
    """
    onto = "Gene ontology (GO)"
    go_list = []
    entry_list = []
    for protein in protein_list:
        go_list.append(protein[onto])
        entry_list.append(protein["Entry"])
    return go_list, entry_list

def getMinMaxMass(protein_list):
    """
    Return the maximum and minimum mass of proteins in a protein list
    """
    mass = "Mass"
    massMax = 0
    massMin = 1000 #arbitrary choosen
    for protein in protein_list:
        if not protein[mass] is None:
            if protein[mass] > massMax:
                massMax = protein[mass]
            if protein[mass] < massMin:
                massMin = protein[mass]
    return massMin, massMax

def getMinMaxLength(protein_list):
    """
    Return the maximum and minimum length of proteins in a protein list
    """
    length = "Length"
    lengthMax = 0
    lengthMin = 1000 #arbitrary choosen
    for protein in protein_list:
        if not protein[length] is None:
            if protein[length] > lengthMax:
                lengthMax = protein[length]
            if protein[length] < lengthMin:
                lengthMin = protein[length]
    return lengthMin, lengthMax


def calcMassDistance(protein1, protein2, massMin, massMax):
    delta = 0
    dist = 0
    if not protein1["Mass"] is None and not protein2["Mass"] is None:
        delta = 1
        dist = abs(protein1["Mass"] - protein2["Mass"]) / (massMax - massMin)
    return delta, dist

def calcLengthDistance(protein1, protein2, lengthMin, lengthMax):
    delta = 0
    dist = 0
    if not protein1["Length"] is None and not protein2["Length"] is None:
        delta = 1
        dist = abs(protein1["Length"] - protein2["Length"]) / (lengthMax - lengthMin)
    return delta, dist
        

def calcGODistance(protein1, protein2):
    """
    A function to calculate the distance between two proteins considering a particular GO. 
    It calculates the number of shared GO (shared_nodes), and unique GO. 
    Operation:
    1 - (number of shared nodes) / (number of shared nodes + number of unique nodes)
    """
    onto = "Gene ontology (GO)"
    delta = 0
    shared_nodes = 0
    unique_nodes = 0
    if not protein1[onto] is None and not protein2[onto] is None:
        delta = 1
        if len(protein1[onto]) >= len(protein2[onto]):
            for go in protein1[onto]:
                if go in protein2[onto]:
                    shared_nodes+=1
                else:
                    unique_nodes+=1
        else:
            for go in protein2[onto]:
                if go in protein1[onto]:
                    shared_nodes+=1
                else:
                    unique_nodes+=1
        return delta, round(1 - (shared_nodes / (shared_nodes + unique_nodes)), 3)
    else:
        return delta, 0

def calcDistance(protein1, protein2, massMin, massMax, lengthMin, lengthMax):
    deltaMass, massDistance = calcMassDistance(protein1, protein2, massMin, massMax)
    deltaLength, lengthDistance = calcLengthDistance(protein1, protein2, lengthMin, lengthMax)
    deltaGO, GODistance = calcGODistance(protein1, protein2)
    try:
        return round(((deltaMass*massDistance + deltaLength*lengthDistance + deltaGO*GODistance) /
                      (deltaMass + deltaLength + deltaGO)), 3)
    except ZeroDivisionError:
        return 1
    
def generateDistanceMatrix(protein_list):
    """
    A function to calculate a distance matrix between proteins.
    """
    massMin, massMax = getMinMaxMass(protein_list)
    lengthMin, lengthMax = getMinMaxLength(protein_list)
    distance_matrix = []
    i=0
    for protein in protein_list:
        distance_list = []
        for protein_checked in protein_list:
            distance_list.append(calcDistance(protein, protein_checked,
                                              massMin, massMax, lengthMin, lengthMax))
        i+=1
        print i
        distance_matrix.append(distance_list)
    return distance_matrix



if __name__ == '__main__':
    protein_list = []
    protein_list = loadfile(protein_list)
    formatOntoData(protein_list)
    distance_matrix = generateDistanceMatrix(protein_list[0:10000])
    
    """
    go_list_mol, entry_list = geneOntoAsList(protein_list[:5])

    distance_matrix = generateDistanceMatrix(go_list_mol)
    print distance_matrix
    """
            
    """
    distArray = ssd.squareform(distance_matrix)
    Z = linkage(distArray, 'ward')
    #c, coph_dist = cophenet(Z, pdist(
    
    
    db = DBSCAN(eps=0.5, min_samples=2, metric="precomputed").fit(distance_matrix)
    
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)


    cluster = {}
    for nb in labels:
        if not nb in cluster.keys(): 
            cluster[nb] = []
    for i in range(len(labels)):
        nb = labels[i]
        cluster[nb].append(entry_list[i])
    
    
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    plt.show()
    """
        
    







