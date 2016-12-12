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

def formatData(protein_list):
    """
    A function to add Gene ontology if it's missing from data (maybe it's a tab format problem)
    then format every GO for every protein as a list (Gene ontology is None if it's empty).
    It also makes sure that mass and length are in the correct format.
    """
    entry_id = "Entry"
    onto = "Gene ontology (GO)"
    mass = "Mass"
    length = "Length"
    new_protein_list = []
    for protein in protein_list:
        new_formatted_protein = []
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
        new_formatted_protein.append(protein[entry_id])
        new_formatted_protein.append(protein[onto])
        new_formatted_protein.append(protein[mass])
        new_formatted_protein.append(protein[length])
        new_protein_list.append(new_formatted_protein)
    return new_protein_list

def getMinMaxMass(protein_list):
    """
    Return the maximum and minimum mass of proteins in a protein list
    """
    mass = 2 #index of mass in a protein list
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
    length = 3 #index of length in a protein list
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
    mass = 2 #index of mass in a protein list
    delta = 0
    dist = 0
    if not protein1[mass] is None and not protein2[mass] is None:
        delta = 1
        dist = abs(protein1[mass] - protein2[mass]) / (massMax - massMin)
    return delta, dist

def calcLengthDistance(protein1, protein2, lengthMin, lengthMax):
    length = 3 #index of length in a protein list
    delta = 0
    dist = 0
    if not protein1[length] is None and not protein2[length] is None:
        delta = 1
        dist = abs(protein1[length] - protein2[length]) / (lengthMax - lengthMin)
    return delta, dist
        

def calcGODistance(protein1, protein2):
    """
    A function to calculate the distance between two proteins considering a particular GO. 
    It calculates the number of shared GO (shared_nodes), and unique GO. 
    Operation:
    1 - (number of shared nodes) / (number of shared nodes + number of unique nodes)
    """
    GO = 1 #index of GO in a protein list
    delta = 0
    shared_nodes = 0
    unique_nodes = 0
    if not protein1[GO] is None and not protein2[GO] is None:
        delta = 1
        if len(protein1[GO]) >= len(protein2[GO]):
            for go in protein1[GO]:
                if go in protein2[GO]:
                    shared_nodes+=1
                else:
                    unique_nodes+=1
        else:
            for go in protein2[GO]:
                if go in protein1[GO]:
                    shared_nodes+=1
                else:
                    unique_nodes+=1
        return delta, round(1 - (shared_nodes / (shared_nodes + unique_nodes)), 3)
    else:
        return delta, 0

def calcDistance(protein1, protein2, massMin, massMax, lengthMin, lengthMax):
    w = 2
    deltaMass, massDistance = calcMassDistance(protein1, protein2, massMin, massMax)
    deltaLength, lengthDistance = calcLengthDistance(protein1, protein2, lengthMin, lengthMax)
    deltaGO, GODistance = calcGODistance(protein1, protein2)
    try:
        return round(((deltaMass*massDistance + deltaLength*lengthDistance + deltaGO*GODistance*w) /
                      (deltaMass + deltaLength + deltaGO*w)), 3)
    except ZeroDivisionError:
        return 1
    
def generateDistanceMatrix(protein_list):
    """
    A function to calculate a distance matrix between proteins.
    """
    massMin, massMax = getMinMaxMass(protein_list)
    lengthMin, lengthMax = getMinMaxLength(protein_list)
    distance_matrix = []
    for protein in protein_list:
        distance_list = []
        for protein_checked in protein_list:
            if protein is protein_checked:
                distance_list.append(0.0)
            else:
                distance_list.append(calcDistance(protein, protein_checked,
                                                  massMin, massMax, lengthMin, lengthMax))
        distance_matrix.append(distance_list)
    return distance_matrix


if __name__ == '__main__':
    protein_list = []
    protein_list = loadfile(protein_list)
    protein_list = formatData(protein_list)
    distance_matrix = generateDistanceMatrix(protein_list)
    
    """   
    db = DBSCAN(eps=0.5, min_samples=5, metric="precomputed").fit(distance_matrix)
    
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)


    cluster = {}
    for nb in labels:
        if not nb in cluster.keys(): 
            cluster[nb] = []
    for i in range(len(labels)):
        nb = labels[i]
        cluster[nb].append(protein_list[i][0])
    print distance_matrix[0]
    """
    
        
    







