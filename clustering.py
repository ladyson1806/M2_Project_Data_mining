#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import division
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import sys
import math
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import sklearn



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
    then format every GO for every protein as a list. It also makes sure that mass and length
    are in the correct format
    """
    onto_bp = "Gene ontology (biological process)"
    onto_mf = "Gene ontology (molecular function)"
    onto_cc = "Gene ontology (cellular component)"
    for protein in protein_list:
        if not protein.has_key(onto_bp):
            protein[onto_bp] = ""
        if not protein.has_key(onto_mf):
            protein[onto_mf] = ""
        if not protein.has_key(onto_cc):
            protein[onto_cc] = ""
        protein[onto_bp] = protein[onto_bp].split(";")   
        protein[onto_mf] = protein[onto_mf].split(";")
        protein[onto_cc] = protein[onto_cc].split(";")
        protein["Mass"] = float(protein["Mass"].replace(",", ""))
        protein["Length"] = float(protein["Length"].replace(",", ""))

def geneOntoAsList(protein_list, gene_ontology):
    """ 
    For every protein, get the list of a particular GO, and its entry ID.
    """
    go_list = []
    entry_list = []
    for protein in protein_list:
        go_list.append(protein[gene_ontology])
        entry_list.append(protein["Protein names"])
    return go_list, entry_list
        

def calcDistance(protein1, protein2):
    """
    A function to calculate the distance between two proteins considering a particular GO. 
    It calculates the number of shared GO (shared_nodes), and unique GO. 
    Operation:
    1 - (number of shared nodes) / (number of shared nodes + number of unique nodes)
    """
    shared_nodes = 0
    unique_nodes = 0
    if len(protein1) >= len(protein2):
        for go in protein1:
            if go in protein2:
                shared_nodes+=1
            else:
                unique_nodes+=1
    else:
        for go in protein2:
            if go in protein1:
                shared_nodes+=1
            else:
                unique_nodes+=1
    try:        
        return round(1 - (shared_nodes / (shared_nodes + unique_nodes)), 3)
    except (ZeroDivisionError):
        return 0
    
def generateDistanceMatrix(protein_list):
    """
    A function to calculate a distance matrix between proteins.
    """
    distance_matrix = []
    for protein in protein_list:
        distance_list = []
        for protein_checked in protein_list:
            distance_list.append(calcDistance(protein, protein_checked))
        distance_matrix.append(distance_list)
    return distance_matrix
            

protein_list = []
protein_list = loadfile(protein_list)
formatOntoData(protein_list)

go_list_mol, entry_list = geneOntoAsList(protein_list, "Gene ontology (molecular function)")
distance_matrix = generateDistanceMatrix(go_list_mol)


db = DBSCAN(eps=0.5, min_samples=10, metric="precomputed").fit(distance_matrix)
labels = db.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
print n_clusters_
cluster = {}
for nb in labels:
    if not nb in cluster.keys(): 
        cluster[nb] = []

for i in range(len(labels)):
    nb = labels[i]
    cluster[nb].append(entry_list[i])

f = open("results_mf.txt", 'w')
for i in cluster.keys():
    f.write("Cluster %s : " % i)
    for entry in cluster[i]:
       f.write(entry)
       f.write("\n")
    f.write("\n\n")
f.close()
     







