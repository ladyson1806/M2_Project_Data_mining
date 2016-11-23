#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import division
import sys
import math


def loadfile(protein_list):
    """
    Load a tab separated file in a list
    """
    try:
        f = open(sys.argv[1], 'r')
    except:
        print "Le fichier ", sys.argv[1], " est introuvable"
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


def calcDistance(protein1, protein2, gene_ontology):
    """
    A function to calculate the distance between two proteins considering a particular GO. 
    It calculates the number of shared GO (shared_nodes), and unique GO. 
    Operation:
    1 - (number of shared nodes) / (number of shared nodes + number of unique nodes)
    """
    shared_nodes = []
    unique_nodes = []
    for GO in protein1[gene_ontology]:
        if GO in protein2[gene_ontology]:
            shared_nodes.append(GO)
        else:
            unique_nodes.append(GO)
    try:        
        return round(1 - (len(shared_nodes) / (len(shared_nodes) + len(unique_nodes))), 3)
    except (ZeroDivisionError):
        return 1
    
def generateDistanceMatrix(protein_list, gene_ontology):
    """
    A function to calculate a distance matrix between proteins considering a particular GO.
    """
    distance_matrix = []
    for protein in protein_list:
        distance_list = []
        for protein_checked in protein_list:
            distance_list.append(calcDistance(protein, protein_checked, gene_ontology))
        distance_matrix.append(distance_list)
    return distance_matrix
            
        
                
protein_list = []
protein_list = loadfile(protein_list)
formatOntoData(protein_list)
print generateDistanceMatrix(protein_list[0:5000], "Gene ontology (cellular component)")









