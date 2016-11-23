#!/usr/bin/env python
import sys

try:
    f = open(sys.argv[1], 'r')
except:
    print "Le fichier ", sys.argv[1], " est introuvable"
    sys.exit()
    
line = f.readline()
line = line.rstrip("\n")
list_key = line.split("\t")

nbArg = []
for arg in line:
    nbArg.append(0.0)
total = 0
for line in f:
    total+=1
    line = line.split("\t")
    i = 0
    for arg in line:
        if (arg!="" and arg!="\n"):
            nbArg[i]+=1
        i+=1
f.close()
for i in range(0, len(list_key)):
    print list_key[i], " : ", (nbArg[i]/total)*100, "%"
print total
