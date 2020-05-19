#!/usr/bin/env python3
import torch
import os
import sys


path = '/media/giacomo/Data/hierarchy_paper/projects/poincare-embeddings/output/'
if len(sys.argv)>1:
    path=sys.argv[1]

files = []
# r=root, d=directories, f = files
for r, d, f in os.walk(path):
   for file in f:
         dirs = torch.load(os.path.join(r, file))
         f = open(os.path.join(r, file)+".txt", "w")
         f.write(str(len(dirs["embeddings"][0]))+"\n")
         for x,y in zip(dirs["objects"],dirs["embeddings"]):
            f.write(str(x)+" "+(" ".join(map(lambda z : str(float(z)), y)))+"\n")
         f.close()
