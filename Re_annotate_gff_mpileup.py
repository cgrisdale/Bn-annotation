#!/usr/bin/python
#June 3rd, 2015
#Cameron Grisdale and Shannon Sibbald
#Input gff annotation file and alignment pileup file and re-annotate gene models based on alignment data

#Imports
import re
import sys
import itertools
import time
from collections import Counter
from operator import itemgetter

### Support Functions ###
def Packd(astring,adict,alist):
  newlist=[]
  if astring in adict: #if gid already in dict
    for i in adict[astring]:
      newlist.append(i)
    newlist.append(alist)
    adict[astring]=newlist
  else:
    adict[astring]=[alist]
  return adict

### Functions ###

#Go through pileup file and create transcripts from contiguous sequence, then check against gene models
#Go through gene models and check against alignments


def covgetter(covfile):
  '''Suvey pileup file for contiguous transcripts and store in dict'''
  pileup,key,currentpos = {},0,0
  with open(covfile) as f:
    for line in f:
      scaf,pos,con,cov,nts,quals = [n for n in line.strip().split('\t')]
      if int(pos) > (int(currentpos)+1): #there is a gap
        key=key+1
        tid=scaf+"_"+str(key) #transript name is scaffold name plus a number
        pileup[tid]=[(pos,con,cov,nts,quals)] #start new transcript
      else:
        pileup[tid]=pileup[tid]+[(pos,con,cov,nts,quals)] #add to existing transcript
      currentpos=pos
      #Packd(pos,pileup,mylist)
  print 'Windows Indexed:', len(pileup)
  return pileup

#Bigelowiella natans
def generate_maps_natans(gff_file): #Chlorarachniophyte.exons.clean.gff3 or Bnatans.exons.clean.gff3
  '''Returns the Gene Map, and the feature map, and a pre-initialized container'''
  print '\033[96mBnatans annotation file\033[0m'
  gene_map,feature_map,collector,stranded = {},{},{},{}
  with open(gff_file,'r') as f:
    for line in f:
      chrm,method,feat_type,start,end,score,strand,frame,attribute = line.strip().split('\t')
      if feat_type != 'exon':
	continue
      gene = attribute.split(';')[1].split('=')[1]
      stranded[gene] = strand
      gene_map.setdefault(chrm, {})[gene] = ''
      if gene not in feature_map:
	feature_map[gene] = [(int(start),int(end),'e')]
      else:
	feature_map[gene] = sorted(feature_map[gene] + [(int(start),int(end),'e')])
  #Resort Exons due to loss no numberical information.
  for gene, fmap in feature_map.items():
    if stranded[gene] == '+':
      newmap = [(feature[0],feature[1],feature[2],fmap.index(feature)+1) for feature in fmap]
      feature_map[gene] = newmap
    elif stranded[gene] =='-':
      newmap = [(feature[0],feature[1],feature[2],len(fmap)-fmap.index(feature)) for feature in fmap]
      feature_map[gene] = newmap
    #elif stranded[gene] == '.':
    #  newmap = [(feature[0],feature[1],feature[2],fmap.index(feature)+1) for feature in fmap]
    #  feature_map[gene] = newmap
  #Add Introns
  for gname,feature_set in feature_map.items():
    if len(feature_set) >=2:
      if feature_set[0][3] ==1:
	introns = [(feature_set[i-1][1]+1,feature_set[i][0]-1,'i',i) for i in range(1,len(feature_set))]
	feature_map[gname] = sorted(feature_map[gname]+introns)
      else:
	introns = [(feature_set[i-1][1]+1,feature_set[i][0]-1,'i',len(feature_set)-i) for i in range(1,len(feature_set))]
	feature_map[gname] = sorted(feature_map[gname]+introns)
  #Create Gene Map, initialize all values in the collector
  for xsomemap in gene_map.values():
    for k, v in xsomemap.items():
      xsomemap[k] = (feature_map[k][0][0],feature_map[k][-1][1])
      collector[k] = Counter()
  ##Fancy Output is a must
  all_loci,all_somes = sum([len(x) for x in gene_map.values()]),len(gene_map)
  print '\033[92mChromosomes Added\033[0m...', str(all_somes).rjust(7)
  print '\033[92mGene Loci Added\033[0m.....\033[0m',str(all_loci).rjust(7)
  #print '\033[96mDistribution of loci\033[0m'
  cols =0
  #for k, v in sorted(gene_map.items(),key=itemgetter(1), reverse=True):
  #  print '\033[92m'+k+'\033[0m'+'.'*(10-len(k)), str(len(v)).rjust(4), str((float(len(v))/all_loci)*100)[0:5]+'%',
  #  if cols == 3:
  #    print ''
  #    cols = 0
  #  else:
  #    cols += 1
  #print '\033[96mFinished with Distribution of loci\033[0m'
  ##fill in phantom contigs
  #for k in ['scaffold_'+ str(x) for x in range(0,900)]:
  # if k in gene_map:
  #   continue
  # else:
  #   gene_map[k] = {}
  #Prime dict for Bisecting
  for k, v in gene_map.items():
    gene_map[k] = sorted((value, key) for key, value in v.iteritems())
  return [gene_map,feature_map,collector,stranded]



if __name__ == '__main__':
  #Workflow
  #Intialized variables...
  mpileup,annotation = sys.argv[1],sys.argv[2]
  txn=covgetter(mpileup)
  gmap,fmap,col,strd=generate_maps_natans(annotation)
  
#  x=0
#  for k,v in txn.items():
#    x+=1
#    if x > 10:
#      break
#    else:
#      print k,v
