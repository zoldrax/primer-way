#!/usr/bin/python
import getopt
import gzip
import math
import os
import re
import subprocess
import sys

import primer3
import pysam
import vcf

try:
              opts,args = getopt.getopt(sys.argv[1:],"hvG:R:V:i:o:p:r:a:n:s:d:")
except:
              pass

textwidth = 60
maxtargetlength = 50
#maxtargetlength = 55000
flanking = 300
#flanking = 500 
overlap = 20
penalty_PCR_product = 100
penalty_base_level = 3
pertarget = 300
pertntblast = 600

p3seq={
              'SEQUENCE_ID': 'example',
              'SEQUENCE_TEMPLATE': ''
}
p3prim={
              'PRIMER_QUALITY_RANGE_MAX': 40,
              'PRIMER_WT_SEQ_QUAL': 0.03,
              'PRIMER_WT_END_QUAL': 0.1,
              'PRIMER_PICK_LEFT_PRIMER': 1,
              'PRIMER_PICK_RIGHT_PRIMER': 1,
              'PRIMER_MAX_END_GC': 3,
              'PRIMER_MAX_POLY_X': 4,
              'PRIMER_PRODUCT_SIZE_RANGE': [[50,300]],
              'PRIMER_MAX_TM': 65,
              'PRIMER_WT_SIZE_LT': 0.5,
              'PRIMER_WT_SIZE_GT': 0.5,
              'PRIMER_WT_TM_LT': 0.1,
              'PRIMER_WT_TM_HT': 0.05,
              'PRIMER_PAIR_WT_DIFF_TM':0.3,
              'PRIMER_WT_HAIRPIN_TH': 2,
              'PRIMER_WT_SELF_ANY_TH': 1,
              'PRIMER_WT_SELF_END_TH': 2,
              'PRIMER_PAIR_WT_COMPL_ANY_TH': 1,
              'RIMER_PAIR_WT_COMPL_END_TH': 2,
              'PRIMER_EXPLAIN_FLAG': 1,
              'PRIMER_NUM_RETURN': pertarget
}

verbose = False
jobname="Exon"
prot_id = ""
region = ""
deletion = ""
ascoord = ""
output_d = ""
input_d = ""
reference_file = ""
reference_gff_file = ""
variant_file = ""
searches=[]
searches_coord=[]
startexon=1

for o, a in opts:
              if o == "-v":
                            verbose = True
              elif o in ("-h", "--help"):
                            print "-R : reference unzipped fasta file"
                            print "-G : reference gzipped GFF description file"
                            print "-V : variants gzipped VCF dbSNP file"
                            print "-v : verbose mode for debug"
                            print "-i : input directory"
                            print "-o : output directory"
                            print "-p : protein ID for automatic exon CDS extraction (example NP_000442.1)"
                            print "-d : deletion region coordinates (example: 11:121212-454545)"                          
                            print "-r : target region coordinates (example: 11:121212-454545)"
                            print "-a : target nucleotide coordinate for allele specific (example: 11:121212)"
                            print "-n : prefix for primer names"
                            print "-s : start from exon number (default 1)"                           
                            sys.exit()
                            
              elif o in ("-s"):
                            try:
                                 startexon = int(a)
                            except:
                                 pass
              elif o in ("-p"):
                            prot_id = a
              elif o in ("-r"):
                            region = a
              elif o in ("-d"):
                            deletion= a
              elif o in ("-a"):
                            ascoord = a                            
              elif o in ("-o"):
                            output_d = a
              elif o in ("-i"):
                            input_d = a
              elif o in ("-R"):
                            reference_file = a
              elif o in ("-G"):
                            reference_gff_file = a
              elif o in ("-V"):
                            variant_file = a
              elif o in ("-n"):
                            jobname = a
                            
if (output_d==""):
              print "Error: Output directory is needed"
              sys.exit()
try:
                  os.stat(output_d)
except:
                  try:
                                os.mkdir(output_d)
                  except:
                                print "Error: Can not create output directory"

                                sys.exit()

if (ascoord!="")&(reference_file!=""):
              exons=[]
              regreg = re.findall("(.+):(\d+)",ascoord)
              with pysam.FastaFile(reference_file) as fh:
                            try:
                                          start=int(regreg[0][1])
                                          end=int(regreg[0][1])
                            except:
                                          start=0
                                          end=0
                            searches.append([fh.fetch(region=regreg[0][0]+":"+str(start-50)+"-"+str(end+500)),regreg[0][0],range(start-50,end+501)])
                            searches.append([fh.fetch(region=regreg[0][0]+":"+str(start-500)+"-"+str(end+50)),regreg[0][0],range(start-500,end+51)])
              prot_id = ""
              region = ""
                                
if (deletion!="")&(reference_file!=""):
              exons=[]
              regreg = re.findall("(.+):(\d+)-(\d+)",deletion)
              with pysam.FastaFile(reference_file) as fh:
                            try:
                                          start=int(regreg[0][1])
                                          end=int(regreg[0][2])
                            except:
                                          start=0
                                          end=0
                            searches.append([fh.fetch(region=regreg[0][0]+":"+str(start-flanking-21)+"-"+str(start-1))+fh.fetch(region=regreg[0][0]+":"+str(end+1)+"-"+str(end+flanking+21)),regreg[0][0],range(start-flanking-21,start)+range(end+1,end+flanking+22)])
              prot_id = ""
              region = ""

if (region!="")&(reference_file!=""):
              exons=[]
              regreg = re.findall("(.+):(\d+)-(\d+)",region)
              with pysam.FastaFile(reference_file) as fh:
                            try:
                                          start=int(regreg[0][1])-flanking-20
                                          end=int(regreg[0][2])+flanking+20
                            except:
                                          start=0
                                          end=0
                            searches.append([fh.fetch(region=regreg[0][0]+":"+str(start)+"-"+str(end)),regreg[0][0],range(start,end+1)])
              prot_id = ""
              
if (prot_id!="")&(reference_gff_file!=""):
              with gzip.open(reference_gff_file,"r") as csvfile:
                            print ("Searching exons for transcript ID = "+prot_id+" ...")
                            exons = re.findall("(chr[\dXY]+)\t.+\tCDS\t(\d+)\t(\d+)\t.+protein_id="+prot_id, csvfile.read())
              print ("Found "+ str(len(exons))+ " exons.")
              fex = open(output_d+"/Exon_list.txt", 'w')
              fex.write("Coord without flanking "+str(flanking+20)+"\n")
              for i in exons:
                            fex.write("\t".join(i)+"\n")
                            with pysam.FastaFile(reference_file) as fh:
                                          try:
                                                        start=int(i[1])-flanking-20
                                                        end=int(i[2])+flanking+20
                                          except:
                                                        start=0
                                                        end=0
                                          searches.append([fh.fetch(region=i[0]+":"+str(start)+"-"+str(end)),i[0],range(start,end+1)])
              fex.flush()
              fex.write("\n")
              
fsint = open(output_d+"/resultforsintes.txt", "w")
exon=0              
for search in searches:
              exon+=1
              if ((region=="")&(ascoord=="")):
                            print "Exon "+str(exon)+" of "+str(len(exons))+" is assaying..."
              if (exon < startexon):
                            print "Skipping..."
                            continue
              pairs={}
              sequence = re.sub('[^ACGTacgt]+', '',search[0])
              lqual=[]
              if (variant_file!=""):
                            vcf_reader = vcf.Reader(filename = variant_file)
                            for i in search[2]:
                                          try:
                                                        variants = vcf_reader.fetch(search[1],i-1,i)
                                          except:
                                                        variants = vcf_reader.fetch(search[1][3:],i-1,i)
                                          qual=0.0001
                                          for variant in variants:
                                                        try:
                                                                      q = 1.0-float(variant.INFO["CAF"][0])
                                                                      if (q>qual):qual=q
                                                        except:
                                                                      pass
                                          lqual.append(int(-10.0*math.log10(qual)))

              rightedge = len(sequence)-flanking
              leftedge = flanking
              if (ascoord!=""): rightedge=leftedge+1
              print len(sequence)
              print len(lqual)
              
              for i in range (leftedge, rightedge, 50):
                            targetlength=rightedge-i
                            if targetlength > maxtargetlength:
                                          targetlength=maxtargetlength

                            p3seq['SEQUENCE_TEMPLATE']=sequence
                            
                            if (variant_file!=""):
                                          p3seq['SEQUENCE_QUALITY']=lqual
                            else:
                                          p3prim['PRIMER_WT_SEQ_QUAL']=0
                                          p3prim['PRIMER_WT_END_QUAL']=0
                                          
                            if (ascoord==""):
                                          if (deletion==""):
                                                        p3seq['SEQUENCE_TARGET']=[i,targetlength]
                                          else:
                                                        p3seq['SEQUENCE_TARGET']=[(flanking+20),2]
                            else:
                                          if (exon<2):
                                                        p3seq['SEQUENCE_FORCE_LEFT_END']=50
                                          else:
                                                        p3seq['SEQUENCE_FORCE_RIGHT_END']=500
                                          


                            result = primer3.bindings.designPrimers(p3seq, p3prim)

                            for i in xrange(pertarget):
                                          si=str(i)
                                          spair=(result['PRIMER_LEFT_'+si+'_SEQUENCE'],result['PRIMER_RIGHT_'+si+'_SEQUENCE'])
                                          pairs[spair]=[
                                                        result['PRIMER_LEFT_'+si][0]+result['PRIMER_LEFT_'+si][1]+1,
                                                        result['PRIMER_RIGHT_'+si][0]-result['PRIMER_RIGHT_'+si][1],
                                                        result['PRIMER_PAIR_'+si+'_PENALTY'],
                                                        0
                                                        ]
 
                            
                            for pair in pairs:
                                          if ((pair[0][-1] in set("ATat")) | (pair[1][-1] in set("ATat"))):
                                                        pairs[pair][2]+=1
              print ("primers: " + str(len(pairs)))
              s=""
              i=0
              complete=0
              for pair in pairs:
                            i+=1
                            complete+=1
                            s+="_".join(pair)+"\t"+pair[0]+"\t"+pair[1]+"\n"
                            if ((i>=pertntblast)|(pair==pairs.keys()[-1])):
                                          if verbose : print s
                                          print ("Tntblast: "+str(complete)+" of "+ str(len(pairs)))
                                          proc = subprocess.Popen(("nice -1 tntblast -i /dev/stdin -d "+ reference_file +" -g -5 -l 1000 --temperature 330 --primer-clamp 1").split(),stdout=subprocess.PIPE, stdin=subprocess.PIPE)
                                          proc.stdin.write(s)
                                          s=""
                                          i=0
                                          proc.stdin.close()
                                          result = proc.stdout.read()
                                          proc.wait()
                                          pcrs = re.findall("name = (\w+)",result)
                                          for pcr in pcrs:
                                                        pairs[tuple(pcr.split('_'))][3]+=1
                                                        pairs[tuple(pcr.split('_'))][2]+=penalty_PCR_product
              print ("Computing best way...")
              ways=[]
              f = open(output_d+"/Exon_"+str(exon)+'_specprimers.txt', 'w')
              f.write("left edge "+str(leftedge)+", right edge "+str(rightedge)+"\n")
              for pair in pairs:
                            if (pairs[pair][3]==0)&(deletion==""):
                                          pairs[pair][2]+=2*penalty_PCR_product
                                          
              pairssort=sorted(pairs, key=lambda d: pairs[d][2])
              for pair in pairssort:
                            f.write(pair[0]+"\t"+pair[1]+"\t"+"\t"+str(pairs[pair][0])+"\t"+str(pairs[pair][1])+"\t"+str(pairs[pair][2])+"\n")
              f.close()

              graph ={}
              graph["START"]=[]
              for pair1 in pairs:
                            if (pairs[pair1][0] < leftedge): graph["START"].append(pair1)
                            if (pairs[pair1][1] > rightedge):
                                          graph[pair1]=graph.get(pair1,[])
                                          graph[pair1].append("END")
                            for pair2 in pairs:
                                          if (pairs[pair2][0] < (pairs[pair1][1] - overlap)):
                                              graph[pair1]=graph.get(pair1,[])
                                              graph[pair1].append(pair2)

              pairs["END"]=[0,0,0,0];
              Visited={}
              ToVisit={"START":0}
              Paths = {"START":["START"]}
              while ToVisit:
                            v=min(ToVisit, key=lambda x: ToVisit[x])
                            Visited[v]=ToVisit[v]; del ToVisit[v]
                            for w in graph.get(v,[]):
                                          if (w not in Visited):
                                                        vwLength = Visited[v] + ((pairs[w][2]-penalty_PCR_product)**2) + penalty_base_level
                                                        if (w not in ToVisit) or (vwLength < ToVisit[w]):
                                                                      ToVisit[w] = vwLength
                                                                      Paths[w] = Paths[v]+[w]
              try:
                            bestway=Paths["END"][1:-1]
              except:
                            bestway=[]

              f = open(output_d+"/Exon_"+str(exon)+'_resultprimers.txt', 'w')
              f.write("left edge "+str(leftedge)+", right edge "+str(rightedge)+"\n")
              f.write("The best way consists "+ str(len(bestway))+" pairs.\n\n")
              print ("The best way consists "+ str(len(bestway))+" pairs.")
              penalty=0
              prpair=0
              specchar=dict()
              specchar[leftedge+20]="{"
              specchar[rightedge-20]="}"
              for pair in bestway:
                            prpair+=1
                            if (deletion!=""):
                                          penalty+=(pairs[pair][2]-0)**2
                            else:
                                          penalty+=(pairs[pair][2]-100)**2
                            f.write(pair[0]+"\t"+pair[1]+"\t"+"\t"+str(pairs[pair][0])+"\t"+str(pairs[pair][1])+"\t"+str(pairs[pair][2])+"\n")
                            f.flush()
                            prname=jobname
                            if (region=="") : prname+="_"+str(exon)
                            if len(bestway)>1:
                                          prname+="_"+str(prpair)
                            fsint.write((prname+"f,"+pair[0]+"\n").lower())  
                            fsint.write((prname+"r,"+pair[1]+"\n").lower())
                            fsint.flush()
                            if (prot_id!="")&(reference_gff_file!=""):
                                          fex.write(prname+"\t"+str(pairs[pair][2])+"\n")
                                          fex.flush()
                            specchar[pairs[pair][0]-len(pair[0])-1]="["
                            specchar[pairs[pair][0]-1]=">"
                            specchar[pairs[pair][1]+1]="<"
                            specchar[pairs[pair][1]+len(pair[1])+1]="]"
              f.write("\n")
              print ("Penalty summ sqr = "+ str(penalty))
              f.write ("Penalty saumm sqr = "+ str(penalty)+"\n\n")    
              n=0
              for i in xrange(0, len(sequence)):
                            f.write(specchar.get(i,"")+sequence[i])
                            if n>textwidth:
                                          f.write("\n")
                                          n=0
                            
                            n+=1
              f.write("\n\n{...} - coding DNA sequence (CDS)\n[...> - forward primer\n<...] - reverse primer")
              f.close()
fsint.close()
if (prot_id!="")&(reference_gff_file!=""): fex.close()
