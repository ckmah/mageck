#!/usr/bin/env python
""" MAGeCK count module
Copyright (c) 2014 Wei Li, Han Xu, Xiaole Liu lab 
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@version: $Revision$
@author:  Wei Li 
@contact: li.david.wei AT gmail.com
"""
from __future__ import print_function

import sys
import argparse
import math
import logging
import string
from mageck.testVisualCount import *
from mageck.mageckCountIO import *
from mageck.mageckCountNorm import *


def mageckcount_parseargs():
  """
  Parse arguments. Only used when mageckCount.py is executed directly.
  """
  parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
  
  parser.add_argument('-l','--list-seq',required=True,help='A file containing the list of sgRNA names, their sequences and associated genes. Support file format: csv and txt.')
  parser.add_argument('--sample-label',default='',help='Sample labels, separated by comma (,). Must be equal to the number of samples provided. Default "sample1,sample2,...".')
  parser.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  parser.add_argument('--trim-5',type=int,default=0,help='Length of trimming the 5\' of the reads. Default 0')
  parser.add_argument('--sgrna-len',type=int,default=20,help='Length of the sgRNA. Default 20')
  parser.add_argument('--count-n',action='store_true',help='Count sgRNAs with Ns. By default, sgRNAs containing N will be discarded.')
  parser.add_argument('--fastq',nargs='+',help='Sample fastq files, separated by space; use comma (,) to indicate technical replicates of the same sample. For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq" indicates two samples with 2 technical replicates for each sample.')

  
  args=parser.parse_args()
  
  
  return args

def mageckcount_checkargs(args):
  """
  Check args
  """
  if hasattr(args,'fastq') and args.fastq :
    if args.list_seq == None:
      logging.error('No library file specified.')
      sys.exit(-1)
  if args.sample_label!='':
    nlabel=args.sample_label.split(',')
    #nfq=args.fastq.split(',')
    nfq=(args.fastq)
    if len(nlabel)!=len(nfq):
      logging.error('The number of labels ('+str(nlabel)+') must be equal to the number of fastq files provided.')
      sys.exit(-1)
  return 0


def mageckcount_revcomp(x):
  '''
  Reverse complement
  '''
  return x.translate(string.maketrans("ACGT","TGCA"))[::-1]

def mageckcount_mergedict(dict0,dict1):
  '''
  Merge all items in dict1 to dict0.
  '''
  nsample=0
  if len(dict0)>0:
    nsample=len(dict0[dict0.keys()[0]])
  for (k,v) in dict0.items():
    if k in dict1:
      v+=[dict1[k]]
    else:
      v+=[0]
  for (k,v) in dict1.items():
    if k not in dict0:
      if nsample>0:
        dict0[k]=[0]*nsample
      else:
        dict0[k]=[]
      dict0[k]+=[v]
  # return dict0



def mageckcount_printdict(dict0,args,ofile,ounmappedfile,sgdict,datastat,sep='\t'):
  '''
  Write the table count to file
  '''
  if args.fastq is not None:
    allfastq=args.fastq
    nsample=len(allfastq)
    slabel=[datastat[f.split(',')[0]]['label'] for f in allfastq]
  elif args.count_table is not None:
    nsample=len(datastat)
    slabel=datastat.keys()
  # print header
  print('sgRNA'+sep+'Gene'+sep+sep.join(slabel),file=ofile)
  # print items
  if len(sgdict)==0:
    for (k,v) in dict0.items():
      print(k+sep+'None'+sep+sep.join([str(x) for x in v]),file=ofile)
  else:
    for (k,v) in dict0.items():
      if k not in sgdict: # only print those in the genedict
        if ounmappedfile != None:
          print(sep.join([k,k])+sep+sep.join([str(x) for x in v]),file=ounmappedfile)
        continue
      sx=sgdict[k]
      print(sep.join([sx[0],sx[1]])+sep+sep.join([str(x) for x in v]),file=ofile)
    # print the remaining counts, fill with 0
    for (k,v) in sgdict.items():
      if k not in dict0:
        print(sep.join([v[0],v[1]])+sep+sep.join(["0"]*nsample),file=ofile)

def mageck_printdict(dict0,args,sgdict,sampledict,sampleids):
  """Write the normalized read counts to file
  
  Parameters
  ----------
  dict0 : dict
    a {sgRNA: [read counts]} structure
  args : class
    a argparse class
  sgdict: dict
    a {sgrna:gene} dictionary
  sampledict: dict
    a {sample name: index} dict
  sampleids: list
    a list of sample index. Should include control+treatment
  
  """
  # print header
  # print items
  dfmt="{:.5g}"
  ofile=open(args.output_prefix+'.normalized.txt','w')
  # headers
  mapres_list=['']*len(sampledict)
  for (k,v) in sampledict.items():
    mapres_list[v]=k
  if len(sampledict)>0:
    cntheader=[mapres_list[x] for x in sampleids]
  else:
    cntheader=None
  logging.info('Writing normalized read counts to '+args.output_prefix+'.normalized.txt')
  if cntheader !=None:
    print('sgRNA\tGene\t'+'\t'.join(cntheader),file=ofile)
  if len(sgdict)==0:
    for (k,v) in dict0.items():
      print(k+'\t'+'None'+'\t'+'\t'.join([str(x) for x in v]),file=ofile)
  else:
    for (k,v) in dict0.items():
      if k not in sgdict: # only print those in the genedict
        logging.warning(k+' not in the sgRNA list')
        continue
      print('\t'.join([k,sgdict[k]])+'\t'+'\t'.join([str(x) for x in v]),file=ofile)
    # print the remaining counts, fill with 0
  ofile.close()




def mageckcount_checklists(args):
  """
  Read sgRNA library file
  Including sgRNAs and associated sequences and lists in csv or txt file
  format: sgRNAid  seq  geneid
  """
  genedict={}
  hascsv=False
  if args.list_seq.upper().endswith('CSV'):
    hascsv=True
  n=0
  seqdict={}
  ndup=0
  for line in open(args.list_seq):
    if hascsv:
      field=line.strip().split(',')
    else:
      field=line.strip().split()
    n+=1
    if field[0] in genedict:
      logging.warning('Duplicated sgRNA label '+field[0]+' in line '+str(n)+'. Skip this record.')
      continue
    if len(field)<3:
      logging.warning('Not enough field in line '+str(n)+'. Skip this record.')
      continue
    sgrnaseq=field[1].upper()
    if n==1:
      import re
      if re.search('[^ATCG]',sgrnaseq) is not None:
        logging.info('Header line of the library file detected; skip the first line ...')
        continue
    if hasattr(args,'reverse_complement') and args.reverse_complement:
      sgrnaseq=mageckcount_revcomp(sgrnaseq)
    if sgrnaseq in seqdict:
      # logging.warning('Duplicated sgRNA sequence '+field[1]+' in line '+str(n)+'. Skip this record.')
      ndup+=1
      continue
    genedict[field[0]]=(sgrnaseq,field[2])
    seqdict[sgrnaseq]=1
  logging.info('Loading '+str(len(genedict))+' predefined sgRNAs.')
  logging.warning('There are '+str(ndup)+' sgRNAs with duplicated sequences.')
  return genedict

def mageckcount_printstat(args,datastat):
  '''
  Write data statistics to PDF file 
  '''
  for (k,v) in datastat.items():
    logging.info('Summary of file '+k+':')
    for (v1,v2) in v.items():
      logging.info(str(v1)+'\t'+str(v2))
  # write to table
  crv=VisualRCount()
  crv.setPrefix(args.output_prefix)
  alllabels={}
  for (fq, fqstat) in datastat.items():
    crv.fastqfile+=[fq]
    if 'label' in fqstat:
      crv.fastqlabels+=[fqstat['label']]
      alllabels[fqstat['label']]=1
    else:
      crv.fastqlabels+=['NA']
    if 'reads' in fqstat:
      crv.reads+=[fqstat['reads']]
    else:
      crv.reads+=[0]
    if 'mappedreads' in fqstat:
      crv.mappedreads+=[fqstat['mappedreads']]
    else:
      crv.mappedreads+=[0]
    if 'totalsgrnas' in fqstat:
      crv.totalsgrnas+=[fqstat['totalsgrnas']]
    else:
      crv.totalsgrnas+=[0]
    if 'zerosgrnas' in fqstat:
      crv.zerocounts+=[fqstat['zerosgrnas']]
    else:
      crv.zerocounts+=[0]
    if 'giniindex' in fqstat:
      crv.gini+=[fqstat['giniindex']]
    else:
      crv.gini+=[0.0]
  #
  crv.startRTemplate()
  crv.writeCountSummary()
  # write to TXT file
  crv.writeCountSummaryToTxt(args.output_prefix+'.countsummary.txt')
  # write to PDF file
  outcsvfile=args.output_prefix+'.count_normalized.txt'
  crv.insertReadCountBoxPlot(os.path.basename(outcsvfile))
  if len(alllabels) > 1:
    crv.insertPCAPlot(os.path.basename(outcsvfile))
    crv.insertClusteringPlot(os.path.basename(outcsvfile))
  crv.closeRTemplate()
  if hasattr(args,"pdf_report") and args.pdf_report:
    if hasattr(args,"keep_tmp") :
      crv.generatePDF(keeptmp=args.keep_tmp)
    else:
      crv.generatePDF()

def mageckcount_processfastq(args,genedict,sgdict):
  """
  Main entry for fastq processing
  """
  # listfq=args.fastq.split(',')
  listfq=[[z for z in x.split(',')] for x in args.fastq]
  nsample=len(listfq)
  datastat={}
  # check labels
  alllabel=args.sample_label
  if alllabel=='':
    slabel=['sample'+str(x) for x in range(1,nsample+1)]
  else:
    slabel=alllabel.split(',')
  for i in range(nsample):
    for fi in listfq[i]:
      datastat[fi]={}
      datastat[fi]['label']=slabel[i]
  alldict={}
  # go through the fastq files
  for filenamelist in listfq:
    dict0={}
    for filename in filenamelist: # technical replicates; should be merged together
      dict00={}
      if filename.upper().endswith('BAM'):
        mageckcount_processonefile_bam(filename,args,dict00,sgdict,datastat[filename])
      elif filename.upper().endswith('SAM'):
        mageckcount_processonefile_sam(filename,args,dict00,sgdict,datastat[filename])
      else:
        mageckcount_processonefile(filename,args,dict00,sgdict,datastat[filename])
      for (k,v) in dict00.items():
        if k not in dict0:
          dict0[k]=0
        dict0[k]+=v
    mageckcount_mergedict(alldict,dict0)
  # write to file
  ofilel=open(args.output_prefix+'.count.txt','w')
  if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
    ounmappedfilel=open(args.output_prefix+'.unmapped.txt','w')
  else:
    ounmappedfilel=None
  mageckcount_printdict(alldict,args,ofilel,ounmappedfilel,sgdict,datastat)
  ofilel.close()
  if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
    ounmappedfilel.close()
  # write the median normalized read counts to csv file
  if len(sgdict)>0:
    allmappeddict={k:v for (k,v) in alldict.items() if k in sgdict} # only keep those with known sgRNAs
  else:
    allmappeddict=alldict
  return (allmappeddict,datastat)

def getcounttablefromfile(filename):
  """
  read count table from file
  Returns:
  ---------------
  x: dict
    {sgrna:[read counts]} 
  y: dict
    {sgrna:gene}
  z: dict
    z={sample_id:index}
  """
  gtab={}
  mapptab={}
  sampleids={}
  nline=0
  nfield=-1
  # if it is CSV file
  hascsv=False
  if filename.upper().endswith('.CSV'):
    hascsv=True
  logging.info('Loading count table from '+filename+' ')
  for line in open(filename):
    nline+=1
    if nline % 100000 == 1:
      logging.info('Processing '+str(nline)+' lines..')
    try:
      if hascsv==False:
        field=line.strip().split()
      else:
        field=line.strip().split(',')
      if len(field)<3:
        logging.warning('Line '+str(nline)+' of the read count table has fewer than 3 columns. Skip this line ...')
      sgid=field[0]
      geneid=field[1]
      # check if duplicate sgRNA IDs are detected
      if sgid in gtab:
        logging.warning('Duplicated sgRNA IDs: '+sgid+' in line '+str(nline)+'. Skip this record.')
        continue
      sgrecs=[float(x) for x in field[2:]]
      # check the number of fields
      if nfield!=-1 and len(sgrecs)!=nfield:
        logging.error('Error: incorrect number of dimensions in line '+str(nline)+'. Please double-check your read count table file.')
        sys.exit(-1)
      nfield=len(sgrecs)
      gtab[sgid]=sgrecs
      mapptab[sgid]=geneid
    except ValueError:
      if nline!=1:
        logging.warning('Parsing error in line '+str(nline)+'. Skip this line.')
      else:
        logging.debug('Parsing error in line '+str(nline)+' (usually the header line). Skip this line.')
        ids=field[2:]
        for i in range(len(ids)):
          sampleids[ids[i]]=i
      continue
  logging.info('Loaded '+str(len(gtab))+' records.')
  return (gtab,mapptab,sampleids)

def mageckcount_processcounttable(args,genedict,sgdict):
  """
  Main entry for count table processing
  Return value:
    gtab
        A {sgRNAID:[counts]} dict
    datastat
        Data statistics
    mapptab
        A {sgRNAID:geneID} dict
  """
  (gtab,mapptab,sampleids)=getcounttablefromfile(args.count_table)
  # check the consistency between provided library and count table
  nmiss=0
  nmissinlib=0
  if len(genedict)>0:
    for sgid in gtab.keys():
      if sgid not in genedict:
        nmiss+=1
    if nmiss>0:
      logging.warning(str(nmiss)+' sgRNA IDs in the count table are not in the library. Please double check.') 
    for sgid in genedict.keys():
      if sgid not in gtab:
        nmissinlib+=1
    if nmissinlib>0:
      logging.warning(str(nmissinlib)+' sgRNA IDs in the library are not in the count table. These sgRNAs will be counted as missing sgRNAs.') 
  # get the statistics of datasets
  datastat={}
  for (sampleid, sampleindex) in sampleids.items():
    datastat[sampleid]={}
    datastat[sampleid]['label']=sampleid
    nzerosg=0
    ntotalsgcount=0
    nrdcnt=[]
    for (sgid, cntv) in gtab.items():
      cval=cntv[sampleindex]
      if cval==0:
        nzerosg+=1
      ntotalsgcount+=cval
      nrdcnt+=[math.log(cval+1.0)]
    datastat[sampleid]['mappedreads']=ntotalsgcount
    datastat[sampleid]['reads']=ntotalsgcount
    datastat[sampleid]['zerosgrnas']=nzerosg+nmissinlib
    if len(genedict)>0:
      datastat[sampleid]['totalsgrnas']=len(genedict)
    else:
      datastat[sampleid]['totalsgrnas']=len(gtab)
    datastat[sampleid]['giniindex']=mageckcount_gini(nrdcnt)
  # end for loop for samples
  
  return (gtab,datastat,mapptab)


def mageckcount_main(args):
  """
  Main entry for mageck count module
  """
  # check arguments
  mageckcount_checkargs(args)
  # process gene dicts
  genedict={}
  if args.list_seq is not None:
    genedict=mageckcount_checklists(args) # sgid:(seq,gene)
  # save sgRNA ID and gene name
  sgdict={} #
  for (k,v) in genedict.items():
    sgdict[v[0]]=(k,v[1]) # {seq:(sgid,gene)
  if hasattr(args,'count_table') and args.count_table != None:
    # treat it as a count table
    (allmappeddict,datastat,mapptab)=mageckcount_processcounttable(args,genedict,sgdict)
    # note that the key of allmappeddict is sgRNA ID
    # if library file is provided, we need to change sgdict to make it consistent with other situations (like fastq file)
    sgdict={k:(k,v) for (k,v) in mapptab.items()}
    # print('sg size:'+str(len(sgdict)))
  else:
    # check the listed files: fastq/sam/bam files provided
    (allmappeddict,datastat)=mageckcount_processfastq(args,genedict,sgdict)
    # note that the key of allmappeddict is sgRNA sequence
  
  # normalize read counts
  if hasattr(args,"norm_method"):
    normmethod=args.norm_method
  else:
    normmethod="median"
  if hasattr(args,"control_sgrna"):
    ctrlsg=args.control_sgrna
  else:
    ctrlsg=None
  medalldict=normalizeCounts(allmappeddict,method=normmethod,controlsgfile=ctrlsg)
  ofilel=open(args.output_prefix+'.count_normalized.txt','w')
  mageckcount_printdict(medalldict,args,ofilel,None,sgdict,datastat,sep='\t')
  ofilel.close()
  # print statistics
  mageckcount_printstat(args,datastat)
  return 0



if __name__ == '__main__':
  try:
    args=mageckcount_parseargs()
    mageckcount_main(args)
  except KeyboardInterrupt:
    sys.stderr.write("Interrupted.\n")
    sys.exit(0)



