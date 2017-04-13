""" Processing various file types for MAGeCK count
"""

from __future__ import print_function

import sys
import math
import logging
import string


def mageckcount_gini(x):
  '''
  Return the Gini index of an array
  Calculation is based on http://en.wikipedia.org/wiki/Gini_coefficient
  '''
  xs=sorted(x)
  n=len(xs)
  gssum=sum([ (i+1.0)*xs[i] for i in range(n)])
  ysum=sum(xs)
  if ysum==0.0:
    ysum=1.0
  gs=1.0-2.0*(n-gssum/ysum)/(n-1)
  return gs



def mageckcount_processonefile(filename,args,ctab,genedict,datastat):
  '''
  Go through one fastq file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A dictionary of sgRNA sequence and count
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})

  Return value
  -----------
  datastat
    a dictionary structure of statistics
  '''
  # ctab={}
  nline=0
  logging.info('Parsing FASTQ file '+filename+'...')
  nreadcount=0
  # checking possible sgRNA length
  lengthpool={}
  for k in genedict.keys():
    if len(k) not in lengthpool:
      lengthpool[len(k)]=0
    lengthpool[len(k)]+=1
  lengthpoolkeys=sorted(lengthpool.keys(),reverse=True)
  logging.info('Possible gRNA lengths:'+','.join([str(t) for t in lengthpoolkeys]))
  if filename.upper().endswith('.GZ'):
    import gzip
    openobj=gzip.open(filename,'rt')
  else:
    openobj=open(filename)
  for line in openobj:
    # line=line.encode('utf-8')
    nline=nline+1
    if nline%1000000==1:
      logging.info('Processing '+str(round(nline/1000000))+ 'M lines..')
    if nline%4 == 2:
      nreadcount+=1
      fseq=line.strip()
      if args.trim_5 >0:
        fseq=fseq[args.trim_5:]
      # check length
      # for l in lengthpool.keys():
      if len(genedict)==0:
        if len(fseq)<args.sgrna_len:
          continue
        fseq=fseq[:args.sgrna_len]
        if fseq.count('N')>0 and args.count_n==False:
          continue
        if fseq not in ctab:
          ctab[fseq]=0
        ctab[fseq]=ctab[fseq]+1
      else:
        findrecord=False
        for l in lengthpoolkeys: # iterate all possible lengths
          testl=l
          if len(fseq)<testl:
            continue
          fseqc=fseq[:testl]
          if fseqc.count('N')>0 and args.count_n==False:
            continue
          if fseqc not in genedict:
            continue
          else:
            if fseqc not in ctab:
              ctab[fseqc]=0
            ctab[fseqc]=ctab[fseqc]+1
            findrecord=True
            break
        # save unmapped file
        if args.unmapped_to_file and findrecord==False: 
          if len(fseq)<args.sgrna_len:
            continue
          fseqc=fseq[:args.sgrna_len]
          if fseqc.count('N')>0 and args.count_n==False:
            continue
          if fseqc not in ctab:
            ctab[fseqc]=0
          ctab[fseqc]=ctab[fseqc]+1
  # 
  openobj.close()
  # calculate statistics
  datastat['reads']=nreadcount
  # check if a library is provided
  if len(genedict)==0:
    datastat['mappedreads']=0
    datastat['totalsgrnas']=0
    datastat['zerosgrnas']=0
    datastat['giniindex']=1
  else:
    nmapped=0
    nrdcnt=[]
    for (k,v) in ctab.items():
      if k in genedict:
        nmapped+=v
        nrdcnt+=[math.log(v+1.0)]
    nzerosg=0
    for (k,v) in genedict.items():
      if k not in ctab:
        nzerosg+=1
        nrdcnt+=[math.log(0.0+1.0)]
    logging.info('mapped:'+str(nmapped))
    datastat['mappedreads']=nmapped
    datastat['totalsgrnas']=len(genedict);
    datastat['zerosgrnas']=nzerosg
    datastat['giniindex']=mageckcount_gini(nrdcnt)
  #return ctab
  return 0


def mageckcount_processonefile_bam(filename,args,ctab,genedict,datastat):
  '''
  Go through bam file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A dictionary of sgRNA sequence and count
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})

  Return value
  -----------
  datastat
    a dictionary structure of statistics
  
  # important note for the alignment
  1. Make sure 5' and 3' adapters are properly removed before mapping. Either use cutadapt or --trim5/--trim3 option in bowtie2.
  2. Make sure no reverse-complement mapping is allowed; otherwise, there will be multiple mappings for sgRNAs whose sequnces are reverse complemented. In bowtie2, --no-rc parameter should be specified.
  3. Carefully check the alignment strategy used in the aligner. Some sequences may map to multiple sgRNAs and are counted multiple times.
  4. When building index, some aligners (like bowtie2) will remove sgRNAs with identical sequence. This will create some warning messages "sgRNA in the BAM file does not match th provided library file".

  Reference: 
  BAM specification 
  https://samtools.github.io/hts-specs/SAMv1.pdf 
  Reference: Tao Liu's MACS2 
  https://github.com/taoliu/MACS/blob/master/MACS2/IO/Parser.pyx
  '''
  import sys 
  import gzip 
  import io 
  import math 
  import struct 
  from struct import unpack 
  """ 
  Encode table  
  ACMGRSVTWYHKDBNN -> [0,15] 
  Code int bit bit(reverse) int(reverse) 
    0  0000    
  A 1  0001    
  C 2  0010    
  M 3  0011    
  G 4  0100    
  R 5  0101    
  S 6  0110    
  V 7  0111    
  T 8  1000    
  W 9  1001    
  Y 10 1010    
  H 11 1011    
  K 12 1100    
  D 13 1101    
  B 14 1110    
  N 15 1111    
  """ 
  encodetable='NACMGRSVTWYHKDBN' 
  nline=0
  logging.info('Parsing BAM file '+filename+'...')
  nreadcount=0
  
  # start processing the bam file
  # open bam file 
  fhd=io.BufferedReader( gzip.open( filename, mode='rb' ) ) 
  # check the first 3 character must be BAM 
  fhd.seek(0) 
  magic_header = fhd.read( 3 ) 
  if magic_header!= "BAM": 
    logging.error('Error: not recognized BAM file: '+filename) 
    sys.exit(-1) 
  # check header 
  fhd.seek( 4 ) 
  header_len =  unpack( '<i', fhd.read( 4 ) )[ 0 ] 
  fhd.seek( header_len + fhd.tell() ) 
  # next, get chromosome, and check whether it matches the given genedict
  genedict_sgid={v[0]:k for (k,v) in genedict.items()}
  nc = unpack( '<i', fhd.read( 4 ) )[ 0 ] 
  refnames=['']*nc 
  refnameslen=[0]*nc 
  for x in range( nc ): 
    # read each chromosome name 
    nlength = unpack( '<i' , fhd.read( 4 ) )[ 0 ] 
    refstr=fhd.read( nlength ) 
    # jump over chromosome size, we don't need it 
    refstrlen = unpack( '<i', fhd.read( 4 ) )[ 0 ] 
    #fhd.seek( fhd.tell() + 4 ) 
    #print(refstr+':'+str(refstrlen)) 
    refstr=refstr[:-1]
    refnames[x]=refstr 
    refnameslen[x]=refstrlen 
    if refstr not in genedict_sgid:
      logging.warning('sgRNA ID '+ refstr+' in the BAM file does not not match the provided library file. Please double check.')
  logging.info(str(nc)+' references detected in the BAM file.') 
  #
  # next, iterate the bam file
  while True: 
    nline=nline+1
    if nline%1000000==1:
      logging.info('Processing '+str(round(nline/1000000))+ 'M lines..')
    # 
    nreadcount+=1
    tmpdata=fhd.read(4) 
    if len(tmpdata)==0: 
      break 
    entrylength = unpack( '<i', tmpdata )[ 0 ] 
    data = fhd.read( entrylength ) 
    # refid, position 
    refid=unpack( '<i', data[:4] )[ 0 ] 
    if refid == -1: 
      # didn't find any match
      refstr='*' 
    else: 
      # find matches
      refstr=refnames[refid] 
      if refstr not in genedict_sgid:
        logging.warning('sgRNA ID: '+refstr+' is not present in the library file. Please double-check the consistency between library file and SAM/BAM file.')
        continue
      fseqc=genedict_sgid[refstr]
      if fseqc not in ctab:
        ctab[fseqc]=0
      ctab[fseqc]=ctab[fseqc]+1
    # other fields in BAM file; not used
    if False:
      position=unpack( '<i', data[4:8] )[ 0 ] 
      bin_mq_nl=unpack( '<i', data[8:12] )[ 0 ] 
      read_name_len=bin_mq_nl&0x000000ff  
      # print('name length:'+str(read_name_len)) 
      flag_nc=unpack( '<i', data[12:16] )[ 0 ] 
      n_cigar_op=flag_nc&0x0000ffff 
      flag=flag_nc>>16 
      # length 
      readlen = unpack( '<i', data[16:20] )[ 0 ] 
      next_refID=unpack( '<i', data[20:24] )[ 0 ] 
      next_pos=unpack( '<i', data[24:28] )[ 0 ] 
      tlen=unpack( '<i', data[28:32] )[ 0 ] 
      read_name=data[32:32+read_name_len] 
      # cigar 
      tstart=32+read_name_len 
      tend=tstart+n_cigar_op*4 
      # the following is the only 1st int of cigar string 
      cigarstr=unpack('<i',data[tstart:tstart+4]) 
      # seq, encoded 
      tstart=tend 
      tend=tstart+int(math.ceil((readlen+1)/2)) 
      seq=data[tstart:tend] 
      # quality 
      tstart=tend 
      tend=tstart+readlen 
      phredqual=data[tstart:tend] 
      if nline<10 and False: # only for debug purposes 
        logging.info('refid: '+str(refid)+' '+refstr+', position:'+str(position)+', name:'+read_name) 
        # decode the sequence 
        nseqcode=int(math.floor((readlen+1)/2)) 
        seqstring='' 
        for i in range(nseqcode): 
          iit=seq[i] 
          iit_int=unpack('<B',iit)[0] 
          k1=iit_int>>4 
          k2=iit_int&0xf 
          seqstring+=encodetable[k1] 
          if readlen %2==1 and i==nseqcode-1: 
            pass 
          else: 
            seqstring+=encodetable[k2] 
        # end for
        logging.info(seqstring) 
      # end if
    # end if FALSE
  # end while
  fhd.close() 
  
  # calculate statistics
  datastat['reads']=nreadcount
  nmapped=0
  nrdcnt=[]
  for (k,v) in ctab.items():
    if k in genedict:
      nmapped+=v
      nrdcnt+=[math.log(v+1.0)]
  nzerosg=0
  for (k,v) in genedict.items():
    if k not in ctab:
      nzerosg+=1
      nrdcnt+=[math.log(0.0+1.0)]
  logging.info('mapped:'+str(nmapped))
  datastat['mappedreads']=nmapped
  datastat['totalsgrnas']=len(genedict);
  datastat['zerosgrnas']=nzerosg
  datastat['giniindex']=mageckcount_gini(nrdcnt)
  #return ctab
  return 0


def mageckcount_processonefile_sam(filename,args,ctab,genedict,datastat):
  '''
  Go through sam file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A dictionary of sgRNA sequence and count
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})

  Return value
  -----------
  datastat
    a dictionary structure of statistics
  
  # important note for the alignment
  # Please see the notes for bam files

  Reference: 
  BAM specification 
  https://samtools.github.io/hts-specs/SAMv1.pdf 
  Reference: Tao Liu's MACS2 
  https://github.com/taoliu/MACS/blob/master/MACS2/IO/Parser.pyx
  '''
  import sys 
  import gzip 
  import io 
  import math 
  import struct 
 
  nline=0
  logging.info('Parsing SAM file '+filename+'...')
  nreadcount=0
  
  # the {sgRNA_id: sequence} directory
  genedict_sgid={v[0]:k for (k,v) in genedict.items()}
  # whether the warning of particular sgRNA id is already present
  sgrnaid_haswarning={}
  # start processing the sam file
  for line in open(filename):
    if line[0]=='@':
      continue
    nline=nline+1
    if nline%1000000==1:
      logging.info('Processing '+str(round(nline/1000000))+ 'M lines..')
    # 
    nreadcount+=1
    field=line.strip().split()
    # refid, position 
    refid=field[2]
    if refid != "*": 
      # find matches
      refstr= refid
      if refstr not in genedict_sgid:
        if refstr not in sgrnaid_haswarning:
          # logging.warning('sgRNA ID: '+refstr+' in the SAM file does not match the provided library file. Please double-check.')
          sgrnaid_haswarning[refstr]=1
        continue
      fseqc=genedict_sgid[refstr]
      if fseqc not in ctab:
        ctab[fseqc]=0
      ctab[fseqc]=ctab[fseqc]+1
  # end while
  # calculate statistics
  datastat['reads']=nreadcount
  nmapped=0
  nrdcnt=[]
  for (k,v) in ctab.items():
    if k in genedict:
      nmapped+=v
      nrdcnt+=[math.log(v+1.0)]
  nzerosg=0
  for (k,v) in genedict.items():
    if k not in ctab:
      nzerosg+=1
      nrdcnt+=[math.log(0.0+1.0)]
  logging.info('mapped:'+str(nmapped))
  logging.warning('sgRNAs not in the library (may be due to duplicated sequences):'+str(len(sgrnaid_haswarning)))
  datastat['mappedreads']=nmapped
  datastat['totalsgrnas']=len(genedict);
  datastat['zerosgrnas']=nzerosg
  datastat['giniindex']=mageckcount_gini(nrdcnt)
  #return ctab
  return 0



