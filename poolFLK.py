#!/usr/bin/python

import numpy as np
import heapq
import argparse
import sys
import logging
from hapflk import popgen

##### GLOBAL VARIABLES #####

VERSION = "0.1"
logging.basicConfig(format='[%(levelname)s]\t%(asctime)s - %(message)s', datefmt='%m/%d/%Y %H:%M:%S')
LOGGER = logging.getLogger("logger")

##### Program options ####
def getCommandLineParser():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # some options
    parser.add_argument('--version', help='print the version and exits', default=False, action='store_true')
    parser.add_argument('-p','--prefix',dest='prefix',help='prefix for output files',default='hapflk')
    parser.add_argument('--only_kinship', help='Perform only the Kinship Matrix estimation', default=False, action='store_true')
    parser.add_argument('--debug',help=argparse.SUPPRESS,default=False,action="store_true") ## for debug purpose
    # input options
    io_opts=parser.add_argument_group('Input Files')
    io_opts.add_argument('--sync',metavar='FILE.sync',help='Synchronized population file (gzipped or not)')
    io_opts.add_argument('--pop_names',metavar='FILE',help='Population names in one column for the synchronized file.', default=None)
    io_opts.add_argument('--kinship',help='Read population kinship from file (if None, kinship is estimated)',metavar='FILE',default=None)
    # flk options
    flk_opts=parser.add_argument_group('Population kinship ','Set parameters for getting the population kinship matrix')
    flk_opts.add_argument('--reynolds-snps',dest='reysnps',type=int,help='Number of SNPs to use to estimate Reynolds distances',default=10000,metavar='L')
    flk_opts.add_argument('--min-freq', dest='minfreq', type=float, help='Minimum allele frequency of SNPs to use to estimate Reynolds distances', default=0.01,metavar='F')
    flk_opts.add_argument('--outgroup',default=None,help='Use population POP as outgroup for tree rooting (if None, use midpoint rooting)',metavar="POP")
    flk_opts.add_argument('--keep-outgroup',dest='keepOG',default=False,help='Keep outgroup in population set',action="store_true")
    pop_group=parser.add_argument_group('SNP/Population selection','Filter SNP/populations')
    pop_group.add_argument('--pops',help='Select by index (1-based) populations in the sync file. Comma-separated list of populations',metavar='POPLIST')
    return parser

##### Simple Sync Parser #####

class SinglePop:

    ## no Ns or deletions
    def __init__(self,a,t,c,g,n,deletion):
        self.__a=int(a)
        self.__t=int(t)
        self.__c=int(c)
        self.__g=int(g)
        
    def issnp(self,mincount):
        alcount=self.count_alleles(mincount)
        if(alcount>1):
            return True
        else:
            return False

    def count_alleles(self,mincount):
        alcount=0
        if self.A>=mincount:
            alcount+=1
        if self.T>=mincount:
            alcount+=1
        if self.C>=mincount:
            alcount+=1
        if self.G>=mincount:
            alcount+=1
        return alcount

    def countForAllele(self,allele):
        return eval("self."+allele)     

    @property
    def A(self):
        return self.__a
    
    @property
    def T(self):
        return self.__t

    @property
    def C(self):
        return self.__c

    @property
    def G(self):
        return self.__g

    @property
    def cov(self):
        return self.A+self.T+self.C+self.G

    def __str__(self):
        return ":".join(map(str,[self.A,self.T,self.C,self.G,0,0]))

class PopLine:
    def __init__(self,chr,pos,refc,populations,pvalue=None):
        self.__chr=chr
        self.__pos=int(pos)
        self.__refc=refc
        self.__populations=populations
        self.__pvalue=pvalue

    @property
    def chr(self):
        """
        >>> PopLine("2L",1,"N",[],0.2).chr
        '2L'
        """
        return self.__chr
    
    @property
    def pos(self):
        return self.__pos
    
    @property
    def refc(self):
        return self.__refc

    @property
    def populations(self):
        return self.__populations
    
    def subpopulations(self,populations):
        tpops=self.populations
        toret=[]
        for i in populations:
            toret.append(tpops[i-1])
        return toret
    
    @property
    def popcount(self):
        return len(self.__populations)
    
    def __str__(self):
        popstr="\t".join(map(str,self.populations))
        tojoin=[self.chr,self.pos,self.refc,popstr]
        return "\t".join([str(x) for x in tojoin])

class SyncReader:
    def __init__(self,fileName):
        self.__filename=fileName
        self.__openFH()

    def __iter__(self):
        return self
    
    def next(self):
        line=""
        while(1):
            line=self.__filehandle.readline()
            if line=="":
                raise StopIteration
            line=line.rstrip('\n')
            if line != "":
                break
        
        a=line.split()
        chr=a.pop(0)
        pos=a.pop(0)
        refc=a.pop(0)
        population=[]
        for p in a:
            po=None
            if p=="-":
                po=SinglePop(0,0,0,0,0,0)
            else:
                s=p.split(":")
                po=SinglePop(*s)
            population.append(po)
        
        return PopLine(chr,pos,refc,population)

    def countSnps(self):
        total = 0
        for line in self.__filehandle:
            total += 1
        self.close()
        self.__openFH()
        return total

    def close(self):
        self.__filehandle.close()

    def __openFH(self):
        if(isinstance(self.__filename,str)):
            if self.__isGZIP(self.__filename):
                self.__filehandle=gzip.open(self.__filename,"r")
            else:
                self.__filehandle=open(self.__filename,"r")
        else:
            self.__filehandle=self.__filename

    def __isGZIP(self, fileName):
        return fileName.endswith("gz") or fileName.endswith("bgz")

class SNP():
    
    def __init__(self, chrom, pos, al1, al2):
        self.__chr = chrom
        self.__pos = pos
        self.__al1 = al1
        self.__al2 = al2

    @property
    def chr(self):
        return str(self.__chr)

    @property
    def pos(self):
        return self.__pos

    @property
    def al1(self):
        return str(self.__al1)

    @property
    def al2(self):
        return str(self.__al2)

##### Script Methods #######

def getPopulationInfo(fileName, popList, popNamesFile = None):
    LOGGER.info("Loading population information")
    syncReader = SyncReader(fileName)
    npops = syncReader.next().popcount
    syncReader.close()
    if popList is None:
        populations = xrange(1, npops+1)
    else:
        populations = list(sorted(set(map(int, popList.split(",")))))
        if 0 in populations:
            LOGGER.error("Population list should be 1-based and found a 0")
            sys.exit(1)
            raise
        elif populations[-1] > npops:
            LOGGER.error("Population list should could not contain more than the number of populations in the file")
            sys.exit(1)
    
    if popNamesFile:
        popNames = []
        with open(popNamesFile) as f:
            line = f.readline().rstrip('\n')
            while line != "":
                popNames.append(line)
                line = f.readline().rstrip('\n')
        if len(popNames) != npops:
            LOGGER.error("Population names file should have the same rows as population columns")
            LOGGER.debug("popNames =  %s | npops = %s", popNames, npops)
            sys.exit(1)
        popNames = [popNames[i-1] for i in populations]
    else:
        popNames = map(str, populations)
    return populations, popNames


def getFreqMatrix(fileName, populations, min_freq):
    LOGGER.info("Loading allele frequencies for %s populations", len(populations))
    reader = SyncReader(fileName)
    sync_bases = ['A', 'T', 'C', 'G']
    myMap = []
    frqs = []
    snpIndx = []
    indx = 0
    # for each record
    for record in reader:
        totalCounts = np.zeros(shape=(len(populations),len(sync_bases)),dtype='int16')
        i = 0
        for pop in record.subpopulations(populations):
            popCounts = []
            for base in sync_bases:
                popCounts.append(pop.countForAllele(base))
            totalCounts[i,] = popCounts
            i += 1
        totalSum = np.sum(totalCounts, 0)
        # two major alleles indexes
        # indexes = np.argpartition(totalSum, 1)[len(sync_bases)-2:]
        indexes = map(totalSum.tolist().index, heapq.nlargest(2, totalSum)[::-1])
        ## get biallelic
        biallelic = totalCounts[:,indexes]
        ## get the alleles
        alleles = [sync_bases[x] for x in indexes]
        alleleFreqs = biallelic.astype(float) / np.sum(biallelic,1)[:,None]
        frqs.append(alleleFreqs[:,0])
        myMap.append(SNP(record.chr, record.pos, alleles[0], alleles[1]))
        # filter the ones that could be used for reynolds distance
        allFreq = totalSum[indexes]/np.sum(a[indexes])[0]
        if allFreq <= min_freq:
            snpIndx.append(indx)
        indx += 1
    return np.transpose(np.vstack(frqs)), myMap, snpIndx

def computeAndWriteKinship(freqs, nsnp, popNames, keepOG, filePrefix, myMap, snpIdx):
    LOGGER.info("Computing Reynolds distances and heterozygosity")
    snp_subset = np.asarray(snpIdx)
    if nsnp < len(snpIdx):
        ## TODO: be careful with the indexes monomorfic SNPs (filter them before)
        pbinom = float(nsnp)/len(snpIdx)
        subset = np.array(np.random.binomial(1,pbinom,len(snpIdx)),dtype=bool)
        snp_subset = snp_subset[subset]
    reynolds_dist = popgen.reynolds(freqs[:,snp_subset])
    heteroZ = popgen.heterozygosity(freqs[:,snp_subset])
    writeFreqs(freqs[:,snp_subset], [myMap[i] for i in snp_subset], popNames, filePrefix+"-reynolds")
    
    LOGGER.info("Computing Kinship Matrix")
    fij = popgen.popKinship_new(reynolds_dist, popNames, outgroup = options.outgroup, fprefix=filePrefix, keep_outgroup = keepOG, hzy = heteroZ)
    filename = filePrefix + "_fij2.txt"
    fout=open(filename,'w')
    for i in range(fij.shape[0]):
        tw=[popNames[i]]
        for j in range(fij.shape[1]):
            tw.append(str(fij[i,j]))
        print >> fout,' '.join(tw)
    fout.close()
    return filename
    
def writeFreqs(freqs, myMap, popNames, filePrefix):
    fout=open(filePrefix+".frq",'w')
    print >> fout,'chr','pos','all_min','all_maj',' '.join(popNames)
    for indx in xrange(len(myMap)):
        snp = myMap[indx]
        tw = [snp.chr, str(snp.pos), snp.al1, snp.al2]
        for ip, nom in enumerate(popNames):
            tw.append(str(freqs[ip, indx]))
        print >> fout,' '.join(tw)
    fout.close()

def computeFLK(kinship, freqs, myMap, filePrefix):
    ## computing
    LOGGER.info("Computing FLK")
    myFLK=popgen.FLK_test(kinship)
    # columns in res = 'pzero','FLK','pval.FLK','eigen.FLK'
    myFLK_res = np.apply_along_axis(myFLK.eval_flk,0,freqs)
    LOGGER.debug("Writing results")
    fout=open(filePrefix+".flk",'w')
    ## TODO: remove pzero
    print >>fout,'chr','pos','pzero','flk','pvalue'
    for indx in xrange(len(myMap)):
        snp = myMap[indx]
        tw = [snp.chr, snp.pos, myFLK_res[0,indx], myFLK_res[1,indx], myFLK_res[2,indx]]
        print >> fout,' '.join([str(x) for x in tw])

## main
if __name__=='__main__':

    myparser = getCommandLineParser()
    options=myparser.parse_args()
    if options.version:
        print VERSION
        sys.exit(0)
    if not options.sync:
        myparser.print_help()
        sys.exit(1)
    if options.debug:
        LOGGER.setLevel(logging.DEBUG)
        LOGGER.debug("debug mode")
    else:
        LOGGER.setLevel(logging.INFO)

    LOGGER.info("Starting")
    populations, popNames = getPopulationInfo(options.sync, options.pops, options.pop_names)
    freqMatrix, myMap, snpIndx = getFreqMatrix(options.sync, populations, options.minfreq)
    ## TODO: check if it is correct
    LOGGER.info("Loaded %s SNPs (%s polymorphic)", freqMatrix.shape[1], len(snpIndx))
    ## get hte kinship matrix
    if not options.kinship:
        options.kinship = computeAndWriteKinship(freqMatrix, options.reysnps, popNames, options.keepOG, options.prefix, myMap, snpIndx)

    LOGGER.info("Loading Kinship Matrix")
    ## decide if keep the outgroup
    if not options.keepOG:
        kinship=popgen.popKinship_fromFile(options.kinship,[pop for pop in popNames if pop != options.outgroup])
    else:
        kinship=popgen.popKinship_fromFile(options.kinship, popNames)

    if not options.only_kinship:
        if options.keepOG:
            filter_outgroup = np.array([x != options.outgroup for x in popNames],dtype=bool)
            computeFLK(freqMatrix[filter_outgroup,])
        else:
            computeFLK(kinship, freqMatrix, myMap, options.prefix)
    LOGGER.info("Finished")