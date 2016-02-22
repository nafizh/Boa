

"""
1. Gather blasted bacteriocins and context genes
2. Cluster bacteriocins and context genes
3. Run naive bayes to classifiy individual bacteriocins/context genes
4. Run majority vote to classify entire clustesr
"""
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from collections import defaultdict

import sys
import os,shutil
import site
import argparse
import string
import numpy
import re
import subprocess
import pickle

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
print base_path
import genbank
import blast
import intergene
import genome
#import annotated_genes
import annotation #Nafiz
import bacteriocin
import nbayes
import rforests
import cdhit
import fasta
import gff

import hmmer #Nafiz
from mafft import MAFFT #Nafiz
import faa #Nafiz

class PipelineHandler(object):
    def __init__(self,
                 rootdir,
                 genome_dir,
                 intergenes,
                 annotated_genes,
                 bacteriocins,
                 bacteriocin_radius,
                 similarity,
                 bac_evalue,
                 training_labels,
                 training_directory,
                 intermediate,
                 output,
                 numThreads,
                 formatdb,
                 verbose                
                 ):
        
        #Declare global vars
        self.rootdir            =     rootdir  			
        self.genome_dir         =     genome_dir
        self.bacteriocins   	=	  bacteriocins  		
        self.bacteriocin_radius =	  bacteriocin_radius	
        self.similarity	        =	  similarity			
        self.bac_evalue	        =	  bac_evalue
        self.training_labels    =     training_labels			
        self.training_directory =     training_directory			
        self.intermediate   	=	  intermediate  		
        self.output	         	=	  output
        self.numThreads         =     numThreads
        self.formatdb           =     formatdb				
        self.verbose			=	  verbose
        if not os.path.exists(self.intermediate):
            #shutil.rmtree(self.intermediate)
            os.mkdir(self.intermediate)
            
        if intergenes!=None and annotated_genes!=None: 
            self.intergenes	        =	  intergenes			
            self.annotated_genes    =	  annotated_genes   	
        else:
            self.intergenes = "%s/intergeneDB.fa"%self.intermediate
            self.annotated_genes = "%s/annotated_genesDB.fa"%self.intermediate
        
        self.all_fasta          ="%s/all.fna"%self.intermediate
        self.all_faidx          ="%s/all.fai"%self.intermediate
        self.six_fasta          ="%s/all_trans.fna"%self.intermediate
        self.six_faidx          ="%s/all_trans.fai"%self.intermediate
        self.faa                ="%s/all.faa"%self.intermediate
        self.faaidx             ="%s/all.faaidx"%self.intermediate
        self.gff                ="%s/all.gff"%self.intermediate
        self.blasted_tab_bacteriocins = "%s/blasted_bacteriocins.txt"%self.intermediate
        self.blasted_fasta_bacteriocins = "%s/blasted_bacteriocins.fa"%self.intermediate
        self.cand_context_genes_tab = "%s/cand_context_genes.txt"%self.intermediate
        self.cand_context_genes_fasta = "%s/cand_context_genes.fa"%self.intermediate
        self.cand_context_cluster = "%s/cand_context_cluster"%self.intermediate
        self.blast_context_out = "%s/classify"%self.intermediate #def blastContextGenes
        self.hmmer_context_out = "%s/hmmer"%self.intermediate

        # copied from quorum_pipeline.py - Nafiz
        #Declare object handlers
        self.clusterer = None                       
        self.textClassifier = None
        self.hmmers = []   
        self.nbpickle = "nb.zip"
        self.clusterpickle = "cluster.zip"
        self.textout = "text_out.txt"
        self.jobs  = []
        self.split_files = [] 
        self.classes = ["toxin","modifier","immunity","transport","regulator"]
        self.class_files = ["%s/%s.fa"%(self.intermediate,s) for s in self.classes]
        self.hmmer_class_out = ["%s/%s.out"%(self.intermediate,s) for s in self.classes] 
        self.operons_out = "%s/operons.txt"%self.intermediate
        self.pred_operons_out = "%s/predicted_operons.txt"%self.intermediate

        

    """ Builds database such as the 
        intergenic database
        annotated genes database """ #changed - Nafiz
    def preprocess(self,buildAnnotations=True):
        

        # print "Preprocessing"
        # if buildAnnotations:
        #     annotation.go(self.rootdir, self.annotated_genes)    #changed to annotation - Nafiz 
        #     intergene.go(self.rootdir, self.intergenes)


        print "Preprocessing"
        #Combine all genome files into a single genome fasta file
        #https://github.com/mortonjt/Boa/blob/master/src/format/fasta.py
        fasta.go(self.genome_dir,
                 self.all_fasta,
                 self.all_faidx,
                 self.six_fasta,
                 self.six_faidx) 
        indexer = fasta.Indexer(self.all_fasta,self.all_faidx) #a class at fasta.py 
        indexer.index()
        indexer.load()

        #https://github.com/mortonjt/Boa/blob/master/src/genome/intergene.py
        intergene.go(self.genome_dir,self.intergenes)

        #https://github.com/mortonjt/Boa/blob/master/src/annotation/annotation.py
        annotation.go(self.genome_dir,self.annotated_genes,index_obj=indexer) 
        
        #Combine all gff files together
        outhandle = open(self.gff,'w')
        for root, subFolders, files in os.walk(self.genome_dir):
            for fname in files:
                genome_dir = []
                organism,ext = os.path.splitext(os.path.basename(fname))
                absfile=os.path.join(root,fname)
                if ext==".gff":
                    shutil.copyfileobj(open(absfile),outhandle)
        outhandle.close()
        
        tmpfile = "tmp%d.faa"%(os.getpid())
        outhandle = open(tmpfile,'w')
        for root, subFolders, files in os.walk(self.genome_dir):
            for fname in files:
                genome_dir = []
                organism,ext = os.path.splitext(os.path.basename(fname))
                absfile=os.path.join(root,fname)
                if ext==".fna":
                    shutil.copyfileobj(open(absfile),outhandle)
        outhandle.close()
        faa.reformat(tmpfile,self.faa)
        os.remove(tmpfile)
        
        faaindex = fasta.Indexer(self.faa,self.faaidx)
        faaindex.index()


    """ Runs blast to identify bacteriocins and context genes
        
        step 2 and 3 in the paper 
        calls bacteriocin.py to do blast which also calls interval_filter.py 
        to 
        1) Set intersection of species/translation sequences
        2) Set intersection of 50 kb regions
        3) Set intersections of all context genes present
        4) Set intersection of key context genes present
        5) Use e-values for tie breakers
    """
    def blast(self, njobs = 5):
        # print "Blasting"
        # bacteriocin.main(self.genome_dir,
        #                  self.bacteriocins,
        #                  self.intergenes,
        #                  self.annotated_genes,
        #                  open(self.blasted_tab_bacteriocins,'w'),
        #                  open(self.cand_context_genes_tab,'w'),
        #                  self.intermediate,
        #                  self.bac_evalue,
        #                  self.numThreads,
        #                  self.formatdb,
        #                  self.bacteriocin_radius,
        #                  self.verbose,
        #                  False)

        print "Blasting"
        
        """ First split up the main bacteriocin file into a bunch of smaller files"""
        split_bacfiles = ["%s/bacteriocin.%d"%(self.intermediate,i)
                          for i in xrange(njobs)]
        self.split_files += split_bacfiles
        split_bachandles = [open(f,'w') for f in split_bacfiles]

        out_fnames = ["%s/blasted.%d"%(self.intermediate,i) for i in xrange(njobs)]
        out_bac = ["%s.bacteriocins.txt"%(out) for out in out_fnames]
        out_genes = ["%s.annotated.txt"%(out) for out in out_fnames]
        
        index=0
        for record in SeqIO.parse(self.bacteriocins,"fasta"):
            split_bachandles[index].write(">%s\n%s\n"%(str(record.id),
                                                       str(record.seq)))
            index=(index+1)%njobs
        #Close files
        for handle in split_bachandles: 
            handle.close()

        print "evalue",self.bac_evalue

        if self.formatdb: 
            blast_cmd = ' '.join([
                        """module load anaconda; module load blast;module load blast+;"""
                        """python %s/src/genome/bacteriocin.py  """,
                        """ --genome-files %s        """,
                        """ --annotated-genes=%s     """,
                        """ --intergenes=%s          """,
                        """ --bacteriocins=%s        """,
                        """ --bacteriocin-radius=%d  """,
                        """ --bac-evalue=%s          """,
                        """ --num-threads=%d         """,
                        """ --intermediate=%s        """,
                        """ --output=%s              """,
                        """ --formatdb               """,
                        """ --verbose                """
                        ])
        else:
            blast_cmd = ' '.join([
                        """module load anaconda; module load blast;module load blast+;"""
                        """python %s/src/genome/bacteriocin.py  """,
                        """ --genome-files %s        """,
                        """ --annotated-genes=%s     """,
                        """ --intergenes=%s          """,
                        """ --bacteriocins=%s        """,
                        """ --bacteriocin-radius=%d  """,
                        """ --bac-evalue=%s          """,
                        """ --num-threads=%d         """,
                        """ --intermediate=%s        """,
                        """ --output=%s              """,
                        """ --verbose                """
                        ])
            
        """ Release jobs """
        jobs = []
        for i in xrange(njobs):
            cmd = blast_cmd%(self.rootdir,
                            self.all_fasta,
                            self.annotated_genes,
                            self.intergenes,
                            split_bacfiles[i],
                            self.bacteriocin_radius,
                            str(self.bac_evalue),
                            self.numThreads,
                            self.intermediate,
                            out_fnames[i])
            
        #    batch_file = "%s/blast%i.%d.job"%(os.getcwd(),i,os.getpid())
        #    self.batch_files.append(batch_file)
            proc = subprocess.Popen(cmd,shell=True, stdin=subprocess.PIPE)
        #    proc.submit()
            #proc.output = out_fnames[i]
            jobs.append(proc)
            self.jobs.append(proc)
        
        for job in jobs: 
            job.wait() 
        
        """ Collect all of the results from the jobs"""
        bacteriocins_out = open(self.blasted_fasta_bacteriocins,'w')
        context_genes_out = open(self.cand_context_genes_fasta,'w')
        out_bac = ["%s.bacteriocins.txt"%(out) for out in out_fnames]
        out_genes = ["%s.annotated.txt"%(out) for out in out_fnames]

        for i in xrange(njobs):
            if os.path.exists(out_bac[i]):
                shutil.copyfileobj(open(out_bac[i]),bacteriocins_out)
            if os.path.exists(out_genes[i]):
                shutil.copyfileobj(open(out_genes[i]),context_genes_out)

        bacteriocins_out.close()
        context_genes_out.close()


    """ Clusters bacteriocins and context genes together

        input --> "cand_context_genes.fa" file 
        output --> "cand_context_cluster" file 
                   "cand_context_cluster.clstr" file 

    """

    # not being used - Nafiz
    def cluster(self,preprocess=False): # preprocess=False by Nafiz
        print "Clustering"
        if preprocess:
            fasta.preprocess(self.cand_context_genes_tab,
                            self.cand_context_genes_fasta)
            
        
        self.clusterer = cdhit.CDHit(self.cand_context_genes_fasta, 
                                     self.cand_context_cluster,
                                     self.similarity)
        self.clusterer.run()
        self.clusterer.parseClusters()
        self.clusterer.dump(self.clusterpickle)


    # function inserted from quorum_pipeline.py - Nafiz
    """ Identifies context genes using BLAST

        input -> "cand_context_genes.fa" files
        output -> "classify" file

        step 4 in the paper -
        classify context genes into five functional categories.

    """
    def blastContextGenes(self,njobs=2):
        print "Blasting Context Genes"
        """ First split up the main bacteriocin file into a bunch of smaller files"""
        split_fastafiles = ["%s/context.%d"%(self.intermediate,i)
                          for i in xrange(njobs)]
        #self.split_files += split_fastafiles
        split_fastahandles = [open(f,'w') for f in split_fastafiles]
        out_classes = ["%s/contextout.%d"%(self.intermediate,i) for i in xrange(njobs)]
        index=0
        for record in SeqIO.parse(self.cand_context_genes_fasta,"fasta"):
            if len(record.seq)<=1: continue #To weed out weird entries
            split_fastahandles[index].write(">%s\n%s\n"%(str(record.id),
                                                         fasta.format(str(record.seq))))
            index=(index+1)%njobs

        #Close files
        for handle in split_fastahandles: 
            handle.close()
        
        context_cmd = ' '.join([
                                 """module load anaconda; module load blast;module load blast+;""",
                                 """python %s/src/genome/context_gene.py""",
                                 """--training-directory=%s""",
                                 """--training-labels=%s""",
                                 """--query=%s""",
                                 """--intermediate=%s""",
                                 """--num-threads=%d""",
                                 """--output=%s"""         
                                 ])    
        
        """ Release jobs """
        jobs = []
        for i in xrange(njobs):
            cmd = context_cmd%(self.rootdir,
                               self.training_directory,
                               self.training_labels,
                               split_fastafiles[i],
                               self.intermediate,
                               self.numThreads,
                               out_classes[i]
                               )
            
        #    batch_file = "%s/context_blast%i.%d.job"%(os.getcwd(),i,os.getpid())
        #    self.batch_files.append(batch_file)
            
            # changed - Nafiz
            proc = subprocess.Popen( cmd, shell=True, stdin=subprocess.PIPE) 
        #    proc.submit()
            #proc.output = out_classes[i]
            jobs.append(proc)
            self.jobs.append(proc)

        for job in jobs: 
            job.wait()

        """ Collect all of the results from the jobs"""
        context_out = open(self.blast_context_out,'w')

        for i in xrange(njobs):
            if os.path.exists(out_classes[i]):
                shutil.copyfileobj(open(out_classes[i]),context_out)
        context_out.close()    
        pass     



    # Added from quorum_pipeline.py - Nafiz
    """ 
    Creates HMMER profiles of all of the blasted context genes and the 
    bacteriocins to find more candidate context genes and bacteriocins

    input -> "classify" file 
    """
    def hmmerGenes(self,min_cluster=10,njobs=2):
        print "Hmmering Genes"
        # Aggregate all of the bacteriocins and context genes
        # into 5 files representing each class 
        class_handles = [open(c,'w') for c in self.class_files ]
        for record in SeqIO.parse(self.blast_context_out,"fasta"):
            label = record.id.split('|')[-1]
            index = self.classes.index(label)
            SeqIO.write(record,class_handles[index],"fasta")

        #Move all bacteriocins into toxins file
        # copy from blasted_bacteriocins.fa -> toxin.fa
        shutil.copyfileobj(open(self.blasted_fasta_bacteriocins),
                           class_handles[0])

        for c in class_handles: 
            c.close()
        for fname in self.class_files:
            assert os.path.getsize(fname)>0
        #Remove all duplicate ids from this files
        for fname in self.class_files:
            f = "%s.%d"%(fname,os.getpid())
            fasta.remove_duplicates(fname,f)
            os.rename(f,fname)
        for fname in self.class_files:
            assert os.path.getsize(fname)>0
        #self.class_files = self.class_files[::-1] #Reverse list for testing]
        # Run separate instances of HMMER for each class
        print "threads",self.numThreads

        self.hmmers = [hmmer.HMMER(f,self.numThreads,min_cluster) 
                       for f in self.class_files]

        for i in xrange(len(self.class_files)):
            H = self.hmmers[i]
            H.writeClusters(similarity=0.7,memory=3000)
            H.HMMspawn(msa=MAFFT,njobs=njobs)
            H.search(self.faa,self.hmmer_class_out[i],maxpower=True,njobs=njobs)
    
    # Added from quorum_pipeline.py
    """ Writes clusters and their corresponding sequences"""
    def writeClusters(self,clusters,seq_dict,outhandle):
        clusternum=0
        for cluster in clusters:
            for gene in cluster:
                acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description,strand,protid=gene.split("|")
                hmm_st,hmm_end,env_st,env_end = map(int,[hmm_st,hmm_end,env_st,env_end])
                full_evalue = float(full_evalue)
                function = clrname.split('.')[0]
                seq = seq_dict[(acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description,strand,protid)]
                outhandle.write(">accession=%s|function=%s|start=%s|end=%s|strand=%s|score=%s|protein_id=%s|cluster_%d|%s\n%s\n"%
                                (acc,function,env_st,env_end,strand,str(full_evalue),protid,clusternum,description,fasta.format(seq))) 
            clusternum+=1



    # copied from quorum_pipeline.py - Nafiz

    """ Finds operons by constructing graphs and finding cliques 
    TODO: Move these parameters to main pipeline handler object"""
    def cliqueFilter(self,clique_radius=50000,threshold=62,functions = ["toxin","modifier","immunity","transport","regulator"]):
        print "Clique filtering","Looking for cliques with",functions
        
        toxin_hits     = hmmer.parse("%s/toxin.out"%self.intermediate)
        modifier_hits  = hmmer.parse("%s/modifier.out"%self.intermediate)
        immunity_hits  = hmmer.parse("%s/immunity.out"%self.intermediate)
        regulator_hits = hmmer.parse("%s/regulator.out"%self.intermediate)
        transport_hits = hmmer.parse("%s/transport.out"%self.intermediate)
        faaindex = fasta.Indexer(self.faa,self.faaidx)
        faaindex.index()
        faaindex.load()
        genefile = gff.GFF(self.gff,fasta_index=self.faaidx)
        genefile.indexdb()
        toxin_hits     = genefile.call_orfs(toxin_hits    ,faaindex)
        modifier_hits  = genefile.call_orfs(modifier_hits ,faaindex)
        immunity_hits  = genefile.call_orfs(immunity_hits ,faaindex)
        regulator_hits = genefile.call_orfs(regulator_hits,faaindex)
        transport_hits = genefile.call_orfs(transport_hits,faaindex)

        toxin_hits = threshold_filter.filter(toxin_hits,threshold)
        modifier_hits = threshold_filter.filter(modifier_hits,threshold)
        immunity_hits = threshold_filter.filter(immunity_hits,threshold)
        regulator_hits = threshold_filter.filter(regulator_hits,threshold)
        transport_hits = threshold_filter.filter(transport_hits,threshold)

        all_hits = toxin_hits+modifier_hits+immunity_hits+regulator_hits+transport_hits
        seq_dict = {x[0]:x[1] for x in all_hits}
   
        del all_hits
        
        toxin_ids,toxin_seqs = zip(*toxin_hits)
        modifier_ids,modifier_seqs = zip(*modifier_hits)
        immunity_ids,immunity_seqs = zip(*immunity_hits)
        regulator_ids,regulator_seqs = zip(*regulator_hits)
        transport_ids,transport_seqs = zip(*transport_hits)
        
        del toxin_hits
        del modifier_hits
        del immunity_hits
        del regulator_hits
        del transport_hits
        
        all_ids = toxin_ids+modifier_ids+immunity_ids+regulator_ids+transport_ids
        all_ids = interval_filter.unique(all_ids)
        # #Sort by start/end position and genome name
        all_ids=sorted(all_ids,key=lambda x: x[6])   
        all_ids=sorted(all_ids,key=lambda x: x[5])
        all_ids=sorted(all_ids,key=lambda x: x[0])
        #all_ids=sorted(all_ids,key=lambda x: x[-3])
        
        del toxin_ids
        del modifier_ids
        del immunity_ids
        del regulator_ids
        del transport_ids
        print "all ids",len(all_ids)
        print '\n'.join(map(str,all_ids[:10]))
        #Find operons with at least a toxin and a transport
        clusters = clique_filter.findContextGeneClusters(all_ids,
                                                         radius=clique_radius,
                                                         backtrans=False,
                                                         functions=["toxin","transport"])
        print "Clusters: ",len(clusters)
        outhandle = open(self.operons_out,'w')
        self.writeClusters(clusters,seq_dict,outhandle) 
        outhandle.close()
        #Predict operons based on just context genes
        clusters = clique_filter.findContextGeneClusters(all_ids,
                                                         radius=clique_radius,
                                                         backtrans=False,   
                                                         functions=["modifier","regulator","immunity","transport"])
        print "Clusters: ",len(clusters)
        outhandle = open(self.pred_operons_out,'w')
        self.writeClusters(clusters,seq_dict,outhandle) 
        outhandle.close()



    """ Classifies individual bacteriocins and context genes based on their text
        not used - Nafiz
    """
    def naiveBayes(self,db):
        self.textClassifier.train()
        sets = self.textClassifier.classify(db,self.cand_context_genes_fasta)
        titles,labels = zip(*sets)
        
        open(self.textout,'w').write('\n'.join(["%s\t%s"%x for x in sets]))
        

    
    """ Classifies entire clusters based on a majority vote"""
    def majorityVote(self):
        pass
    def main(self):
        pass

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins and context genes')
    parser.add_argument(\
        '--pipeline-section', type=str, required=False, default="all",
        help='Section of the pipeline to run (all, preprocess, blast, cluster, nbayes, mvote)')
    parser.add_argument(\
        '--root-dir',type=str, required=False,
        help='Root directory')
    parser.add_argument(\
        '--genome-dir', type=str, required=False,
        help='Directory of all of the fasta/genbank files')
    parser.add_argument(\
        '--intergenes', type=str, required=False,default=None,
        help='FASTA files containing intergenic regions')
    parser.add_argument(\
        '--annotated-genes', type=str, required=False,default=None,
        help='FASTA files containing annotated genetic regions')
    parser.add_argument(\
        '--bacteriocins', type=str, required=False,default=None,
        help='The bacteriocin proteins that are to be blasted')
    parser.add_argument(\
        '--bacteriocin-radius', type=int, required=False, default=5000,
        help='The search radius around every specified bacteriocin')
    parser.add_argument(\
        '--similarity', type=int, required=False, default=0.7,
        help='Clustering similarity')    
    parser.add_argument(\
        '--bac-evalue', type=float, required=False, default=0.00001,
        help='The evalue for bacteriocin hits')
    parser.add_argument(\
        '--training-labels',type=str,required=False,default=None,
        help='Training labels for Naive Bayes')
    parser.add_argument(\
        '--training-directory',type=str,required=False,default=None,
        help='''Training directory containing all 
                genbank files required for training the Naive Bayes model''')
    parser.add_argument(\
        '--intermediate', type=str, required=False,default='.',
        help='Directory for storing intermediate files')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='The output file basename for filtered annotationed regions and bacteriocins')
    parser.add_argument(\
        '--formatdb', action='store_const', const=True, default=False,
        help='Indicates if formatdb should be run')
    parser.add_argument(\
        '--num-threads', type=int, required=False, default=1,
        help='The number of threads to be run by BLAST')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Messages for debugging')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    
    if not args.test:
        proc = PipelineHandler(  args.root_dir,
                                 args.genome_dir,
                                 args.intergenes,
                                 args.annotated_genes,
                                 args.bacteriocins,
                                 args.bacteriocin_radius,
                                 args.similarity,
                                 args.bac_evalue,
                                 args.training_labels,
                                 args.training_directory,
                                 args.intermediate,
                                 args.output,
                                 args.num_threads,
                                 args.formatdb,
                                 args.verbose                
                                ) 
    
        proc.preprocess()
        proc.blast()
    #    proc.cluster()
        proc.blastContextGenes()
        proc.hmmerGenes()
    #    proc.cliqueFilter()
    #    proc.naiveBayes()

        
    else:
        del sys.argv[1:]
        import unittest
        import test_modules
        class TestRun(unittest.TestCase):
            def setUp(self):
                self.root = "/home/nafizh/Boa_project/Boa"
                self.exampledir = "%s/example/Streptococcus_pyogenes"%self.root
                self.bacdir = "%s/bacteriocins"%self.root
                self.genome_dir = test_modules.getFNA(self.exampledir)
                self.bacteriocins = "%s/bagel.fa"%self.bacdir
                self.intergenes = "%s/src/genome/test_intergenes.fa"%self.root #changed - Nafiz
                self.annotated_genes = "%s/src/annotation/test_genes.fa"%self.root #changed - Nafiz
                self.intermediate = "intermediate"
                #self.training_labels = "%s/data/training/training.txt"%self.root
                self.training_directory = "%s/data/training/protein"%self.root
                self.training_labels = "%s/data/training/training_proteins.txt"%self.root
                if not os.path.exists(self.intermediate):
                    os.mkdir(self.intermediate)
                self.bac_evalue = 0.000001
                self.formatdb = True
                self.bacteriocin_radius = 50000
                self.verbose = True
                self.keep_tmp = False
                self.similarity = 0.65
                self.numThreads = 1
                self.output = "out"
                self.keep_tmp = True
            def testrun(self):
                print "Test Run"
                proc = PipelineHandler(  self.root,
                                         self.genome_dir,
                                         self.intergenes,
                                         self.annotated_genes,
                                         self.bacteriocins,
                                         self.bacteriocin_radius,
                                         self.similarity,
                                         self.bac_evalue,
                                         self.training_labels,
                                         self.training_directory,
                                         self.intermediate,
                                         self.output,
                                         self.numThreads,
                                         self.formatdb,
                                         self.verbose                
                                        ) 
                
            #    proc.preprocess()
                self.assertTrue(os.path.getsize(self.annotated_genes) > 0)
                self.assertTrue(os.path.getsize(self.intergenes) > 0)
                
                #proc.blast()
                
                self.assertTrue(os.path.getsize(proc.blasted_tab_bacteriocins) > 0)
                self.assertTrue(os.path.getsize(proc.cand_context_genes_tab) > 0)       
                
                #proc.cluster()
                
                self.assertTrue(os.path.getsize("%s"%(proc.cand_context_genes_fasta)) > 0)                
                self.assertTrue(os.path.getsize(proc.cand_context_cluster) > 0)
                
                proc.load()
                proc.naiveBayes()
                self.assertTrue(os.path.getsize(proc.textout) > 0)
                
        unittest.main()


