##How to run Boa:

First, you need to install the following things:

* Python2.7 or greater
* biopython
* matplotlib
* numpy
* panda
* bx-python
* blastall
* hmmer
* cdhit
* clustalw
* nltk
* mafft
* transeq

Most of these can be installed with the following command on Ubuntu. You can run it with shift+enter, it will run on your terminal.

```sh
sudo apt-get install python python-biopython python-matplotlib python-panda python-numpy nltk clustalw cd-hit hmmer
```

You need to also install blastall, mafft, transeq and bx-python.

blast can be downloaded from here- 
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

mafft can be downloaded from here-
http://mafft.cbrc.jp/alignment/software/source.html

For transeq, you need to install emboss, which can be downloaded from here (you might need to also install x11 dev
files for this, it will show error while installing if you need this)-
http://emboss.open-bio.org/html/adm/ch01s01.html

bx-python can be installed through the following link: 
https://bitbucket.org/james_taylor/bx-python/wiki/Home

If everything is installed properly, now you can run Boa with the pipeline.py file(inside src folder). Following is the command that you can run from your terminal.

```sh
    python pipeline.py 
    --root-dir="/home/nafizh/Boa_project/Boa" 
    --genome-dir="/home/nafizh/Boa_project/Boa/boa_genbank_data/" 
    --bacteriocins="/home/nafizh/Boa_project/Boa/bacteriocins/bagel.fa" 
    --intermediate="/home/nafizh/Boa_project/Boa/intermediate_me" 
    --training-labels="/home/nafizh/Boa_project/Boa/data/training/training_proteins.txt" 
    --training-directory="/home/nafizh/Boa_project/Boa/data/training/protein" 
    --formatdb 
```

Here

* --root-dir is the folder inside which all code and all experimental data reside. My root directory here is 
  "/home/nafizh/Boa_project/Boa", all my code and data is inside this Boa folder.
* --genome-dir is the folder where your bacterial genome files reside.
* --bacteriocin is the file(/bacteriocins/bagel.fa) with bacteriocins that are identified by BAGEL3. It is one part of the "gold standard" data set.
* --intermediate is a directory of your choice where all the files during the process and the result files will be generated.
* --training-labels("data/training/training_proteins.txt") is the literature curated dataset.
* --traning-directory is the directory that contains all of the genbank files for each curated protein.
* --formatdb is required to form a database from experimental bacterial genome files inside --genome-dir. This databased will be run against the --bacteriocin file with blast. This is only required once. So, once run, if you run pipeline.py again, you can ignore it unless your experimental bacterial genome files have changed.

### How pipeline.py is working:

I will try to give a brief overview of what pipeline.py is doing under the hood. It has 5 distinct steps inside the
```python
    if __name__ == __main__:
```
function. They are:

+ preprocess
+ blast
+ blastContextGenes
+ hmmerGenes
+ cliqueFilter


1. During the **preprocess** step, these files will be generated using the files inside --genome-dir:
   all.faa, all.faaidx, all.fai, all.fna, all.gff, all_trans.fai, all_trans.fna, annotated_genesDB.fa, 
   intergeneDB.fa.
   
   For **annotated_genesDB.fa**, preprocess is taking the "amino acid translation corresponding to the nucleotide 
   coding  sequence(CDS)" from each .gbk(genbank - http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) file and 
   putting it inside **annotated_genesDB.fa** with the filename, base span, whether it is in complementary strand 
   or not, locus tag and protein id.
   
   The preprocess step also finds the integenic regions and corresponding nucleotides from the .fna(fasta) files
   and puts it inside the **intergeneDB.fa** file.
   
   The other files are result of necessary formatting procedures for future steps.
   
   preprocess function is using the **annotation.py** and **intergene.py** files to do these.
   <br><br>
2. The **blast** step is creating a database from the all.fna file with a command -
   **"formatdb -i /home/nafizh/temp/all.fna -p F"**
   
   It is then running that database against bagel.fa file with blast(tblastn command specifically). blast function 
   inside pipeline.py is using **bacteriocin.py** to do all these.
   
   In this step, files with names like **"blasted.0.bacteriocins.txt"** and **"blasted.0.annotated.txt"** are
   produced.
   
   The **"blasted.0.bacteriocins.txt"** type files are merged into **blasted_bacteriocins.fa** file. It contains 
   the list of bacteriocins aligned against all of the bacterial genomes provided. It contains the following 
   things(shown with an example from the file)-
   
   > \>80.3|gi|256632183|dbj|AP011121.1|2340005|2340748|-|gene
       MNNLHKHLAPISHAAWAEIEQEASRTIRRNLAGRRVVDTPEPKGTAFSSVGTGRNKQI--QSPSDGIQAV---- 
       QREVLPVIELRVPFTLSRAEIDAVERGSLDSDWQPVKDAAQKIAFAEDRAIFDGYTAAGITGIRQGSSNPHTKLPTSAKDYPRAITNALDTLRLAGVNGPYALV
       LGTKAY-QAVSGGDDVGYPVLKHIESLIEGEPIWAPAIEGAFVISKRGGDLQLDIGQDFSIGYLSHTAETVELYLQE
      
      Here, 
            80.3 - bacteriocin name
            gi|256632183|dbj|AP011121.1 - ncbi id of species blasted against(from .fna file of the species)
            2340005 - bacteriocin start
            2340748 - bacteriocin end
            '-' -  bacteriocin strand
            gene - overlaps intergene or gene
            MNNLHKHLAPISHAAWAEIEQEASRTIRRNLAGR....... - blasted bacteriocin sequence
            
   The **"blasted.0.annotated.txt"** type files are merged into **cand_context_genes.fa** file. It contains the
   list of annotated genes within a radius around all of the blasted bacteriocins. This search radius can be 
   specified in the bacteriocin.py script. It contains the following things(shown with an example from the file)-
   
   > \>80.3|gi|256632183|dbj|AP011121.1|2340005|2340748|-|AP011121|APA01_21470|BAI00262.1|2334839|2335451|+
      MVAGCAAFHPLGEARRMGAAARTGSGTTGSGIRNDFSGWHKHQGSPQSGGSPKKGASFEERDHREALGRSRGGYGTKVCVIADGHGKAFGFALAPGQAHELPLAPAM
      LDSLPATPLWVVADKGYASNAMRERIWDMGARPAIPAKRRDGPVACPKWAYRCRHLVENLWARLKEWRAVATRYEKTATSFLAVIHIAAAADWIKP
      
      Here,
            80.3 - bacteriocin name
            gi|256632183|dbj|AP011121.1 - ncbi id of anchor gene
            2340005 - blast bacteriocin start
            2340748 - blast bacteriocin end
            '-' - blast bacteriocin strand
            AP011121|APA01_21470|BAI00262.1 - accession id of whole genome
            2334839 - anchor gene start
            2335451 - anchor gene end
            '+' - anchor gene strand
            MVAGCAAFHPLGEARRMGAAARTGSG....  - sequence of bacteriocin
   <br><br>         
3. The **blastContextGenes** step identifies context genes using BLAST.
   
   It first splits up the **"cand_context_genes.fa"** file to a bunch of **context.0** type of files depending on 
   the number of njobs provided to the function. Then thr results are put into contextout.0 type files and then put 
   together into the classify file. So, basically for the blastContextGenes() function-
       input: cand_context_genes.fa -> conext.0, context.1 ...
       output: contextout.0, contextout.1 ... -> classify
   <br><br>    
4. The **hmmerGenes** step creates HMMER profiles of all of the blasted context genes and the bacteriocins to find 
more candidate context genes and bacteriocins.

   It is taking the **"classify"** file and putting the context genes in "modifier.fa", "immunity.fa",   
   "transport.fa", "regulator.fa" files.
   
   Then It is taking the **"blasted_bacteriocins.fa"** file and putting the toxins into the **"toxin.fa"** file.
   
   For each type of context genes and bacteriocin, if the number of clusters is equal to or more than the set
   **"min_cluster"** variable, then **"toxin.fa.cluster0.fa"** type files are created with the following line.
   
   ```python
    H.writeClusters(similarity=0.7,memory=3000)
   ```
   
   Then with the following code, it is taking each **"toxin.fa.cluster0.fa"** type file, and first, performing
   multiple alignment with Muscle on each cluster with MAFFT. Then builds an HMM for each cluster. In this step,
   the **"toxin.fa.cluster0.fa.hmm"**, **"toxin.fa.cluster0.fa.sto"** and **"toxin.fa.cluster0.fa_hmmbuild.log"**
   type files are created.

    
   ```python
    H.HMMspawn(msa=MAFFT,njobs=njobs)

    #this function is within the hmmer.py file
    def HMMspawn(self,msa=MAFFT,njobs=4,maxiters=20):
        procs = []
        i = 0

        # number of transport.fa.cluster0.fa type files for each functional category of context genes
        print len(self.clusterfas)

        #looping over the transport.fa.cluster0.fa type files 
        for clrfa in self.clusterfas:
            hmm = HMM(clrfa,module=self.module,threads=self.threads)
            self.hmms.append(hmm)
        for hmm in self.hmms:
            proc = hmm.multipleAlignment(msa=msa,maxiters=maxiters)
            if proc!=None:
                procs.append(proc)
            i+=1
            if i==njobs: #make sure jobs don't overload
                self.wait(procs,fasta=True)
                procs,i = [],0
        self.wait(procs,fasta=True)

        procs,i = [],0
        for hmm in self.hmms:
            proc = hmm.hmmbuild()
            procs.append(proc)
            i+=1
            if i==njobs: #make sure jobs don't overload
                self.wait(procs)
                procs,i = [],0
        self.wait(procs)

    #this function is within the hmmer.py file
    """Perform multiple alignment with Muscle on a single cluster"""
    def multipleAlignment(self,module=subprocess,msa=MAFFT,maxiters=20):
        if self.seqcount(self.clusterfa)==1:
            self.sto = self.clusterfa
            return None
        else:
            #mafft class in mafft.py is being called - Nafiz
            cw = msa(self.clusterfa,self.sto,module)
            self.clustal = cw
            if msa==Muscle:
                proc = cw.run(fasta=True,maxiters=maxiters)
            else:
                proc = cw.run(fasta=True,maxiters=maxiters,threads=self.threads)
            #cw.outputSTO()
            return cw

    #this function is within the hmmer.py file        
    """Builds an HMM for a cluster"""
    def hmmbuild(self,module=subprocess):
        cmd = "hmmbuild %s %s "%(self.hmm,self.sto)
        #Changed from self.module.Popen -> module.Popen - Nafiz
        proc = subprocess.Popen(cmd, stderr=open(self.logs[0],'w+'),shell=True)
        if self.module==quorum: proc.submit()
        #proc.wait()
        return proc
   ``` 
   
   Then with the following code, it is taking the **all.faa** file and the **\*.out** type files, and performs
   HMMER using all clusters on the all.faa file. The **"toxin.fa.cluster0.fa.table"** type files are being
   created here also. It is writing the results in the **\*.out** files.
   
    
   ```python
    H.search(self.faa,self.hmmer_class_out[i],maxpower=True,njobs=njobs)

    #this function is in the hmmer.py file
    """Performs HMMER using all clusters on infasta"""
    def search(self,infasta,out,maxpower=False,njobs=4):
        procs = []
        i = 0
        for hmm in self.hmms:
            proc = hmm.hmmsearch(infasta,maxpower)
            self.tables.append(hmm.table)
            procs.append(proc)
            i+=1
            if i==njobs: #make sure jobs don't overload
                self.wait(procs)
                procs,i = [],0
        self.wait(procs)
        procs,i = [],0     

        with open(out, 'w') as outfile:
            for fname in self.tables:
                if os.path.exists(fname):
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
    
    #this function is in the hmmer.py file
    """Runs Viterbi algorithm on an input protein fasta"""
    def hmmsearch(self,infasta,maxpower=False):
        if maxpower:
            cmd = "hmmsearch --noali --notextw --max --domtblout %s %s %s"%(self.table,self.hmm,infasta)
        else:
            cmd = "hmmsearch --noali --notextw --domtblout %s %s %s"%(self.table,self.hmm,infasta)
        print cmd
        # changed - Nafiz
        proc = subprocess.Popen(cmd,stderr=open(self.logs[0],'w+'),shell=True)
    #    if self.module==quorum: proc.submit() #Nafiz
        return proc
   ```
 <br></br>
5.The **cliqueFilter** step finds operons by constructing graphs and finding cliques. It takes the \*.out files as input and outputs the results in the **"operons.txt"** and **"predicted_operons.txt"** files.


