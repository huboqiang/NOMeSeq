A pipeline which could processing from raw fastq reads to result of NORM-seq for both DNA methylation(WCG) and the accessibility of DNA (GCH).

First, before this pipeline in a server, make sure the required modules were installed. If not, running the following scripts for deploying.

Remember, DO USE the **right version** of the software listed below, or some bugs would be introduced.

```bash
mkdir install_packages

### install python anaconda 2.2.0
cd software/
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.2.0-Linux-x86_64.sh
bash Anaconda-2.2.0-Linux-x86_64.sh  # prefix=/path/for/anaconda
mv Anaconda-2.2.0-Linux-x86_64.sh install_packages

### install samtools 0.1.18
### using old version because the latest one could have somewhat trouble with
### other software like tophat.
wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
tar -jxvf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
make
cd ..
mv samtools-0.1.18.tar.bz2 install_packages


### install bowtie1 1.0.0
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip
unzip bowtie-1.0.0-linux-x86_64.zip
mv bowtie-1.0.0-linux-x86_64.zip install_packages

### install trim_galore
wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip
unzip trim_galore_v0.4.1.zip
mv trim_galore_v0.4.1.zip install_packages

### install bismark 0.7.6
wget http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.7.6.tar.gz
tar -zxvf bismark_v0.7.6.tar.gz
mv bismark_v0.7.6.tar.gz install_packages

### install bedtools 2.24.0
wget https://github.com/arq5x/bedtools2/releases/download/v2.24.0/bedtools-2.24.0.tar.gz
tar -zxvf bedtools-2.24.0.tar.gz
cd bedtools2
make
cd ..
mv bedtools-2.24.0.tar.gz install_packages

### install homer
# please Follow http://homer.salk.edu/homer/

### install HTSeq
pip install HTSeq

###install tabix and pytabix
wget http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
tar -jxvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6
make
cd ..
mv tabix-0.2.6.tar.bz2 install_packages
pip install pytabix
```

After that, download this script:
```
cd $PYTHONPATH  # path for put the python packages. path/to/anaconda/lib/python2.7/site-packages/ for default
git clone https://github.com/hubqoaing/NormSeq
```


Secondly, go to the ./setting file, and change the following values to your own path:
```python
self.Database       = "DIR/TO/DATABASE"          #line 43
self.sftw_py        = "DIR/TO/SOFTWARE_EXE_FILE" #line 65
self.sftw_pl        = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bgzip     = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_tabix     = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_ucsc_dir  = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_igvtools  = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_homer     = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_trim      = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bismark   = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bowtie_dir= "DIR/TO/SOFTWARE_EXE_FILE"
```

Go to the analysis dictionary and copy the bin file here.
``` bash
cd PATH/FOR/ANALYSIS   # go to
copy $PYTHONPATH/NormSeq/run_meth.py ./
```

Next, make the input files. You can download these files in UCSC or so on and then using own-scripts to merge the ERCC information, and generate files in this format.
``` bash
vim sample_input.xls
==> sample_input.xls <==
sample                                  stage   type        tissue  brief_name      merge_name
Sample_PD10_TFC_150713-mES-gWBS1-1      c       5mC_scBS    c       mESC_gWBS1_1    mESC_gWBS1
```
Notice that only NAME_FOR_RAW_FQ were required that this NAME should be the same as 00.0.raw_fq/NAME.
NAME_FOR_PROCESSING will be the name for the rest analysis's results.
NAME_FOR_READING    will be the name for files in statinfo.
stage and sample_group could be writen as anything. It was here only for make the downstream analysis easily.

Before running this pipeline, put the fastq reads in the ./00.0.raw_data dictionary.
```bash
mkdir 00.0.raw_data
for i in `tail -n +2 sample_input.xls | awk '{print $1}`
do
    mkdir 00.0.raw_data/$i && ln -s PATH/TO/RAW_DATA/$i/*gz 00.0.raw_data/$i
done
```

After that, running this pipeline:
```bash
python run_meth.py --ref YOUR_REF --cutSites GCA.GCC.GCT,ACG.TCG sample_input

```

Wait for the results.
Notice if you have to run it in a cluster, please do not running this scripts directly.
For example, if SGE system used, then:

Comments this command
```python
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
```
and using this command in modules in ./frame/*py
```
       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)
```
Method for submit jobs in other system were still developing.
