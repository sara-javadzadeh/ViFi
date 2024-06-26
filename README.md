# ViFi

ViFi is a tool for detecting viral integration and fusion mRNA sequences from Next Generation Sequencing data.  Unlike standard approaches that use reference-based read mapping for identification of viral reads, ViFi uses both reference-based read mapping and a phylogenetic-based approach to identify viral reads.  ViFi also incorporates mappability scores of the reads to filter out false positive integration detection.  The end result is a tool that can accurately and precisely detect integrated viruses, even if the viruses are highly mutated or novel strains.

ViFi is currently in alpha testing, is is constantly undergoing revisions.  High on the priority list is an easier installation process, as well as improve user interface.  Please report any problems/bugs to Nam Nguyen (ndn006@eng.ucsd.edu) so that ViFi can be improved and problems can be quickly corrected.  

The human reference dataset is updated in this forked repository as the original human reference was not accessible anymore. The updated human reference is based on hg38 assembly and is directly adapted from the Amplicon Architect project accessible at https://github.com/virajbdeshpande/AmpliconArchitect.

## ***UPDATE****

Due to major issues with incompatibilities between versions of Pysam and Samtools, Python versions, as well as issues with software compatibility between different platforms, we highly recommend that users  discontinue the use the Python version of ViFi, and instead, use the Dockerized version of ViFi.  The Dockerized version is platform independent and only requires Python (either version 2.7 or 3.0) and Docker to be installed, and no other software package is needed.  We outline below how to set up and install the Dockerized version, and how to run the Dockerized version.  

In addition, we include a [Tutorial] for all the different options within ViFi below.  We will include instructures
on how to run ViFi from the source code, but again, strongly discourage against this usage.

## Installation of ViFi for use in Docker

We provide instructions for preparing ViFi to be used for Docker below.  If Perl is installed,
the setup.sh script can be run that will automatically perform steps 3-7.  Note that ViFi requires a large
amount of diskspace to setup and run (10 Gb) due to the large size of the initial reference repositories.

1. Install Dependencies:
  1. Python (2.7 or 3.0; instructions for 2.7 is shown) 
  2. Docker (https://docs.docker.com/install/)
  
2. Download and run setup.sh (If Perl is installed and on Mac/Linux system).
Running this script will automatically download ViFi from GitHub, automatically download
the repositories from Google Drive, pull the latest ViFi docker image, set all the
environmental variables for ViFi, build the BWA index for
hg19+HPV via Docker, and run a test run of ViFi via Docker.  It can take up to an hour for the full
set of tests to complete and run.  Make sure you have at least 10 Gb of space free for the process to complete.

```
wget https://raw.githubusercontent.com/sara-javadzadeh/ViFi/master/setup_linux_mac.sh
sh setup_linux_mac.sh
```
  
Run steps 3-7 are only necessary if Perl is not installed on the machine or on Windows machine.  If Perl is on the machine, then setup_linux_mac.sh can be run to automatically set up ViFi (see to Step 2).

3. Clone the ViFi repository
```
git clone https://github.com/sara-javadzadeh/ViFi.git
```

4. Set the ViFi directory and include the python source to your Python path
```
echo export VIFI_DIR=/path/to/ViFi >> ~/.bashrc
echo export PYTHONPATH=/path/to/ViFi:/path/to/ViFi/src:$PYTHONPATH >> ~/.bashrc
```
5. Download the data repositories:
While we include some annotations, we are unable to host some large files in the git repository. These may be downloaded from https://drive.google.com/file/d/1XBZbwgcV1n2AWWAyt2RWfSKKxzssRFBo/view?usp=sharing. Thanks to Peter Ulz and Shiting Li for noticing incorrect link earlier.
```
tar -zxvf data_repo.tar.gz
echo "GRCh38" > ./data_repo/reference.txt
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
source ~/.bashrc
```
6. Download the HMM models:
We have pre-build HMM models for HPV, HBV, HCV and EBV.  They are included in the GitHub repository in a compressed format or alternatively can be downloaded from https://drive.google.com/file/d/1VB0qNHRM--CLgZAMmn8sFcNqjzE7QbCn/view?usp=sharing.
```
tar -zxvf viral_data.tar.gz
echo export REFERENCE_REPO=$PWD/viral_data >> ~/.bashrc
source ~/.bashrc
```
7.  Build a BWA index on the reference sequences from human+viral sequences:
We show an example of building an index of human+viral sequences using Hg19 and **HPV** and **HBV** below.  However
any reference organism+viral family could be used.
```
cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hpv/hpv.unaligned.fas > $REFERENCE_REPO/hpv/hg19_hpv.fas
bwa index $REFERENCE_REPO/hpv/hg19_hpv.fas
cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hbv/hbv.unaligned.fas > $REFERENCE_REPO/hbv/hg19_hbv.fas
bwa index $REFERENCE_REPO/hbv/hg19_hbv.fas
```  

## Running ViFi using Docker (RECOMMENDED)

We have also created a dockerized version of ViFi to enable easier time running (see previous section for installation and setup).  To get the latest version of the Dockerized ViFo, run:
```
docker pull namphuon/vifi
```

To run the dockerized version of ViFi, first create the data repositories as above, including setting the environmental variables.  Next, run the following command:

`python $VIFI_DIR/scripts/run_vifi.py -f <READ1> -r <READ2> --docker`

where <READ1> and <READ2> are the FASTQ files (gzipped or unzipped).  Note that the $VIFI_DIR, $AA_DATA_REPO and $REFERENCE_REPO variables must be set in order for the script to find the necessary files.  

Example (assuming that $VIFI_DIR is set):

```
python $VIFI_DIR/scripts/run_vifi.py -f $VIFI_DIR/test/data/test_R1.fq.gz -r $VIFI_DIR/test/data/test_R2.fq.gz  --docker
```


## ViFi Output

The output of ViFi is the list of read clusters discovered, and for each read cluster, the relaxed, stringent, and exact (if split reads are present) ranges are reported, aswell as the read names of the reads in the cluster.

The main output files of interest are
- \<prefix\>.clusters.txt
- \<prefix\>.clusters.txt.range

\<prefix\>.clusters.txt is a tab delimited file that reports the human integration range, the number of reads supporting the integration, and the number of reads mapped to the forward/reverse strand of the human region, as well as the number of viral reads mapping to the virus sequence.  It also includes the names of each discordant read supporting the integration.

Below is the sample:
```
#chr    minpos  maxpos  #reads  #forward        #reverse
##================================================================
chr19   36212224        36212932        7       4       3
##ERR093797.9977893     chr19   36212224        True    False
##ERR093797.7073606     chr19   36212403        True    True
```

The first line is the header information.  Afterward, each integration cluster is separated by a line
containing **=**.  The first line of an integration cluster describes the following:

1. Reference chromosome (chr19)
2. Minimum reference position of all mapped reads belonging to that cluster (36212224)
3. Maximum reference positions of all mapped reads belonging to that cluster (36212932)
4. Number of read pairs belonging to this cluster (7)
5. Number of reads mapped to the forward reference strand (4)
6. Number of reads mapped to the forward reference strand (3)

After this line, each read pair that mapped to this cluster is displayed.  The information is
1. Read name (ERR093797.9977893)
2. Reference chromosome (chr19)
3. Starting read map location (36212224)
4. Read is on the reverse strand (True)
5. Read is read1 (False)

\<prefix\>.clusters.txt.range is a much more condensed summary of the results, showing just the integration range on the
human reference (based upon discordant reads) and attempts to identify the exact integration point if split reads are available.

Below is a sample:
```
Chr,Min,Max,Split1,Split2
chr19,36212564,36212564,-1,-1
```

The first line is header information.  Afterward, each line is information about the cluster.  For example,

1. Reference chromosome (chr19)
2. Minimum reference position of all mapped reads belonging to that cluster (36212224)
3. Maximum reference positions of all mapped reads belonging to that cluster (36212932)
4. If split read exists, minimum split read mapped range, -1 if no split read exists (-1)
5. If split read exists, maximum split read mapped range, -1 if no split read exists (-1)

Finally, ViFi outputs several working files that can be deleted after a run.  These are:
1. hmms.txt - The list of HMM files used during the run
2. \<prefix\>.bam - The aligned (name-sorted order) BAM file containing the input reads
3. \<prefix\>.unknown.bam - A BAM file containing all paired reads in which one or both paired end reads that did not align to any known reference.  ViFi will then search these reads against the HMMs to identify any viral reads.
4. \<prefix\>.viral.bam - A BAM file containing all paired reads that only aligned to viral references
5. \<prefix\>.viral.cs.bam - A coordinate sorted BAM file containing all paired reads that only aligned to viral references
6. \<prefix\>.trans.bam - A BAM file containing all paired reads in which one read aligned to the human and the other aligned to the viral reference.
7. \<prefix\>.fixed.trans.bam - A BAM file created by merging 6. and any human/viral paired end reads discovered by running the viral HMMs on 3.
8. \<prefix\>.fixed.trans.cs.bam - A coordinate sorted BAM file of 7.

## References
1. Nguyen ND, Deshpande V, Luebeck J, Mischel PS, Bafna V (2018) ViFi: accurate detection of viral integration and mRNA fusion reveals indiscriminate and unregulated transcription in proximal genomic regions in cervical cancer. Nucleic Acids Res (April):1–17.

# [Advanced Notes](#advanced_notes)

## Building evolutionary models

ViFi can be run with and without evolutionary models (i.e., the HMMs).  We outline the steps in building the HMMs below.  However, we also include a Docker pipeline that will automatically build the HMMs for the users to use.  The pipeline only requires docker to be installed for use.

## Using Docker pipeline to build HMMs for use in ViFi
The following command will create HMMs from a set of unaligned sequences.  The sequences are assumed to share a common viral ancestor (i.e., don't mix viral families together when running the pipeline).  

```
bash $VIFI_DIR/scripts/build_references.sh <INPUT_SEQ> <OUTPUT_DIR> <PREFIX>
```

The output in the OUTPUT_DIR folder will be a set of HMMs (suffix with *.hmmbuild) and a file containing the list of HMMs.
 
## Using customized reference 

If you want to use a customized reference or a reference for a different organism, you can inform ViFi 
of the reference sequences by supplying a chromosome file to ViFi using the **--chromosome_list**.  The file 
format is a single line that has the sequence names delimited by spaces.  For example:

```
mouse_chr1 mouse_chr2
```

would inform ViFi that any other sequences found in the BAM file that does not match mouse_chr1 and mouse_chr2 are
considered viral sequences.

## Installation from source code (Depreciated):
We provide instructions for installing ViFi on Linux below.  

1. ViFi download (if you have not already cloned this source code):
```
git clone https://github.com/sara-javadzadeh/ViFi.git
```
2. Install Dependencies:
   1. Python 2.7
   ```
   sudo dnf install python2
   ```
   2. Pysam verion 0.9.0 or higher (https://github.com/pysam-developers/pysam):
   ```
   sudo pip install pysam
   ```
   3. Samtools 1.3.1 or higher (www.htslib.org/)
   ```
   sudo apt-get install samtools
   ```
   4. BWA 0.7.15 or higher (bio-bwa.sourceforge.net/)
   ```
   sudo apt-get install bwa
   ```
   5. Install HMMER v3.1b2 and have it on the path (http://hmmer.org/)
   ```
   sudo apt-get install hmmer
   ```
3. Set the ViFi directory and include the python source to your Python path
```
echo export VIFI_DIR=/path/to/ViFi >> ~/.bashrc
echo export PYTHONPATH=/path/to/ViFi:/path/to/ViFi/src:$PYTHONPATH >> ~/.bashrc
```
4. Download the data repositories:
While we include some annotations, we are unable to host some large files in the git repository.  These may be downloaded from https://drive.google.com/file/d/1il10KUxJ5Q5JvR5pHJB4GUMBlBPgTjrj/view?usp=sharing. Thanks to Peter Ulz and Shiting Li for noticing incorrect link earlier.
```
tar -zxvf data_repo.tar.gz
echo "GRCh38" > ./data_repo/reference.txt
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
source ~/.bashrc
```
5. Download the HMM models:
We have pre-build HMM models for HPV, HBV, HCV and EBV. They are included in the GitHub repository in a compressed format or alternatively can be downloaded from https://drive.google.com/file/d/1VB0qNHRM--CLgZAMmn8sFcNqjzE7QbCn/view?usp=sharing.
```
tar -zxvf viral_data.tar.gz
echo export REFERENCE_REPO=$PWD/viral_data >> ~/.bashrc
source ~/.bashrc
```
6.  Build a BWA index on the reference sequences from human+viral sequences:
We show an example of building an index of human+viral sequences using Hg19 and **HPV** and **HBV** below.  However
any reference organism+viral family could be used.
```
cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hpv/hpv.unaligned.fas > $REFERENCE_REPO/hpv/hg19_hpv.fas
bwa index $REFERENCE_REPO/hpv/hg19_hpv.fas

cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hbv/hbv.unaligned.fas > $REFERENCE_REPO/hbv/hg19_hbv.fas
bwa index $REFERENCE_REPO/hbv/hg19_hbv.fas
```
## Running ViFi  (Depreciated)

We show the most basic example of running ViFi below.  This version assumes that the user has
followed all the previous steps.  More advanced options, such as using a customized reference organism/viral
family is provided in the [Advanced Notes](#advanced_notes) section.  
```
python run_vifi.py -f <input_R1.fq.gz> -r <input_R2.fq.gz> -o <output_dir>
```

Note that this version defaults to searching for **HPV**.  To search for HBV, run the following command.
```
python run_vifi.py -f <input_R1.fq.gz> -r <input_R2.fq.gz> -o <output_dir> -v hbv
```
