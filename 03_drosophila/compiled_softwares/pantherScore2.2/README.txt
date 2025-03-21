Copyright (C) 2007 Paul Thomas
This file may be copied and redistributed freely, without advance permission,
provided that this Copyright statement is reproduced with each copy. 

LIMITATION OF WARRANTY
NOTHING IN THIS AGREEMENT WILL BE CONSTRUED AS A REPRESENTATION MADE OR
WARRANTY GIVEN BY PAUL THOMAS OR ANY THIRD PARTY THAT THE USE OF
DATA PROVIDED HEREUNDER WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK
OR OTHER RIGHTS OF ANY THIRD PARTY. DATA IS PROVIDED "AS IS" WITHOUT
WARRANTY OF ANY KIND WHATSOEVER, EXPRESS OR IMPLIED, INCLUDING IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. PAUL
THOMAS MAKES NO WARRANTY THAT ITS DATA DOES NOT CONTAIN ERRORS.

############################################################################

PANTHER HMM scoring tools - Version 2.2
http://www.pantherdb.org/downloads/
6/18/2019

##########
Introduction: 

This tool is used to score protein sequences against the PANTHER HMM library, to help infer protein function.

If you have any questions, please contact us at: pantherfeedback@yahoo.com
##########
Requirements:

1. PANTHER HMM library (ftp://ftp.pantherdb.org/panther_library/14.1/PANTHER14.1_hmmscoring.tgz)
2. UNIX
3. Perl
4. HMMER3 - Download from   
http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz

the location to Perl, and HMMER binaries must be defined in your $PATH variable or specified by the users in the arguments.  If you have any questions on how to set up $PATH, please contact your UNIX system administrator.

##########
Usage:
% tcsh
% cd pantherScore2.2
% source panther.cshrc
% ./pantherScore2.2.pl -l <panther_hmm_library> -D B -V -i <fasta file> -o <output file> -n
or % ./pantherScore2.1.pl -l <panther_hmm_library> -D B -V -i <fasta file> -o <output file> -n -s

-panther.cshrc is a sample cshrc file.  You just need to make sure the lib directory is in your PERL5LIB variable path 
-alternatively, for the -D option, use A, to get back additional hits, rather than just the best hit
-an example fasta file was included called test.fasta

##########
Sample output:
test    PTHR19264:SF179 G-PROTEIN COUPLED RECEPTOR HOMOLOG      3e-246  828.9
1-370,

-tab delimited file in the following format:
  col 1 - sequence ID
  col 2 - PANTHER accession if (containts :SF, is a subfamily HMM)
  col 3 - PANTHER family or subfamily name
  col 4 - HMM evalue score, as reported by HMMER
  col 5 - HMM bit score, as reported by HMMER. It is used as the primary score to select the best hit. 
  col 6 - alignment range of protein for this particular HMM

##########

Notes: 
-a subset of the PANTHER HMMs have been released to the InterPro.  If you have obtained this README file through the InterProScan download and would like to download all the PANTHER HMMs to score your proteins against, please visit: http://pantherdb.org/downloads/

-pantherScore2.2 use hmmscan program within HMMER3 as default for scoring input sequences against PANTHER14.1 hmm models built with HMMER3 hmmbuild program. You can overwrite this default by use flag -s to use hmmsearch program. hmmscan run faster for small number of sequences (<100), for large number of sequences use hmmsearch.

########## 
Interpretation of scores:

-The HMM bit score is used as the primary score for selecting the best hit. The larger the number, the better the score.
-E-values can be used as a reference for relatedness also as below:
	-closely related: if the score is better than E-23 (very likely to be a correct functional assignment)  
	-related : if the score is better than E-11, but worse than E-23 (molecular function likely to be the correct but biological process/pathway less certain) 
	-distantly related : if the score is better than E-3, but worse than E-11 (protein is evolutionarily related but function may have diverged) 


##########

Troubleshooting:

Ultimately, if you have any problems, please contact us at: feedback@pantherdb.org

but, before you do that it would be helpful if you can try each of the the following commands:
% hmmsearch PANTHER14.1/books/PTHR18976/hmmer.hmm test.fasta
% hmmsearch -Z 10000 PANTHER14.1/books/PTHR18976/hmmer.hmm test.fasta
% hmmsearch --cpu 1 PANTHER14.1/books/PTHR18976/hmmer.hmm test.fasta

All of these commands should run properly and generate results.  If they do not, this means that you have a problem with the way you installed hmmsearch, or more likely, you have a problem with how you compiled HMMER.

In particular, if you have problems with the --cpu option (if you have problems with this option you might see a POSIX or threads error), you should download and recompile HMMER so that it properly works with threads.  HMMER can be downloaded from http://hmmer.wustl.edu/.

If you send an email to PANTHER feedback, please tell us the command you are using, and send us the fasta file you are using.

##########
Version History:

version 2.2 - fixed bugs for sequence range output, and changed to use bit score as the primary score for selecting the best hit.
version 2.1 - fixed bugs for sequence range output
version 2.0 - hmm model build using HMMER3, update this script to score using hmmscan or hmmsearch of HMMER3
version 1.03 - update hmm.pm and blast.pm modules in the lib/ directory to support the changes in the HMMER2 and BLAST algorithms.
version 1.02 - updated document, so that users can better trouble shoot errors with running the program (specifically, the threads issue)  8/23/07
version 1.01 - minor change.  if gaps between alignment ranges are <= 15 positions, then the alignment ranges are joined.


If you have any questions, please contact us at: pantherfeedback@yahoo.com

