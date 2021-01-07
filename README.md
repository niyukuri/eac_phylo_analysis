# HIV viral diversity and transmission clusters in East Africa

Analysis and visualisation of results for HIV viral diversity and transmsission clusters across east African countries, with a simulation exercise to compare phylogenetic trees and transmission clusters built from protease and polymerase sequence data.



## CONTENTS

This repo contains the information necessary to reproduce the analysis and the simulation exercise:

* [Code and data files](#code-and-data-files)
   * Code files for the analysis
   * Data files 
   * Results files
   * Simulation files for the simulation exercise
* [System and software requirements](#system-and-software-requirements)
* [Contact information](#contact-information)

## CODE AND DATA FILES 


### Code files

All code is written in R. R is a statistical programming language and software package that is distributed under a GNU General Public License. R documentation and software is available for free download through the R Project for Statistical Computing website at http://www.r-project.org. The software is available as an executable file for a wide variety of operating systems and computer architectures, and a compilable binary is also available should you need it.

  ***code*** -- contains R scripts to analyse and visualise results

 
 ### Data files
  
  ***data*** -- contains a file with accession numbers at [LAN](<https://www.hiv.lanl.gov>), and Fasta files of sequence data.
  
  
### Results files

  ***results*** -- contains data objects of analysis outputs (phyloegentic trees), subfolder (***transmission_clusters***) for transmission clusters outputs.   


### Simulation files

  ***simpact_simulation*** -- contains R simulation script (***simpact_simulation.R***), root sequences files (hiv.seq.C.pol.j.fasta and hiv.seq.C.protease.j.fasta) to generate the sequences (C.Epidemic.pol.fas and C.Epidemic.protease.fas), and outputs when comparing phylogenetic trees built using Pol and protease sequences and transmission clusters.



## SYSTEM AND SOFTWARE REQUIREMENTS

### Operating system


  We ran the analysis on a personal computer (Linux Ubuntu Version 16.04).

### Required software

  **R version 3.4.4** <www.r-project.org> For statistical computing. To Install R, do the following:
  
  1. Open an internet browser and go to www.r-project.org.
  2.  Click the "download R" link in the middle of the page under "Getting Started."
  3. Select a CRAN location (a mirror site) and click the corresponding link.
  4. Click on the "Download R for ***your OS***" link at the top of the page.
  
  

  **FastTree version 2.1.10** <http://www.microbesonline.org/fasttree/#Install> Reconstructs a phylogenetic tree from a large alignment dataset. To install FastTree, do the following:
  
  1. Visit the website for downloading instructions: <http://www.microbesonline.org/fasttree/#Install>
  2. If you have a Linux operating system, you can directly download the executable files that are linked on that website. Those downloaded files can then be placed in your R working directory
  3. If you are using an OS X operating system, open the link "FastTree.c" in a new browser window
  4. Right-click on the program and click "Save as"
  5. Save anywhere on your computer
  6. Open the Terminal on your computer and change your working directory to the folder that contains "FastTree.c". After the prompt type:  `cd "file/path/here"`
  7. After the directory has been changed, after the prompt type: `gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm`
  8. Now check to see if a new executable file has been created in that folder
  9. Copy that file and paste it into your R working directory




  **Seq-Gen version 1.3.4** <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4> Simulates viral evolution across a transmission network and produces a sequence alignment. To install Seq-Gen, do the following:
  
  1. Visit the following Github repository to download the latest version of Seq-Gen: <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4>
  2. Click on the "Source Code" zip file to download
  3. Click on the zip file to unzip the folder
  4. Navigate to the source folder to confirm there is a file called "Makefile"
  5. Now you will need to compile the program using the Terminal on your computer
  6. Via the Terminal, change your working directory to the source folder by typing after the prompt: `cd "file/path/here/Seq-Gen-1.3.4 2/source"`
  7. Once your working directory has been set to the source folder, type after the prompt: `make`
  8. Now open the source folder and verify that a new file is present called "seq-gen"
  9. Copy that file and paste it into your R working directory
  
  If installing on HPC facility, you may follow the instructions from 1 up to 5. And you will load the tool via the the PBS file, for example `module add /apps/chpc/scripts/modules/bio/app/Seq-Gen/1.3.4`.
  

 **ClusterPicker version 1.2.3** <http://hiv.bio.ed.ac.uk/software.html> Cluster Picker identifies clusters in newick-formatted phyogenetic trees containing thousands of sequences. Cut-offs for within cluster genetic distance and bootstrap support are selected by the user.

  To use ClusterPicker, do the following:
  
  1. Install Java 1.6.0 or higher
  2. Visit the website for downloading instructions for ClusterPicker: <http://hiv.bio.ed.ac.uk/software.html>
  3. Donwload the ClusterPicker command line version


**FigTree version 1.4.4** <http://tree.bio.ed.ac.uk/software/figtree/> Graphical viewer of phylogenetic trees. 

  1. Visit the website for downloading instructions: <http://tree.bio.ed.ac.uk/software/figtree/>
  2. Compared binaries (for Mac, Windows and Linux) are available at <https://github.com/rambaut/figtree/releases>


**Note:** An executable file for FastTree and ClusterPicker must be in the working directory or in any other directory which will require appropriate resourcing them in the working directory. 


## CONTACT INFORMATION

David Niyukuri
Email: <niyukuri@sun.ac.za>




