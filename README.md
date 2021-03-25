# CKI
The central kinase inference (CKI) algorithm, a trans-omics-based computational method, was developed for computationally identifying central protein kinases (PKs) in a biological process.
<br>

# The description of each source code
### CKI.pl
The program infers central PKs from the integration of mRNA expression, substrate p-site intensity and network state for each PK. The usage of the code is shown as below: <br><br>
perl CKI.pl
<br>

### CKI-mRNA.pl
The program identifies PKs with differentially expressed mRNAs (DEMs) from the transcriptomic data. The usage of the code is shown as below: <br><br>
perl CKI-mRNA.pl
<br>

### CKI-Intensity.pl
Identifies PKs with differentially regulated p-site (DRP) profiles in substrates, based on an intensity-based approach and phosphoproteomic data. The usage of the code is shown as below: <br><br>
perl CKI-Intensity.pl
<br>

### CKI-Network.pl
Identifies PKs potentially associated with up- or down-regulated network modules, based on the network-based method and phosphoproteomic data. The usage of the code is shown as below: <br><br>
perl CKI-Network.pl
<br>

### Example_data
This folder contains example data and files for testing. 

# Computation Requirements
### OS Requirements
Above codes have been tested on the following systems: <br>
Windows: Windows 7, Windos 10<br>
Linux: CentOS linux 7.8.2003

### Software Requirements
Perl (v5.26.3 or later) program with modules (Statistics::Distributions).

### Hardware Requirements
All codes and softwares could run on a "normal" desktop computer, no non-standard hardware is needed.<br>
<br>

# Contact
Dr. Yu Xue: xueyu@hust.edu.cn<br>
Shaofeng Lin: linshaofeng@hust.edu.cn
