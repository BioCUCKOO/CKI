# CKI
The central kinase inference (CKI) algorithm, a trans-omics-based computational method, was developed for computationally identifying central protein kinases (PKs) in a biological process.
<br>

# The description of each source code
### CKI.pl
The program infers central PKs from the integration of mRNA expression, substrate p-site intensity and kinase-substrate network for each PK. The usage of the code is shown as below: <br><br>
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
All codes could run on a "normal" desktop computer, no non-standard hardware is needed.<br>

# Instruction
For users who want to run CKI in own computer, files contain the information of mRNAs and p-sites should be first generated as the form demonstrated in Example_data/.<br>
"All.igps.txt" contains the upstream regulatory PKs of p-sites predicted by iGPS (http://igps.biocuckoo.org/). <br>
"mRNA-DOX-NT.txt" contains the information of PKs in transcriptional level. <br>
"Site-Raw.txt" contains all p-site with intensities of samples, and split into single file of individual sample.<br>
At last, once you have all above files ready, you can open the command prompt in your computer with Perl v5.26.3 installed, and open the directory in which you hava saved perl codes. After run the command "perl CKI.pl" to cpmpile the perl program, the result file "CKI-kinase.txt" was generated which contains the outputs of central protein kinases.

# Additional information
Expected run time is depended on the number of p-sites quantified, it will take about 3 minutes for 5,000 sites.

# Contact
Dr. Yu Xue: xueyu@hust.edu.cn<br>
Dr. Pengyu Huang: huangpengyu@yeah.net<br>
Dr. Yilai Shu: yilai_shu@fudan.edu.cn<br>
Shaofeng Lin: linshaofeng@hust.edu.cn
