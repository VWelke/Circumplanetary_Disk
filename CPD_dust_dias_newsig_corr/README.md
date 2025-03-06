This model follows from CPD_dust_dias_newsigma, because I just realised I commented outt the rhobg[masked]




#############################################################################
UCL cluster

ssh xzcapvwe@dias.hpc.phys.ucl.ac.uk 
password:  Oowei5ohbiSh

scp -r "C:\Users\LHEM\Desktop\Van_Code_Projects\Circumplanetary_Disk\CPD_dust_dias_newsigma" xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_newsigma 


cd MSci_project/models/CPD_newsigma

sbatch trial2.sh

scp  xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_sig2127/mockmodel.py "C:\Users\LHEM\Desktop\Van_Code_Projects\Circumplanetary_Disk\CPD_dust_dias_newsigma" 

scp xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_sig2127/image.out "C:/Users/LHEM/Desktop/Van_Code_Projects/Circumplanetary_Disk/CPD_dust_dias_newsigma"


scp xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_sig2127/radmc3d.out "C:/Users/LHEM/Desktop/Van_Code_Projects/Circumplanetary_Disk/CPD_dust_dias_newsigma"


scp xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_sig2127/dust_temperature.dat "C:/Users/LHEM/Desktop/Van_Code_Projects/Circumplanetary_Disk/CPD_dust_dias_newsigma"





-rw-rw-r-- 1 xzcapvwe xzcapvwe  2138732 Feb 27 01:00 image.out
-rw-rw-r-- 1 xzcapvwe xzcapvwe       35 Feb 26 23:49 radmc3d.out
-rw-rw-r-- 1 xzcapvwe xzcapvwe 86443239 Feb 26 23:49 dust_temperature.dat



Fuck I just realised I still have the masked rhobg commented out , damn  gotta run this again
- So now, given that the CPD is located far enough, the big dust follows the drop off from PPD is should be very small locally, can I say that because of dust filtration, only small grains are here? is 100 microns dust big enough to be filtered out?

Apparently if calculate the Stroke's number , it will not be filtered, so I need to rerun
#1. Copy it to new directory
cp -r CPD_sig2127 CPD_rhobg
ran, hai
2. remove the big files and move here the directory
scp -r xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_rhobg "C:\Users\LHEM\Desktop\Van_Code_Projects\Circumplanetary_Disk\CPD_dust_dias_newsig_corr"  


scp xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_rhobg/radmc3d.out "C:/Users/LHEM/Desktop/Van_Code_Projects/Circumplanetary_Disk/CPD_dust_dias_newsig_corr"


scp xzcapvwe@dias.hpc.phys.ucl.ac.uk:~/MSci_project/models/CPD_rhobg/dust_temperature.dat "C:/Users/LHEM/Desktop/Van_Code_Projects/Circumplanetary_Disk/CPD_dust_dias_newsig_corr"