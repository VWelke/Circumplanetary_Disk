This model follows from CPD_dust_dias_ne8_nphi116


1. try nphi = 120 to see if it breaks the barrier

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
#############################################################################
// run in interactive bash
srun --mem-per-cpu=4G -c 16 --pty bash -i
bash trial2.sh


1. try reduce nr from 150 to 100  -> fail
2. reduce nr to 70  -> fail
3. nr go back to 150, reduce nphi from 116 to 80 -> tiny bit faster, still very very slow, fail
4. nphi go back to 116, now try r array use back log space -> fail 
5. Keep r at log space -> open a new file to cofirm the surface density plot -> Yes it seems correct
6. try to change other parameters: r_min, r_max,phi_min, phi_max, sigmad02, 
7. Try increase the mask region to 5 au 
    r_min = 32.2, r_max = 42.2, phi_max = np.arctan(5/37.2)-> failure
8. Go back to r_min = 35.2, r_max = 49.2, phi_max = np.arctan(2/37.2), try to copy the CPD_dust_dias_ne7 which works and change the sigmad02 to 2127 
    cp -r CPD_dust_dias_ne7 CPD_sig2127   -> succeed still very fast
9. add the plot_sigmad snippet and change r to linspace -> works, still very fast, the plot doesnt show up though
10. Change the slope to -1.2  -> works still very fast
11. Change the mask region to 2 au CPD radius -> fail, become very slow
12. Change back the mask region to 1 au CPD radius first, and rerun 
13. let mask region size be 0 -> commented out all the masked sigma, hh, rhod
Plot r_shifted vs r_masked and pp?
14. Change back mask region and run the simulation at 10^8 photons
#############################################################################

Code here to obtain radmc3d image
And what to check after that


cd /c/Users/LHEM/Desktop/Van_Code_Projects/Circumplanetary_Disk/CPD_dust_dias_newsigma

python mockmodel.py

radmc3d mctherm setthreads 4

radmc3d image incl 0 phi 0 loadlambda

radmc3d sed incl 0 phi 0 zoomau 35.2 39.2 -2 +2
#radmc3d spectrum incl 0 phi 0 zoomau 35.2 39.2 -2 +2 lambdarange 350. 7500. nlam 100

#################################################################################

check the temperature plots, sed, image, 
