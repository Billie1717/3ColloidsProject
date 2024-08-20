from __future__ import division, print_function
import numpy as np
#import matplotlib.pyplot as plt
import math
import time
import random
import sys, getopt
import os

#'''
#
#'''

def positions(arc1,number, radius):
    x = []
    y = []
    z = []
    #taurads = tau*(np.pi/180)
    for i in range(number):
        #if i == 2:
        #    l = arc2
        #    phi = (l/radius)
        #    x.append(radius*np.sin(phi)*np.cos(taurads))
        #    y.append(radius*np.sin(phi)*np.sin(taurads))
        #    z.append(radius*np.cos(phi))
        
        l = i * arc1 - 0.5*arc1
        phi = (l/radius)
        x.append(radius*np.sin(phi))
        y.append(0)
        z.append(radius*np.cos(phi))
    
    
    return x,y,z
    
def main():
    ################ Spherical Membrane ################
    Nmem = 5882
    arc1 = float(sys.argv[1])
    #arc2 = float(sys.argv[2])
    #tau = int(sys.argv[3])
    sigma_c = float(sys.argv[2])
    rc_c = 1.2 #(sigma_c+2)/(sigma_c+1) #ensures rcut ends up being the radius of the colloid plus 1 membrane bead
    kappa = 15
    press  = 1.0
    D0 = float(sys.argv[3])
    seed = int(sys.argv[4])
    
    print('intput params', arc1,D0)
    
    #filepattern = sys.argv[1] #'colloids_arc1_%d_arc2_%d_tau_%d_D0_%d'%(arc1,arc2,tau,D0)
    #print('filepattern',filepattern)
    workingdir = '/nfs/scistore15/saricgrp/bmeadowc/Scratch/Teather/MC_mem_clean/'
    fname_in = workingdir +'VesicalFiles/icos_%d_0_8.dat'%(Nmem) 			# read in positions of membrane particles
    f_in = open(fname_in,'rU')
    #tmp = f_in.readline()
    #tmp = f_in.readline()					#? two times necessary? YES - skip 2 lines
    #lbox = int(tmp.split()[1])
    #lbox = int(tmp.split()[1])+12
    x_mem = []
    y_mem = []
    z_mem = []
    r_mem = []
    while True:
        tmp = f_in.readline()
        if tmp:
            x_mem.append(float(tmp.split()[1]))
            y_mem.append(float(tmp.split()[2]))
            z_mem.append(float(tmp.split()[3]))
            r_mem.append(np.sqrt(float(tmp.split()[1])**2+float(tmp.split()[2])**2+float(tmp.split()[3])**2))
        else:
            break

    cell_radius = np.mean(r_mem)				# radius of the cell (-1 such that the generated filament definitely is inside the cell)
    #Nmem = len(x_mem)
    zMAX = max(z_mem)
    print(zMAX)
    xCM = np.sum(x_mem)/float(Nmem)
    yCM = np.sum(y_mem)/float(Nmem)
    zCM = np.sum(z_mem)/float(Nmem)
    print('Center of mass: [%.3f,%.3f,%.3f]'%(xCM,yCM,zCM))
    x_mem = x_mem - xCM
    y_mem = y_mem - yCM
    z_mem = z_mem - zCM
    xCM = np.sum(x_mem)/float(Nmem)
    yCM = np.sum(y_mem)/float(Nmem)
    zCM = np.sum(z_mem)/float(Nmem)
    print('Cell radius: %.3f'%(cell_radius))
    print('Center of mass: [%.3f,%.3f,%.3f]'%(xCM,yCM,zCM))

    ################ Cargo particles ################
    Ndumb = 2
    radius1 = cell_radius + sigma_c*0.8
    radius2 = 22.695 +sigma_c*0.2 #cell_radius + sigma_c*0.5 23.195
    x_bead, y_bead, z_bead = positions(arc1, Ndumb, radius1)
    x_beadSp, y_beadSp, z_beadSp = positions(arc1, Ndumb, radius2)
    Ntot = Nmem + Ndumb
    
    print('Total # of particles in the system: %d\n\tMembrane: %d\n\tDumbbell: %d'%(Ntot,Nmem,Ndumb))

    ################  Write particles.in ################
    print('Writing icos_%d_%d_%d.dat\n'%(Nmem,Ndumb,sigma_c))
    #os.mkdir(workingdir + 'runs_Rshift/teth_{:d}_mass_{:.1f}_kF_{:.2f}_kS_{:.2f}_arc1_{:.2f}_arc2_{:.2f}_tau_{:.2f}'.format(seed,mass,k_fixed,k_samp,arc1,arc2,tau))
    #outfile=workingdir + 'runs_Rshift/teth_{:d}_mass_{:.1f}_kF_{:.2f}_kS_{:.2f}_arc1_{:.2f}_arc2_{:.2f}_tau_{:.2f}/particles.in'.format(seed,mass,k_fixed,k_samp,arc1,arc2,tau) 
    outfile = 'icos_%d_%d_%d.dat'%(Nmem,Ndumb,sigma_c)
    f=open(outfile,'w')
    for ii in range(Nmem):
        r = float(np.sqrt(x_mem[ii]**2+y_mem[ii]**2+z_mem[ii]**2))
        f.write('1' + '\t' + str(x_mem[ii]) + '\t' + str(y_mem[ii]) + '\t' + str(z_mem[ii]) +'\n')
    for ii in range(Ndumb):
        f.write(str(2) + '\t' + str(x_bead[ii]) + '\t' + str(y_bead[ii]) + '\t' + str(z_bead[ii]) + '\n')
    f.write('\n')
    f.close()
    ################  Write in.ves ################
    
    print('Writing in.ves\n')
    infile = 'in.ves'
    fin=open(infile,'w')
    fin.write(str(kappa)+' '+str(press)+' ' +str(Nmem)+' '+str(Ndumb)+' '+str(sigma_c)+' '+str(D0)+' %.1f'%(rc_c)+' '+str(seed)+'\n')
    fin.write('kappa\npressure\nNmem\nNcolloid\nsigma_colloid\nD0\nrc_colloid\nseed\n ' )
    fin.close
    
    print('Writing submission file\n')
    runscriptname = 'runscript.sh'
    fsub = open(runscriptname, "w")
    fsub.write(
        '''#!/bin/bash
#SBATCH --job-name=Job_%s
#SBATCH --output=OUTFILE_%s.dat
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --mail-user=bmeadowc@ist.ac.at
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1
\n'''%(outfile, outfile))    #cd $TMPDIR
    fsub.write("/nfs/scistore15/saricgrp/bmeadowc/Scratch/Teather/MC_mem_clean/scripts/Vescicle_EPB_averaged2Beads in.ves\n")
    #fsub.write("/nfs/scistore15/saricgrp/bmeadowc/Scratch/Teather/MC_mem_clean/scripts/Vesc_EPerBead in.ves\n") 
    #fsub.write("/nfs/scistore15/saricgrp/bmeadowc/Scratch/Teather/MC_mem_clean/scripts/Vescicle_EnPerBead in.ves\n") #change to lammps dr

    fsub.close()
    
    
if __name__ == "__main__":
    main()
