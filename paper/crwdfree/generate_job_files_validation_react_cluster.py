#!/usr/bin/python

"""
Generate job files for the geek_model
"""

from pandas import DataFrame
import numpy as np
from scipy.integrate import quad
import time
import os
import sys

args = sys.argv


def mass2rad(mass):
   radius = 0.0515*(mass*1000)**(0.393)
   return radius

def rad2mass(radius):
   M = (radius/0.0515)**(1./0.393)/1000.0
   return M

def radius_distribution(radius,mu,sigma):
   # radius in nm
   M = rad2mass(radius)
   p = 1/(np.sqrt(2.0*np.pi)*sigma*M)*np.exp(-1.*(np.log(M/mu)**2)/(2.0*sigma**2))
   return p



AVOGADRO_NUMBER = 6.2e23


volume = 10e-18 # (0.1 mum)^3 in L

"""Parameter to be varied """
mu_sigma = [(21.1, 0),]
realizations = 10
volume_fractions = [0.0,0.2,0.1,0.3,0.4]

# saturations =  np.arange(0.1,1.0,0.1)


# Reference concentrations
A_ref = 50e-6
B_ref = 50e-6
C_ref = 50e-6

A_concentrations = np.array([0.25, 0.5, 1., 2.0, 4.0 ])*1.0
B_concentrations = np.array([0.25, 0.5, 1., 2.0, 4.0 ])*1.0
C_concentrations = np.array([0.25, 0.5, 1., 2.0, 4.0 ])*1.0



""" Parameters for the cluster """
n_cpu = 28
nodes = 1
debug = False

# Time stamp to identifiy all input files
timestamp = time.ctime().replace(" ", "_").replace(":","_")
# Write all to single data
cwd = os.getcwd()

input_file = cwd+"/validation_difflim_cluster.py"
folder = cwd+"/input_case_study_difflim_"+timestamp+"/"
log_folder = args[1]+"/out_"+timestamp
os.mkdir(log_folder)
os.mkdir(folder)

output = args[1]+"/run_"+timestamp
os.mkdir(output)

job_id = "$SLURM_JOB_ID"



def round_to_e5(x):
    return round(x, 5-np.int(np.floor(np.log10(np.abs(x)))))

# Generate a data frame with all inputs

# Coummns ordered according to simulation file
# docer/work/geek_model/pgm_cluster.py
columns =[ "simulation_type" ,
           "volume_fraction",
           "realization",
           "A_concentration",
           "B_concentration",
           "C_concentration",
           "data"]



class GetOutOfLoop( Exception ):
    pass

counter = 0
chunksize = n_cpu*nodes
complete_data_frame = DataFrame(columns = columns)
data_frames = [DataFrame(columns = columns),]
cpu_times = []
n_particles = []
this_job_file_index = 0

try:
   for this_mu, thi_sigma in mu_sigma:
      for this_volume_fraction in volume_fractions:
         for this_A in A_concentrations:
             for this_B in B_concentrations:
                 for this_C in C_concentrations:
                         for this_realization in range(realizations):

                            # Create job File
                            this_job = {}
                            this_job["simulation_type"] = 'diff'
                            this_job["realization"]     = this_realization
                            this_job["volume_fraction"] = this_volume_fraction


                            if this_A > 0:
                               this_job["A_concentration"]  = round_to_e5(this_A)
                            else:
                               this_job["A_concentration"]  = this_A

                            if this_B > 0:
                               this_job["B_concentration"]  = round_to_e5(this_B)
                            else:
                               this_job["B_concentration"]  = this_B

                            if this_C > 0:
                               this_job["C_concentration"]  = round_to_e5(this_C)
                            else:
                               this_job["C_concentration"]  = this_C


                            this_job["data"] = output
                            #this_job["job_id"] = job_id
                            complete_data_frame = complete_data_frame.append(this_job, ignore_index=True)

                            data_frames[this_job_file_index] = data_frames[this_job_file_index].append(this_job, ignore_index=True)

                            # Calcualte the number of particles
                            n = volume*(this_A+this_B+this_C)*AVOGADRO_NUMBER
                            if thi_sigma == 0:
                               var = this_volume_fraction*volume/(4./3.*np.pi*(mass2rad(this_mu)*1e-8)**3)
                               n += var
                               # Beaware of scling nm to dm

                            else:
                               f2 = lambda r: radius_distribution(r,this_mu,thi_sigma)
                               const = quad(f2, 0.0 , 10.0 )[0]
                               f = lambda r: radius_distribution(r,this_mu,thi_sigma)/const*(r*1e-8)
                               mean_r = quad(f, 0.0 , 10.0 )[0]

                               var = this_volume_fraction*volume/(mean_r**3*4/3*np.pi)
                               n += var

                            print('Number of particles {} \t {} \t sigma {}'.format(n,var,thi_sigma))
                            n_particles.append(n)

                            counter += 1

                            if debug and (this_job_file_index > 3):
                                  raise(GetOutOfLoop)

                            if counter >= chunksize:

                               counter = 0

                               cpu_time = 60.0 * 70.0 # Seconds
                               cpu_times.append(cpu_time)

                               data_frames.append(DataFrame(columns = columns))
                               this_job_file_index += 1
                               n_particles = []

except GetOutOfLoop:
   pass


#complete_data_frame.to

def chunk_data_frame(this_data_frame,chunksize):
   counter = 0
   this_job_file_index = 0
   data_frames = [DataFrame(columns = columns),]
   for index, row in this_data_frame.iterrows():
      counter += 1
      data_frames[this_job_file_index] = data_frames[this_job_file_index].append(row, ignore_index=True)

      if counter >= chunksize:
        counter = 0
        data_frames.append(DataFrame(columns = columns))
        this_job_file_index += 1
   return data_frames


from  datetime import timedelta

# Cut the dataframe into input files with n_node*n_cpu jobs
def write_sbatch_file(data_frame,input_file,folder,this_id,cpu_time):
   conf_file_name = "particle_model_cluster_"+str(this_id)+".sbatch"
   file_obj = open(folder+conf_file_name,'w')

   file_obj.write("#!/bin/bash \n")
   file_obj.write("#SBATCH --workdir "+log_folder+"\n")
   file_obj.write("#SBATCH --nodes 1 \n")
   file_obj.write("#SBATCH --ntasks 28 \n")
   file_obj.write("#SBATCH --mem 64G \n")

   file_obj.write("#SBATCH --time "+str(timedelta(seconds=round(cpu_time)))+'\n')
   if debug:
      file_obj.write("#SBATCH --partition debug \n")
   else:
      file_obj.write("#SBATCH --partition parallel \n")
   file_obj.write("#SBATCH --account lcsb \n")
   file_obj.write("\n")
   file_obj.write("module load gcc/7.3.0 \n")
   file_obj.write("module load python/3.6.5 \n")
   file_obj.write("module load openblas \n")
   file_obj.write("\n")
   file_obj.write("source /home/weilandt/virtual_env/openbread/bin/activate \n")
   file_obj.write("\n")
   
   for index, row in data_frame.iterrows():
      this_program = "python3"
      this_input_args = ""
      for arg in row:
          if isinstance(arg, float):
              this_input_args += format(arg,'.5e')+" "
          else:
              this_input_args += str(arg)+" "

      this_line = this_program+" "+input_file+" "
      this_line = this_line + " " +this_input_args +" "+str(index)+" $SLURM_JOBID & \n" + "sleep 0.1 \n"
      file_obj.write(this_line)

   file_obj.write("wait\n")
   file_obj.write("deactivate\n")
   file_obj.close()



#from multiprocessing import Pool,cpu_count
#from pymes.utils.parallel_processing import apply_async
#pool = Pool(10)

this_id = 0
for this_df,this_cpu_time in zip(data_frames,cpu_times):
#   apply_async(pool,write_conf_file, args = (this_df,input_file,folder,this_id))
   write_sbatch_file(this_df,input_file,folder,this_id,this_cpu_time)
   this_id += 1
