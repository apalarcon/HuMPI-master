from  humpi.HuMPI_paramters import * 
from scipy.integrate import odeint 
from scipy.optimize import curve_fit 
import numpy as np 
import scipy.optimize as opt
import random
import sys, time
from netCDF4 import Dataset
from humpi.HuMPI_functions import * 
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
import imp
import argparse
import warnings
import functools
print = functools.partial(print, flush=True)

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning) 


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def HuMPI_main(pathfile):

	if rank==0:
		disclaimer()
	
		print("\n")
		check_parameters_files(pathfile)
	
		print("---> | Using parameters from: " + pathfile)
		print("       "+program_name() +" will run using " + str(int(size)) + " CPUs")
		print("============================================================================================================\n")
	
	content = imp.load_source("", pathfile) 
	dates = check_paths(content, "date")
	input_path = check_paths(content, "input_path")
	input_formats = check_paths(content, "input_format")
	sst_netcdf_files = check_paths(content, "sst_netcdf_file")
	sst_vars = check_paths(content, "sst_var")
	sst_lons = check_paths(content, "sst_lon")
	sst_lats = check_paths(content, "sst_lat")
	sst_unitss = check_paths(content, "sst_units")
	i_terrain_masks = check_paths(content, "i_terrain_mask")
	j_terrain_masks = check_paths(content, "j_terrain_mask")
	output_path = check_paths(content, "output_path")
	input_data = check_paths(content, "input_data")
	

	if input_data=="":
		dates=[dates]
		input_formats=[input_formats]
		sst_netcdf_files=[sst_netcdf_files]
		sst_vars=[sst_vars]
		sst_lons=[sst_lons]
		sst_lats=[sst_lats]
		sst_unitss=[sst_unitss]
		i_terrain_masks=[i_terrain_masks]
		j_terrain_masks=[j_terrain_masks]
	else:
		dates,input_formats, sst_netcdf_files, sst_vars, sst_lons, sst_lats, sst_unitss, i_terrain_masks, j_terrain_masks= get_input_data(input_data=input_data)
	
	if rank==0:
			print ("\n---> | Computing the tropical cyclones maximum potential intensity")
			print("============================================================================================================")
	for index in range(0,len(dates)):
		date=dates[index]
		input_format=input_formats[index]
		sst_netcdf_file=sst_netcdf_files[index]
		sst_var=sst_vars[index]
		sst_lon=sst_lons[index]
		sst_lat=sst_lats[index]
		sst_units=sst_unitss[index]
		i_terrain_mask=i_terrain_masks[index]
		j_terrain_mask=j_terrain_masks[index]
		
		latitude,longitude,temp=get_values(input_path=input_path,
											input_format=input_format,
											sst_netcdf_file=sst_netcdf_file,
											sst_var=sst_var,
											sst_lon=sst_lon,
											sst_lat=sst_lat,
											sst_units=sst_units,
											i_terrain_mask=i_terrain_mask,
											j_terrain_mask=j_terrain_mask)
		

		if rank==0:
			if input_format.lower()=="netcdf":
				procfile=sst_netcdf_file
			else:
				procfile=sst_var
			print ("     + Processing --> | "+input_path + "/"+procfile)

		lon=np.copy(longitude)
		lat=np.copy(latitude)
		
		lon=convert_lons(lon)
		
		sendcounts=np.zeros(size)
		displacements=np.zeros(size)
		
		a,b=temp.shape
		div=(a*b)/int(size)
		rest=(a*b)%int(size)

		temp=temp.reshape(1,a*b)
		lat=lat.reshape(1,a*b)
		sendcounts[:]=int(div)
		sendcounts[0]=int(div+rest)
		

		displacements= np.insert(np.cumsum(sendcounts), 0, 0)[0:-1]
		

		temp=np.array(temp,np.float64)
		lat=np.array(lat,np.float64)

		if rank==0:
			output_temp = np.zeros(int(div+rest))
			output_lat = np.zeros(int(div+rest))
		else: 
			output_temp = np.zeros(int(div))
			output_lat = np.zeros(int(div))
		

		output_lat=np.array(output_lat,np.float64)
		output_temp=np.array(output_temp,np.float64)

		comm.Scatterv([temp,sendcounts, displacements,MPI.DOUBLE],[output_temp,MPI.DOUBLE], root=0)
		comm.Scatterv([lat,sendcounts, displacements,MPI.DOUBLE],[output_lat,MPI.DOUBLE], root=0)

		comm.Barrier()
		##########################################################################################
		
		

		for ndic in n_pres:
			comm.Barrier()
			if rank==0: start_time = MPI.Wtime()
			pscfinal,vmaxfinal=vfunc(output_temp,output_lat,ndic,ndic_vmax)
				
		p_min=np.zeros(a*b)
		v_max=np.zeros(a*b)
		
		comm.Gatherv(pscfinal,[p_min,sendcounts,displacements,MPI.DOUBLE], root=0)
		comm.Gatherv(vmaxfinal,[v_max,sendcounts,displacements,MPI.DOUBLE], root=0)
		comm.Barrier()
		p_min=p_min.reshape(a,b)
		v_max=v_max.reshape(a,b)

		if rank==0:
			end_time = MPI.Wtime()
			runtime = end_time - start_time
			runtime=time.strftime("%H:%M:%S", time.gmtime(runtime))
			
			
			print ("     --> Saving the tropical cyclones maximum potential intensity to a netcdf file\n")
			
			write_nc(latitude,longitude,v_max,p_min/100,date,filename=output_path+"/"+program_name()+"out_"+date)
		
	if rank==0:
		print("\n============================================================================================================")
		print ("\n---> | "+program_name() + " run time: ", runtime)
		
		
		ending_credits()


