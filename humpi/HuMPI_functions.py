from scipy.integrate import odeint 
from scipy.optimize import curve_fit
import numpy as np 
from  humpi.HuMPI_paramters import * 
from netCDF4 import Dataset,num2date,date2num
from datetime import datetime, timedelta
import os
import imp
import argparse
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning) 


@np.vectorize
def vfunc(t,l,n,ndic):
	if not np.isnan(t):
		return pi_function(t,l,n,ndic)
	return t,t


def program_name():
	return "HuMPI"

def program_fullname():
	return "Hurricane Maximum Potential Intensity model"

def get_currentversion():
	pathpkg = os.path.dirname(__file__)
	version_file = pathpkg+"/VERSION"
	with open(version_file) as vfile:
		version = vfile.readlines()[0].strip()
	return(version)

def disclaimer():
	print(
		"\n============================================================================================================"
	)
	print("||                                +     +           +       +  +++++++  +++++++                           ||")
	print("||                                +     +  +     +  + +   + +  +     +     +                              ||")
	print("||                                +++++++  +     +  +  +++  +  +++++++     +                              ||")
	print("||                                +     +  +     +  +   +   +  +           +                              ||")
	print("||                                +     +  +++++++  +   +   +  +        +++++++                           ||")
	print("||                           <--------------------------------------------------->                        ||")
	print("||                          " + program_fullname() +" Version " +str(get_currentversion())+"                       ||")
	print("||                                             Copyright 2022                                             ||")
	print("||                                                                                                        ||")
	print("||                          Albenis Pérez-Alarcón, José C. Fernández-Alvarez                              ||")
	print("||                                        Oscar Onoe Diaz Rodriguez                                       ||")
	print("||                                                                                                        ||")
	print("||                 " +program_name() + " Version " +str(get_currentversion())+ " is free under the terms of the GNU General Public license            ||")
	print("||                                         Departamento de Meteorologia                                   ||")
	print("||                         Instituto Superior de Tecnologias y Ciencias Aplicadas                         ||")
	print("||                                          Universidad de La Habana                                      ||")
	print("||                                                 Cuba                                                   ||")
	print("||                    contact: albenisp@instec.cu, jcfernandez@instec.cu, oodr72@gmail.com                ||")
	print("||                                                                                                        ||")
	print("============================================================================================================")

def ending_credits():
	print("\n============================================================================================================")
	print(program_name() +" Version " +str(get_currentversion()) + " has successfully finished")
	print("Bye :)")
	print("============================================================================================================")



def write_nc(latitude,longitude,vvar,pvar,date,filename="output"):
	ncout = Dataset(filename+".nc", 'w', format='NETCDF4')
	ncout.history='HuMPI outputs'
	ncout.title = "Tropical cyclones potential maximum intensity from " + program_name() + " model"
	today = datetime.now()
	ncout.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using "+program_name()+"."
	ncout.institution = (" Departamento de Meteorología, Instituto Superior de Tecnologías y Ciencias Aplicadas, Universidad de La Habana, Cuba")
	
	ncout.source = (
        program_name() + " version "
        + str(get_currentversion())
        + " ((c) Albenis Pérez-Alarcón; José C. Fernández-Alvarez; Oscar Diaz Rodriguez): albenisp@instec.cu, jcfernandez@instec.cu, oodr72@gmail.com"
    )
	
	# define axis size
	ncout.createDimension('time', 1)
	ncout.createDimension('nx', latitude.shape[1])
	ncout.createDimension('ny', latitude.shape[0])
	
	time = ncout.createVariable('time', 'f8', ('time',))
	time.long_name = 'time'
	time.units = 'seconds since 1970-01-01 00:00:00'
	time.axis = 'T'
	time.calendar = "Standard"

	# create latitude axis
	lat = ncout.createVariable('lat', 'f4', ('ny','nx'))
	lat.standard_name = 'latitude'
	lat.long_name = 'latitude)'
	lat.units = 'degrees'
	lat.axis = 'Y'

	# create longitude axis
	lon = ncout.createVariable('lon', 'f4', ('ny','nx'))
	lon.standard_name = 'longitude'
	lon.long_name = 'longitude'
	lon.units = 'degrees'
	lon.axis = 'X'

	# create variable array
	vout = ncout.createVariable('vmax', 'f4', ('ny','nx'),zlib=True)
	vout.long_name = 'vmax'
	vout.units = 'm/s'
	vout.standard_name = "Potential maximum wind speed";
	vout.coordinates = "ny,nx" ;
	vout.original_name = "vmax"
	
	vout1 = ncout.createVariable('pmin', 'f4', ('ny','nx'),zlib=True)
	vout1.long_name = 'pmin'
	vout1.units = 'hPa'
	vout1.standard_name = "Potential minimum central pressure";
	vout1.coordinates = "ny,nx" ;
	vout1.original_name = "pmin"
	
	dates=datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]),0,0,0)
	time[:]=date2num(dates, units=time.units,calendar=time.calendar)
	
	lon[:] = longitude[:]
	lat[:]= latitude[:]
	vout[:] = vvar[:]
	vout1[:] = pvar[:]
	ncout.close()
##################################################################################

def print_error_message(message):
	print("=============================================================================================================")
	print ("ERROR: "+message)
	print("=============================================================================================================")
	raise SystemExit("Bye :)")


def check_paths(pfile, path):
	try: 
		fpath = getattr(pfile, path)
	except:
		fpath = ""
	return fpath

def check_parameters_files(pathfile):
	if os.path.exists(pathfile):
		pass
	else:
		print_error_message(" --parameterfile or -pf are missing or "+ pathfile + " is not in the directory")

def str2boolean(arg):
    if isinstance(arg, bool):
        return arg
    if arg.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif arg.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")
	
def read_args():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"--parameterfile",
		"-pf",
		help="name of parameters file. To use "+program_name()+" run: python run_"+program_name()+".py -pf <your_input_file>",
		metavar="",
		type=str,
		default="humpi_inputs.cfg",
	)
	parser.add_argument(
		"--HuMPI_help",
		"-hh",
		help="Help for " + program_name() + " input parameters. Run: python run_"+program_name()+".py -hh t or python run_"+program_name()+".py -humpi_help t",
		metavar="",
		type=str2boolean,
		default=False,
	)
	parser.add_argument(
		"--get_template",
		"-gt",
		help="Get a template for the " + program_name() + " input parameters. Run: python run_"+program_name()+".py -gt t or python run_"+program_name()+".py --get_template t",
		metavar="",
		type=str2boolean,
		default=False,
	)
	
	parser.add_argument(
		"--get_input_data",
		"-id",
		help="Get input_data file template for " + program_name() + " multiples runs. Run: python run_"+program_name()+".py -id t or python run_"+program_name()+".py --get_input_data t",
		metavar="",
		type=str2boolean,
		default=False,
	)
	
	args = parser.parse_args()
	
	return args

def get_HuMPI_inputs_template():
	disclaimer()
	wpath=os.getcwd()
	pathpkg = os.path.dirname(__file__)
	os.system("cp "+pathpkg+"/humpi_input.cfg " + wpath)
	print("\n============================================================================================================")
	print("humpi_inputs template was successfully copied to thr working directory: " + wpath)
	print("Bye :)")
	print("============================================================================================================")

def get_HuMPI_input_data():
	disclaimer()
	wpath=os.getcwd()
	pathpkg = os.path.dirname(__file__)
	os.system("cp "+pathpkg+"/humpi_input_data.cfg " + wpath)
	print("\n============================================================================================================")
	print("humpi_input_data template was successfully copied to thr working directory: " + wpath)
	print("Bye :)")
	print("============================================================================================================")

def help():
	disclaimer()
	print("")
	print(program_name() + " HELP FOR INPUT PARAMETERS")
	print("------------------------------------------------------------------------------------------------------------")
	print("+ date              = 'yyyymmdd'                  -> date")
	print("+ input_path        = 'path'                      -> Path to input data")
	print("+ input_format      = 'netcdf' or 'txt'           -> Data input format")
	print("+ sst_netcdf_file   = 'sst_filename'              -> SST netcdf filename. E.g. sst.nc")
	print("+ sst_var           = 'sst'                       -> SST var name.")
	print("+ sst_lon           = 'lon'                       -> longitude var name.")
	print("+ sst_lat           = 'lat'                       -> latitude var name.")
	print("+ sst_units         = 'C' or 'K'                  -> Units for SST")
	print("+ i_terrain_mask    = 0                           -> row position of masked value for terrain")
	print("+ j_terrain_mask    = 0                           -> column position of masked value for terrain")
	print("+ output_path       = 'path'                      -> Path for " + program_name() + "output")
	print("+ input_data        = path_to_input_data_file     -> Paramters for multiples runs")
	print("")
	print("NOTES")
	print("--------------------")
	print("* If input_format='txt' set sst_var as full txt filename for sst. E.g. sst.txt")
	print("* If input_format='txt' set sst_lat as full txt filename for latitudes. E.g. lat.txt")
	print("* If input_format='txt' set sst_lon as full txt filename for longitudes. E.g. lon.txt")
	print("* The input data file must be follow the structure below")
	print("   input_format,  sst_filename,     sst_var,      sst_lat,     sst_lon,   sst_units,  date,     i_terrain_mask,   j_terrain_mask ")
	print("")
	print("* To get a template for "+program_name() +" input parameters  run_"+program_name()+".py -gt t")
	print("* To get a template for "+program_name() +" input data file   run_"+program_name()+".py -id t")
	
def get_input_data(input_data=""):
	input_d=open(input_data)
	input_d=input_d.readlines()
	dates=[]
	input_formats=[]
	sst_netcdf_files=[]
	sst_vars=[]
	sst_lons=[]
	sst_lats=[]
	sst_unitss=[]
	i_terrain_masks=[]
	j_terrain_masks=[]
	
	for i in range(0,len(input_d)):
		if input_d[i][0]=="#" or input_d[i][0]==" " or input_d[i][0]=="\n":
			pass
		else:
			line=input_d[i].split(",")
			input_format=line[0].replace(" ", "")
			sst_netcdf_file=line[1].replace(" ", "")
			sst_var=line[2].replace(" ", "")
			sst_lat=line[3].replace(" ", "")
			sst_lon=line[4].replace(" ", "")
			sst_units=line[5].replace(" ", "")
			date=line[6].replace(" ", "")
			i_terrain_mask=int(line[7])
			j_terrain_mask=int(line[8].split("\n")[0])
			
			dates=np.append(dates,date)
			input_formats=np.append(input_formats,input_format)
			sst_netcdf_files=np.append(sst_netcdf_files,sst_netcdf_file)
			sst_vars=np.append(sst_vars,sst_var)
			sst_lons=np.append(sst_lons,sst_lon)
			sst_lats=np.append(sst_lats,sst_lat)
			sst_unitss=np.append(sst_unitss,sst_units)
			i_terrain_masks=np.append(i_terrain_masks,i_terrain_mask)
			j_terrain_masks=np.append(j_terrain_masks,j_terrain_mask)
	
	return dates,input_formats, sst_netcdf_files, sst_vars, sst_lons, sst_lats, sst_unitss, i_terrain_masks, j_terrain_masks
		
		
		
	
def get_values(input_path="./",
				input_format="netcdf",
				sst_netcdf_file="file",
				sst_var="sst",
				sst_lon="lon",
				sst_lat="lat",
				sst_units="C",
				i_terrain_mask=0,
				j_terrain_mask=0):
	
	if input_format.lower()=="netcdf":
		check_iffile_exist(input_path+"/"+sst_netcdf_file)
		nc=Dataset(input_path+"/"+sst_netcdf_file)
		check_netcdfvar(nc,sst_var,sst_netcdf_file)
		check_netcdfvar(nc,sst_lat,sst_netcdf_file)
		check_netcdfvar(nc,sst_lon,sst_netcdf_file)
		
		sst=nc.variables[sst_var][:]
		lat=nc.variables[sst_lat][:]
		lon=nc.variables[sst_lon][:]
		
		if len(sst.shape)>2:
			sst_len=len(sst.shape)
			while sst_len>2:
				sst=sst[0,:]
				sst_len=len(sst.shape)
	
	elif input_format.lower() in ('txt','dat'):
		check_iffile_exist(input_path+"/"+sst_var)
		check_iffile_exist(input_path+"/"+sst_lat)
		check_iffile_exist(input_path+"/"+sst_lon)
		
		sst=np.loadtxt(input_path+"/"+sst_var)
		lat=np.loadtxt(input_path+"/"+sst_lat)
		lon=np.loadtxt(input_path+"/"+sst_lon)
		
	else:
		print_error_message("input format "+ input_format + " is unrecognized")
	
	if len(lat.shape)==1:
		lon,lat=np.meshgrid(lon,lat)
	
	sst[sst==sst[int(i_terrain_mask),int(j_terrain_mask)]]=np.nan




	if sst_units.lower() in ("celsius", "c"):
		sst=sst+273.15

	return lat,lon,sst
		
def convert_lons(lons):
	for i in range(0,lons.shape[0]):
		for j in range(0,lons.shape[1]):
			if lons[i,j]>=180:
				lons[i,j]=lons[i,j]-360
	
	return lons
	
	
def check_netcdfvar(nc_fid,key,ncfile):
	try:
		for ncattr in nc_fid.variables[key].ncattrs():
			pass
	except KeyError:
		print_error_message(key + " does not contain variable attributes in --> " + ncfile +" file")	

def check_iffile_exist(pathfile):
	if os.path.exists(pathfile):
		pass
	else:
		print_error_message( pathfile+ " is not in the directory")

	
def evgrh(r,rm,rho,pc,pa,fcor,x,b):
	vgr=-fcor*r/2 + ( (fcor*r/2)**2+r/rho*(x*rm/r**2 * (pa-pc)*np.exp(x*(b - rm/r))) )**1/2
	return vgr/100
	

def rvgrh(r,rm,rho,pc,pa,fcor,B):
	vgr=np.sqrt(rm**B*B*(pa-pc)*np.exp(-(rm/r)**B)/(rho*r**B)+(fcor*r)**2/4)-fcor/2*r
	return vgr
	
def pholland(r,rm,pc,pa,B):
	p = pc+(pa-pc)*np.exp(-(rm/r)**B)
	return p
  
def pi_emanuel(ts,tc,ta,t0,pc,pa,fcor,ebl,ck,cd):
	b=1.
	print ('calculo de las humedades especificas') 
	qc=(3.802/pc)*np.exp(17.67*(tc-273.15)/(243.5 + tc-273.15))*100
	qa=(3.802/pa)*np.exp(17.67*(ta-273.15)/(243.5 + ta-273.15))*100

	ha= cp*ta + 9.8*ebl + lv*qa
	hc= cp*tc + 9.8*ebl + lv*qc
	vmax= np.sqrt((ck/cd)*((((ts-t0)/(ts))*(hc-ha)-(1./4.)*(fcor**2)*(ra**2))/(1. -(ck/cd)*((0.5)-b))))
	return vmax
  
def pi_emanuel_sustitucion(ts,tc,ta,t0,pc,pa,fcor,ebl,ck,cd):
	qc=(3.802/pc)*np.exp(17.67*(tc-273.15)/(243.5 + tc-273.15))*100
	qa=(3.802/pa)*np.exp(17.67*(ta-273.15)/(243.5 + ta-273.15))*100
	ha= cp*ta + 9.8*ebl + lv*qa
	hc= cp*tc + 9.8*ebl + lv*qc
	print ('ha',ha)
	print ('hc',hc)

	vmax= np.sqrt((ck/cd)*((((ts-t0)/(t0))*(hc-ha)-(1./4.)*(fcor**2)*(ra**2))/(1 -(ck/cd)*(ts/(2*t0)-1))))
	return vmax
  
  
def pi_emanuel_new(fcor,rm,rg,ra,ts,tb,tc,ta,t00,pc,pa,rho,ck,cd):
	qc=(3.802/pc)*np.exp(17.67*(tc-273.15)/(243.5 + tc-273.15))*100
	qa=(3.802/pa)*np.exp(17.67*(ta-273.15)/(243.5 + ta-273.15))*100
	ha= cp*ta + 9.8*ebl + lv*qa
	hc= cp*tc + 9.8*ebl + lv*qc
	pm=pholland(rm,rm,psc,psa,B)
	alpha1 = 1 - (1./2.)*(ck/cd)*(tb/t00)
	alpha2 = -0.5*(ck/cd)*fcor*rm*(tb/t00)
	alpha3 = (-(tb-t00)/t00)*(ck/cd)*(hc-ha) - (ck/cd)*tb*Rd*np.log(pm/psa) + 0.25*(ck/cd)*(fcor**2)*(ra**2)*(tb/t00)
	if alpha2**2 - 4*alpha1*alpha3<0:
		vmax=-9999
	else:
		vmax = (-alpha2 + np.sqrt(alpha2**2 - 4*alpha1*alpha3))/(2*alpha1)
	
	return vmax
  
  
def recalculate_pc_ciclostrofic(pa,rho,vmax,B):
	pc= pa - (vmax**2)*(1./B * rho * np.exp(1))
	return pc
	

def vgrholland(r,rm,pc,pa,rho,fcor,B):
	vgr=np.sqrt(rm**B*B*(pa-pc)*np.exp(-(rm/r)**B)/(rho*r**B)+(fcor*r)**2/4)-fcor/2*r
	return vgr
 

def pres_holland(r,rm,rg,rho,fcor,pc,pa):
	pm = pc+(pa-pc)*np.exp(-(rm/r)**B)
	return pm
	
  
def ph(r):
	p = pc+(pg-pc)*np.exp(-(rm/r)**B)
	return p

  
  
def vgrh(r):
	vgr=np.sqrt(rm**B*B*(pg-pc)*np.exp(-(rm/r)**B)/(rho*r**B)+(fcor*r)**2/4)-fcor/2*r
	return vgr
	
  

def cdragblack(vnet):
	
	cd0=0.7e-3
	cd1=6.5e-5
	
	if vnet <= 20:
		cd=cd0+cd1*vnet
	else:
		cd=2.0e-3
	return cd


def cdragline(vnet):
	cd0=0.7e-3
	cd1=6.5e-5
	cd=cd0+cd1*vnet
	return cd

def cdraggarrant(vnet):
	if vnet <= 4:
		cd=1.0e-3
	else:
		cd=(0.75+0.067*vnet)*1.0e-3
	return cd

r = np.arange(re,r0,rstep)

for i, ind in enumerate(r):
	if ind==rm-15000: r1=r[i]
	if ind==rm+15000: r2=r[i]
   

def pi_function(ts, lat,nudic,ndic_vmax):
	from  humpi.HuMPI_paramters import   psc,depth,n_pres,r0,re,rm,ra,rg,psg,psa,B,t0,t00,nciclo,tolerancia,ebl,wh,wsc,whneg,pressprofile,cd,ck,cdparametrization,windprofile,depthmin,vst,pst,depth_correction,omega,lv,cp,Rd ,rho,gammha,sigmha


	lat=np.abs(lat)
	fcor= 2*omega*np.sin(lat*pi/180)
	tc= ts - sigmha*ebl
	ta= ts - gammha*ebl
	tb= ts -  sigmha*ebl
  
	if sigmha == 0:
		pc=psc
	else:
		pc=psc*(1+sigmha*ebl/ts)**(-grav/(sigmha*Rd)) 

	if gammha == 0:
		pa=psa
		pg=psg
	else:
		pa=psa*(1+gammha*ebl/ts)**(-grav/(gammha*Rd)) 
		pg=psg*(1+gammha*ebl/ts)**(-grav/(gammha*Rd)) 


	def rvg_pc_B(r,pc,B):   
		vgr=np.sqrt(rm**B*B*(pa-pc)*np.exp(-(rm/r)**B)/(rho*r**B)+(fcor*r)**2/4)-fcor/2*r
		return vgr


	def rvg_pc(r,pc):
		vgr=np.sqrt(rm**B*B*(pa-pc)*np.exp(-(rm/r)**B)/(rho*r**B)+(fcor*r)**2/4)-fcor/2*r
		return vgr

	pres=pres_holland 
  
 
	dic_vmax={1:pi_emanuel,2:pi_emanuel_sustitucion,3:pi_emanuel_new}
	if ndic_vmax==0 or ndic_vmax==1 or ndic_vmax==2: 
		vmax=dic_vmax[ndic_vmax](ts,tc,ta,t0,pc,pa,fcor,ebl,ck,cd)
	if ndic_vmax==3:
		vmax=dic_vmax[ndic_vmax](fcor,rm,rg,ra,ts,tb,tc,ta,t00,pc,pa,rho,ck,cd)
	if vmax>0:
		rm=(46.6*np.exp(-0.015*vmax+0.0169*lat))*1000.
	else:
		rm=-9999
	 
	x2=45
	def vwill(r):
   
		x1= 287.5 - 1.942*vmax + 7.799*np.log(rm/1000.) +  1.819*lat
		n=2.1340 + 0.0077*vmax - 0.4522*np.log(rm/1000.) - 0.0038*lat
		A = 0.5913 + 0.0029*vmax-0.1361*np.log(rm/1000.) - 0.0042*lat
   
		if r >= rm:
			vw= vmax*((1-A)*np.exp((-r/1000.+rm/1000.)/x1) + A*np.exp((-r/1000.+rm/1000.)/x2))
		else:
			vw= vmax*((r/rm)**n)
		return vw



	def w(e):
		return 126.0*e**5-420.0*e**6+540.0*e**7-315.0*e**8+70.0*e**9

	def E(r, r1, r2): 
		return (r-r1)/(r2-r1)

	def vwill_new(r):
		rmax=rm/1000
		r1=rmax-0.0001
		r2=rmax+0.0001
		latt=lat
		x2=45
		x1= 287.6 - 1.942*vmax + 7.799*np.log(rmax)+1.819*np.abs(latt)
		n=2.1340 + 0.0077*vmax - 0.4522*np.log(rmax)-0.0038*np.abs(latt)
		A = 0.5913 + 0.0029*vmax-0.1361*np.log(rmax)-0.0042*np.abs(latt)
		if  r <= r1: 
			vgr=vmax*(r/rm)**n
		elif  r1 < r <= r2: 
			Vi = vmax*(r1/rm)**n
			V0 = vmax*((1-A)*np.exp(-(r2/1000-rmax)/x1))+A*np.exp(-(r2/1000-rmax)/x2)
			vgr= Vi*(1-w(E(r/1000, r1, r2))) + V0*w(E(r, r1, r2))
		else:
			vgr=vmax*((1-A)*np.exp((-r/1000.+rm/1000.)/x1) + A*np.exp((-r/1000.+rm/1000.)/x2))

		return vgr


	if windprofile==1:
		vperfil=vwill  
	elif windprofile==2:
		vperfil=vgrh  
	elif windprofile==3:
		vperfil=vwill_new 


	def recalculated_perss_0():
		pc = recalculate_pc_ciclostrofic(pa,rho,vmax,B)
		return pc
  
	def recalculated_perss_1(pc):
		pm=pholland(rm,rm,pc,pa,B)
		pc =pm*np.exp(-0.5*(vmax**2+vmax*fcor*rm)/(cp*tc))
		return pc
  
	def recalculated_perss_3(pc):
		pm=pres(rm,rm,rg,rho,fcor,pc,pa)
   
		n=2.1340 + 0.0077*vmax - 0.4522*np.log(rm/1000.) - 0.0038*lat
		pc = pm - rho*vmax*(0.5*vmax/n + fcor*rm/(n+1))
		return pc
    
	def recalculated_perss_4():
		perfilv=np.zeros(np.shape(r))
		for index, radio in enumerate(r):
			perfilv[index]=vperfil(radio)
		(pc1, B1), pcov = curve_fit(rvg_pc_B, r, perfilv,method='trf', bounds=([80000,1], [pa, 2.8]))
		B = (vmax**2+vmax*rm*fcor)*(rho * 2.71828182846 )/(pa-pc1)
		return pc1 
  
	def recalculated_perss_5():
		perfilv=np.zeros(np.shape(r))
		for index, radio in enumerate(r):
			perfilv[index]=vperfil(radio)
		(pc1, B), pcov = curve_fit(rvg_pc_B, r, perfilv,method='trf', bounds=([80000,1], [pa, 2.8]))
		pm=pholland(rm,rm,pc1,pa,B)
		pc =pm*np.exp(-0.5*(vmax**2+vmax*fcor*rm)/(cp*tc))
		return pc
	
	def recalculated_perss_6(pa):
 
		def three(y,r,vmax,rm,rho,fcor,r1,r2):
			latt=lat
			n=2.1340+0.0077*vmax-0.4522*np.log(rm)-0.0038*latt
			x1=287.6-1.942+7.799*np.log(rm)+1.819*latt
			x2=45.0
			A=0.5913+0.0029*vmax-0.1361*np.log(rm)-0.0042*latt
		
			y2 =(rho/r)*((vmax*(1-A)*np.exp(-((r-rm)/x1))+A*np.exp(-((r-rm)/x2))+(fcor*r)/2)**2)-(1./4.)*rho*r*(fcor**2)
			return np.array([y2])
		ic1 = np.array([pa])
	
		inte1 = np.arange(re/1000.,rm/1000.,-1.)
	
		
		sol1 = odeint(three, ic1, inte1,args=(vmax,rm/1000.,rho,fcor,r1/1000.,r2/1000.),rtol=1e-4, atol=1e-4)
		y11 = sol1[:,0]
		ps=y11[-1]

	
  
		global pn 
		def solAnalitica(r,vmax,rm,rho,fcor,pn=ps):
			latt=lat
			n=2.1340+0.0077*vmax-0.4522*np.log(rm)-0.0038*latt
			a=pn+(rho*(vmax**2)/(rm**(2*n)))*((r**(2*n))/2*n)-(rho*(vmax**2)/(rm**(2*n)))*((rm**(2*n))/2*n)+rho*fcor*(vmax/(rm**n))*(r**(n+1))/(n+1)
			b=rho*fcor*(vmax/(rm**n))*(rm**(n+1))/(n+1)
			c=a-b
			return c

		inte3 = np.arange(r0/1000.,rm/1000.,1)
		a=solAnalitica(inte3,vmax,rm/1000.,rho,fcor,pn=ps)
		
		a=list(a)
		a.append(ps)
		a=np.array(a)
		inte3=list(inte3)
		inte3.append(rm/1000.)
		inte3=np.array(inte3)
		pc=a[0]
		return pc

	def recalculated_perss_7(pc):
		pm=pres(rm,rm,rg,rho,fcor,pc,pa)
		pc=pm*np.exp(-(vmax**2)/(2*cp*ts))
		return pc
	def recalculated_perss_8():
		pc = recalculate_pc_ciclostrofic1(pa,rho,vmax,B)
		return pc
	dic={0:recalculated_perss_0,1:recalculated_perss_1,3:recalculated_perss_3,4:recalculated_perss_4,5:recalculated_perss_5,6:recalculated_perss_6,7:recalculated_perss_7,8:recalculated_perss_8}

	cont=1
	error=100
 
	while (error >= tolerancia and cont <= nciclo) and vmax>0:
		
		
		if nudic==0 or nudic==4 or nudic==5 or nudic==8: pc=dic[nudic]()
		if nudic==1 or nudic==3 or nudic==7: pc=dic[nudic](pc)
		if nudic==6: pc=dic[nudic](pa)
   
		if sigmha == 0:
			psc=pc
		else:
			psc=pc/(1+sigmha*ebl/ts)**(-grav/(sigmha*Rd))
    
		if ndic_vmax==0 or ndic_vmax==1 or ndic_vmax==2:
			vmax_rec=dic_vmax[ndic_vmax](ts,tc,ta,t0,pc,pa,fcor,ebl,ck,cd)
		if ndic_vmax==3:
			vmax_rec=dic_vmax[ndic_vmax](fcor,rm,rg,ra,ts,tb,tc,ta,t00,pc,pa,rho,ck,cd)

		error=abs((vmax-vmax_rec)/vmax)
		vmax=vmax_rec 
		cont+=1

	nciclo= cont
 

	def system(y,r):

		u,v=y[0],y[1] 
		vgr=vperfil(r)
		vnet=np.sqrt(u**2+v**2)
		cd=cdrag(vnet)
		friction=-cd*vnet/ebl
		uinv=(1./u)
		vvg=(v-vgr)
		su= uinv*((vgr**2-v**2)/r - fcor*vvg - friction*u)
		sw= su - u/r
		if (sw > 0):
			alfa=1
		else:
			alfa=0
		wh=(ebl/(1+alfa))*sw 
		whneg=(wh-abs(wh))/2.
		wv=(whneg+wsc)/ebl
		dudr=  wv - su
		dvdr=  uinv*(wv*vvg - (v/r +fcor)*u + friction*v)
		return [dudr,dvdr]


	def whtop(r,u,v):
		vgr=vperfil(r)
		sw=(1./u)*((vgr**2-v**2)/r + fcor*(vgr-v)+ (cd/ebl)*np.sqrt(u**2+v**2)*u) -u/r
		if (sw > 0):
			alfa=1
		else:
			alfa=0
		wh=(ebl/(1+alfa))*sw
		return wh


	if cdparametrization==1:
		cdrag=cdragblack 
	elif cdparametrization==2:
		cdrag=cdraggarrant   
	elif cdparametrization==3:	
		cdrag=cdragline 

	if vmax>0:
		vgr=vperfil(re) 
		cd=cdrag(vgr) 
		a=(ebl*fcor)/(cd*vgr)
		v1=-0.5*a**2 + np.sqrt(0.25*a**4 + a**2)
		u1=-np.sqrt((1-v1)*v1)
		u0=u1*vgr 
		v0=v1*vgr

		sol=odeint(system,[u0,v0],r) 
		u,v=sol[:,0],sol[:,1] 
		vnet=np.sqrt(u**2+v**2) 
		vmaxfinal=max(vnet)

		vmaxfinal=vmaxfinal-vmaxfinal*(25./100.)

		def correct_depth(vmax,vst,psc,pst,depth,depthmin):
			if abs(depth) < abs(depthmin):
				if vmax > vst:
					vmax=vst+(vmax-vst)*abs(depth/depthmin)
				if psc < pst:
					psc=pst-(pst-psc)*abs(depth/depthmin)
			return vmax,psc

		if depth_correction == 1:
			vmaxfinal,pscfinal= correct_depth(vmaxfinal,vst,psc,pst,depth,depthmin)
	else:
		vmaxfinal=-9999
   
	if vmaxfinal<0:
		pscfinal,vmaxfinal=np.nan,np.nan

	return (pscfinal,vmaxfinal)


