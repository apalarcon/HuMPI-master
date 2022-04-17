---
title: 'humpi: The python code for the Hurricane Maximum Potential Intensity (HuMPI) model'
tags:
  - Python
  - meteorology
  - programming  languages and software
  - atmospheric sciences
  - scientific programming
  - tropical cyclones maximum potential intensity
  - development software
authors:
  - name: Albenis Pérez-Alarcón
    orcid: 0000-0002-9454-2331
    affiliation: "1,2" 
  - name: José Carlos Fernández-Alvarez
    orcid: 0000-0003-3409-6138
    affiliation: "1,2"
  - name: Oscar Díaz-Rodríguez
    orcid: 0000-0002-9635-1184
    affiliation: "3"
affiliations:
 - name: Departamento de Meteorología, Instituto Superior de Tecnologías y Ciencias Aplicadas, Universidad de La Habana, La Habana, Cuba
   index: 1
 - name: Centro de Investigación Mariña, Universidade de Vigo, Environmental Physics Laboratory (EPhysLab), Campus As Lagoass/n, Ourense, 32004, Spain
   index: 2 
 - name: Centro de Ciencias de la Atmósfera, Universidad Autónoma de México, Ciudad de México, México
   index: 3 
date: 17 April 2022
bibliography: paper.bib
---
# Summary
The potential intensity (PI) of tropical cyclones (TCs) is the maximum surface wind speed and minimum central pressure limits found by representing the storm as a thermal heat engine [@Gilford2021,@Emanuel1986,@Hollan1997]. The PI theory proposed by @Emanuel1986 (hereafter E-PI) has been widely accepted as the upper bound for TCs intensity [@garner2015, @kieu2016, @kowaleski2016]. In the E-PI theory, the dynamic and thermodynamic processes of the TC are described as an energy cycle like a Carnot engine, absorbing heat from the ocean, giving it up at the tropopause. Nevertheless, the "superintensity" phenomenon, which occurs when the observed or modelled TC intensity is higher than the E-PI prediction, is a research challenge nowadays [@persing2003,@rousseau2019,@li2020]. In a recent attempt to avoid the "superintensity" phenomenon, @Perezalarcon2021 proposed a new hurricane maximum potential intensity (HuMPI) model based on the E-PI theory. HuMPI describes the TC thermo-energetic cycle as a generalized Carnot cycle and includes a TC model for the atmospheric boundary layer [@smith2003,@smith2008]. For further details of HuMPI physics description, see  @Perezalarcon2021.

@bister2002 coded the E-PI as a FORTRAN subroutine, while Kerry Emanuel later converted it for use as a MATLAB function. Despite the widespread use of the E-PI theory, the codes in FORTRAN and MATLAB have not been well documented [@Gilford2021]. Due to the advantages of the Python programming language, @Gilford2021 recently developed the Tropical Cyclone Potential Intensity Calculations in Python (piPy) to implement the E-PI theory. Therefore, this work aims to implement the HuMPI model formulation in Python.

# Statement of Need
The humpi Python package implements the HuMPI model for its extensive use in scientific research to understand the changes in TC intensity due to climate change.

# Python implementation
humpi (v1.0) is written in Python v3.8 and uses the mpi4py package for parallel runs. The humpi package requires netCDF4, numpy, scipy, mpi4py, os, time and datetime packages. Similar to piPy, the run times of humpi will depend on the user's particular implementation and computing resources. Computing the maximum intensity of TCs with humpi requires the sea surface temperature as input. Below we provided the basic commands for humpi usage

* <b>For help</b>
```
import humpi
humpi.help()
```
* <b>To get HuMPI input parameters templete file</b>
```
import humpi
humpi.get_HuMPI_inputs_template()
```

* <b>To get HuMPI input data file for multiple runs</b>
```
import humpi
humpi.get_HuMPI_input_data()
```

Additionally, you can use the basic implementation to run the HuMPI model, as indicated in the run_HuMPI.py script in the Github repository (https://github.com/apalarcon/HuMPI-master):
```
import humpi

args = humpi.read_args()
if args.HuMPI_help:
	humpi.help()
elif args.get_template:
	humpi.get_HuMPI_inputs_template()
elif args.get_input_data:
	humpi.get_HuMPI_input_data()
else:
	humpi.HuMPI_main(args.parameterfile) 
```

* <b>For help</b>
```
python run_HuMPI.py -hh t
```
* <b> For getting input paramters template</b>
```
python run_HuMPI.py -gt t
```
* <b> For getting input data file for multiples runs </b>
```
python run_HuMPI.py -id t
```
* <b> For running using MPI </b>
```
mpiexec -n N python run_HuMPI.py -pf input_paramters_file
```


# Acknowledgements
We acknowledge the supporting from Departamento de Meteorología, Instituto Superior de Tecnologías y Ciencias Aplicadas, Universidad de La Habana and Centro de Física de la Atmósfera del Instituto de Meteorología de Cuba.

# References
 
