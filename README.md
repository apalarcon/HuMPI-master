# HuMPI: Hurricane Maximum Potential Intensity model
The python code for the Hurricane Maximum Potential Intensity model

[![Anaconda-Server Badge](https://anaconda.org/apa_py/humpi/badges/version.svg)](https://anaconda.org/apa_py/humpi)

[![Anaconda-Server Badge](https://anaconda.org/apa_py/humpi/badges/license.svg)](https://anaconda.org/apa_py/humpi)

For HuMPI model physics description, see

Pérez-Alarcón, A.; Fernández-Alvarez, J.C.; Díaz-Rodríguez, O; (2021). Hurricane Maximum Potential Intensity Model, Revista Cubana de Física, 38(2), 85-93, http://www.revistacubanadefisica.org/index.php/rcf/article/view/2021v38p085



# HuMPI requirements


<table>
<thead>
<tr>
<th>python</th>
<th>status</th>
</tr>
</thead>
<tbody>
<tr>
<td>Python2</td>
<td> tested</td>
</tr>
<tr>
<td>Python3</td>
<td> tested</td>
</tr>
</tbody>
</table>

*  netCDF4
*  numpy 
*  scipy 
*  mpi4py

# Installation

We recommend using the HuMPI package under the Anaconda environment.

[![Anaconda-Server Badge](https://anaconda.org/apa_py/humpi/badges/installer/conda.svg)](https://conda.anaconda.org/apa_py)
```
conda install -c apa_py humpi
```
Using PIP
```
pip install humpi
```

# Usage
<b>For HuMPI help</b>:
```
import humpi
humpi.help()
```

<b> To get HuMPI input parameters templete file</b>:
```
import humpi
humpi.get_HuMPI_inputs_template()
```

<b> To get HuMPI input data file for multiple runs</b>:
```
import humpi
humpi.get_HuMPI_input_data()
```

<b> To run HuMPI model we recomend the using of run_HuMPI.py script as follow:</b>
 * For help
```
python run_HuMPI.py -hh t
```
 * For getting input paramters template
```
python run_HuMPI.py -gt t
```
 * For getting input data file for multiples runs
```
python run_HuMPI.py -id t
```
 * For running using MPI
```
mpiexec -n N python run_HuMPI.py -pf input_paramters_file
```
 
# Input data
Sea sruface temperature in netcdf or ascii format. In both cases the latitude and longitude arrays must be provided

# HuMPI outputs
A netcdf file, which contains the maximum potential surface wind speed and minimum central pressure for tropical cyclones 
 
# Contact and Support
- Albenis Pérez Alarcón: apalarcon1991[a]gmail.com

- José Carlos Fernández Alvarez: fortunajcfa[a]gmail.com

- Oscar Díaz Rodríguez: oodr71[a]gmail.com

### Referencing
If you use HuMPI, please cite:
Pérez-Alarcón, A.; Fernández-Alvarez, J.C.; Díaz-Rodríguez, O; (2021). Hurricane Maximum Potential Intensity Model, Revista Cubana de Física, 38(2), 85-93, http://www.revistacubanadefisica.org/index.php/rcf/article/view/2021v38p085

### License
Copyright 2022 Albenis Pérez-Alarcón, José C. Fernández-Alvarez, Oscar Díaz Rodríguez. 

This software is published under the GPLv3 license. This means: 
1. Anyone can copy, modify and distribute this software. 
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. If you dare build your business solely from this code, you risk open-sourcing the whole code base.
6. If you modify it, you have to indicate changes made to the code.
7. Any modifications of this code base MUST be distributed with the same license, GPLv3.
8. This software is provided without warranty.
9. The software author or license can not be held liable for any damages inflicted by the software.

