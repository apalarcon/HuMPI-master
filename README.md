# HuMPI: Hurricane Maximum Potential Intensity model
The python code for the Hurricane Maximum Potential Intensity model

[![Anaconda-Server Badge](https://anaconda.org/apa_py/humpi/badges/version.svg)](https://anaconda.org/apa_py/humpi) [![Anaconda-Server Badge](https://anaconda.org/apa_py/humpi/badges/license.svg)](https://anaconda.org/apa_py/humpi)

For HuMPI model descripction, see

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

# Installation

We recommend using the HuMPI package under the Anaconda environment.

To install humpi use  
```
conda install -c apa_py humpi
```

# Usage
<b>For HuMPI input parameters use</b>:
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
  
# Contact
- Albenis Pérez Alarcón: apalarcon1991[a]gmail.com

- José Carlos Fernández Alvarez: fortunajcfa[a]gmail.com

- Oscar Díaz Rodríguez: oodr71[a]gmail.com


