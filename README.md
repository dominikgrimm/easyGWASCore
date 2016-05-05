*********************
CEasyGWAS Framework
*********************

Core Framework of easyGAS (http://easygwas.ethz.ch) for computing genome-wide association studies and meta-analysis. This is a C/C++ Framework with Python interfaces. The code includes several standard methods for performing GWAS, such as linear regression, logistic regression and popular linear mixed models to also account for population stratification. In addition, the package contains code for the network guided multi-locus mapping method SConES (http://bioinformatics.oxfordjournals.org/content/29/13/i171.short).


Install Dependencies
------------------------

- Install SWIG: http://www.swig.org/download.html
- Install SCONS by typing in the following lines into the terminal
  $:> sudo pip install --egg scons



Compiling the Code
----------------------

Compile C++ API and Binaries Only
*********************************

To compile the code you have to go to the root directory of the Framework and type:

$>: scons -j 4 --build=release

This command builds the code using 4 CPUs and the release command uses optimization techniques to make the code faster.


Compile C++ API and Python Interfaces
**************************************

To compile the C/C++ code and build the python interfaces, SWIG has to be installed. To compile the interface you have to type:

$>: scons -j 4 --build=release --interface=python


*******
License
*******

Code by: Dominik Gerhard Grimm
Year: 2011-2016
Group: Machine Learning and Computational Biology Research Group
Insitute: Max Planck Institute for Intelligent Systems and Max Planck Institute for Developmental Biology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

