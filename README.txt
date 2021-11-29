
Before running any of the files, run the script 'create_mex_functions.m'. This only needs to be done once.

The following scripts compute volumes of some interesting shapes along with diagnostic plots
1. calculating_volumes_simple_shpaes.m -- calculates volumes and produces diagnostic plots for some simple shapes
2. calculating_disks_partition_function.m -- computes the partition function for the configuration space of n hard disks radius r in a torus.
3. calculating_disks_configuration_components_volumes.m --computes the volume of components of the configuration space for n disks radius r in a box. Also produces some diagnostic plots
4. calculating_volumes_birkhoff_polytope.m--computes the volume of the birkhoff polytope
5. contingency_tables.m-- the number of nxm contingency tables with specified row and column sums can be used to approximate the volume of transport polytopes. Examples
    2,3,4 show that this is surprisingly accurate! we can turn this idea on it's head: we can use our volume computation algorithms to approximate the number of contingency tables with given row
    and column sums (currently not implemented)




The code includes a redistribution of the files 'getCororSet.m'  downloaded from the MATLAB file exchange. The citation for this file 

Diana (2021). Color blind friendly colormap (https://www.mathworks.com/matlabcentral/fileexchange/46802-color-blind-friendly-colormap), MATLAB Central File Exchange. Retrieved November 26, 2021. 



The file 'getCororSet.m' comes with the following copyright notice:
---------------------------------------------------------------------------------------------------------------------------------
Copyright (c) 2017, Massimo Ciacci
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

