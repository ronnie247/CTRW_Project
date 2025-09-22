# CTRW_Project
Codes and SI for the first CTRW paper in JCP (2024)

Update - as of March 19, 2024.
The following codes have been uploaded to the repository:
1. makepore.py - makes the pore for the Nav 1.7 ion channel
2. Dyn.py - runs the CTRW on any of the lattices using the "dynamic" assignment of rates. Uses the Log-Normal Distribution by default, to use the shifted gaussian, change the use of the function as given at the top of the code.
3. Stat.py - runs the CTRW on any of the lattices using the "static" assignment of rates. Uses the Log-Normal Distribution by default, to use the shifted gaussian, change the use of the function as given at the top of the code.

They need the following input files:
1. 6j8i.pdb - pdb file for the Nav1.7, untouched.
2. ionpore3.pdb - pdb file for the Nav1.7, only the pore domain.
3. 1c4dIC1.pdb - pdb file for ONE gramicidin channel with no water/ions.
4. GramPore.xyz - pore coordinates for the gramicidin channel.
5. ILPore.xyz - pore coordinates for the Nav 1.7 ion channel.
6. PoreCoordinates1.xyz - also pore coordinates, output of the makepore.py file.
