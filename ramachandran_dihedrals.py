# Program for calculation of all dihedral (phi, psi,omega) for a given protein pdb file

#This program requires python2.7 or later with math and pylab modules installed


import math,sys
import pylab as pl

try:

	pdb_file = sys.argv[1]
	out_file = sys.argv[2]

except:

	print("Usage: all_dihedral.py {input.pdb} {out_filename}")
	raise SystemExit


# A list varibale to store atom_name, atom_number, x,y,z co-ordinates  
dihed_atoms = []

# opening and Iterating through file line by line

with open(pdb_file,'r') as fopen:
	
  for line in fopen:

# Parsing atomic co-ordinates line  
		
		if line[0:4] == 'ATOM':
			if (line[12:16]).strip() in ['N','CA','C']:
				dihed_atoms.append(((int(line[5:11])),(line[12:16]).strip(),(line[22:26]).strip(),float(line[30:38]),float(line[38:46]),float(line[46:54])))


# Function for finding vector form 3D cartesian co-ordinates
def cart_to_vec(x,y):

	vec = [(x[0] - y[0]), (x[1] - y[1]), (x[2] - y[2])]
		
	return vec

# Function for calculation of cross_product 

def cross_prod(x,y):
	
	xy = (((x[1] * y[2]) - (x[2] * y[1])), -((x[0] * y[2]) - (x[2] * y[0])), ((x[0] * y[1]) - (x[1] * y[0])))
	
	return xy

# Function for normalization of a vector
def norm_vec(x):
	
	m = (((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]))**0.5)
	
	n = [(x[0]/m), (x[1]/m), (x[2]/m)]
	
	return n

# Function for dot product calculation
def dot_prod(x,y):
	
	dp = ((x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]))

	return dp


# Function for calculation of dihedral angle
def dihed(atoms,start,stop):
	for i in range(start,len(atoms)-stop,3):

	
		A = cart_to_vec(atoms[i][3:6],atoms[i+1][3:6])

		B = cart_to_vec(atoms[i+1][3:6],atoms[i+2][3:6])

		C = cart_to_vec(atoms[i+3][3:6],atoms[i+2][3:6])
	
	
		n1 = cross_prod(A,B)
		n2 = cross_prod(C,B)
		
		
		N1 = norm_vec(n1)
		N2 = norm_vec(n2)
		
		sign = dot_prod(cross_prod(n2,n1),B)
		if sign < 0:

			dihed_angle = -(math.degrees(math.acos(dot_prod(N1,N2))))
		else:
			dihed_angle = (math.degrees(math.acos(dot_prod(N1,N2)))) 	
			
		yield "dihed between {} {} {} {} of residues {} and {}  is {}".format (atoms[i][0],atoms[i+1][0],atoms[i+2][0],atoms[i+3][0],atoms[i][2],atoms[i+3][2],dihed_angle)
	
fout = open(out_file,'w')

fout.write("# Psi dihedral angles between Ni - CAi - Ci - Ni+1 #\n\n")

# calculation and writing psi dihedral angle to file
x = []
for psi in dihed(dihed_atoms,0,3):
	
	s = psi.split()
	x.append(s[12])
	fout.write("Psi {}\n".format (psi))


# calculating and writing omega dihedral angle to file

fout.write("\n\n# Omega dihedral angles between CAi - Ci - Ni+1 - CAi+1 #\n\n")

for omega in dihed(dihed_atoms,1,2):
	
	fout.write("Omega {}\n".format(omega))
	
# Calculating and writing phi dihedral atoms to file

fout.write("\n\n# Phi dihedral angles between Ci - Ni+1 - CAi+1 - Ci+1 #\n\n")

y = [] 

for phi in dihed(dihed_atoms,2,1):
	
	s = phi.split() 
	y.append(s[12])
	fout.write("Phi {}\n".format(phi))

fout.close()
# Ploting ramachndran plot
pl.plot(y,x,'ro')

# make axis labels
pl.xlabel('phi angle',fontsize=40)
pl.ylabel('psi angle',fontsize=40)
pl.xticks(fontsize=20)
pl.yticks(fontsize=20)
# set axis limits
pl.xlim(-180, 180)
pl.ylim(-180, 180)

pl.show()
