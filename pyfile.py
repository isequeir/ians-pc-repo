###############
#  testing git stuff/repository and basic commands, pycharm automatically saveswhich is nice


# Toy graphene model

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function
from pythtb import * # import TB model class
from math import sqrt, pi
import numpy as np
import matplotlib.pyplot as plt

# define lattice vectors in units of nm
a = .25
a_cc = .14
c0 = .3

# set model parameters in units of eV
delta = 0 #models potental energy from eletric field, top layer gets +v/2 bottom layer -v/2
t0 = -3.1 #intralayer hopping amplitude
t1 = -.5 #interlayer hopping amplitude

# # define lattice vectors
#
# lat=[[1.0,0.0],[0.5,np.sqrt(3.0)/2.0]]
# define lattice vectors
lat = [[a/2.0, a*sqrt(3.0)/2.0], [-a/2.0, a*sqrt(3.0)/2.0]]
# define coordinates of orbitals/ or basis wrt lattice points
lat_p = [[0, 0], [0,0]]
orb = [[0.0, 0.0], [4*a_cc/sqrt(3), 4*a_cc/sqrt(3.0)]]


# make two dimensional tight-binding graphene model
my_model=tb_model(2,2,lat,orb)
latticepoints = tb_model(2,2,lat,lat_p)
latticepoints.set_hop(t0, 0, 1, [0, 0])
latticepoints.set_hop(t0, 1, 0, [1, 0])
latticepoints.set_hop(t0, 1, 0, [0, 1])

# set on-site energies
my_model.set_onsite([-delta,delta])
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t0, 0, 1, [0, 0])
my_model.set_hop(t0, 1, 0, [1, 0])
my_model.set_hop(t0, 1, 0, [0, 1])

(fig,ax) = my_model.visualize(0,1)
line_cut=my_model.cut_piece(8,0,glue_edgs=False)
area_cut = line_cut.cut_piece(8,1, glue_edgs=False)
area_cut_lp = line_cut.cut_piece(8,1, glue_edgs=False)
(fig, ax) = area_cut.visualize(dir_first=0,dir_second=1, draw_hoppings=True)

#plt.show()




# print tight-binding model
my_model.display()

# generate list of k-points following a segmented path in the BZ
# list of nodes (high-symmetry points) that will be connected
# Gamma = [0, 0]
# K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
# M = [0, 2*pi / (3*a_cc)]
# K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]

path=[[2./3.,1./3.],[.5,.5],[0.,0.],[1./3.,2./3.]]
# labels of the nodes
label=( r'$K1$', r'$M$', r'$\Gamma $',r'$K2$' )
# total number of interpolated k-points along the path
nk=121

# call function k_path to construct the actual path
(k_vec,k_dist,k_node)=my_model.k_path(path,nk)
# inputs:
#   path, nk: see above
#   my_model: the pythtb model
# outputs:
#   k_vec: list of interpolated k-points
#   k_dist: horizontal axis position of each k-point in the list
#   k_node: horizontal axis position of each original node

print('---------------------------------------')
print('starting calculation')
print('---------------------------------------')
print('Calculating bands...')

# obtain eigenvalues to be plotted
evals=my_model.solve_all(k_vec)

# figure for bandstructure

fig, (ax, ax1) = plt.subplots(2, 1)
# specify horizontal axis details
# set range of horizontal axis
ax.set_xlim(k_node[0],k_node[-1])
# put tickmarks and labels at node positions
ax.set_xticks(k_node)
ax.set_xticklabels(label)
# add vertical lines at node positions
for n in range(len(k_node)):
  ax.axvline(x=k_node[n],linewidth=0.5, color='k')
# put title
#ax.set_title("Graphene band structure")
#ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

# plot first and second band
ax.plot(k_dist,evals[0], label="$\pi$")
ax.plot(k_dist,evals[1], label="$\pi*$")
ax.legend()
print('Done.\n')


print()
print('---------------------------------------')
print('starting DOS calculation')
print('---------------------------------------')
print('Calculating DOS...')

# calculate density of states
# first solve the model on a mesh and return all energies
kmesh=23*23 # needs to be square because than the dos will look symmetric
kpts=[]
for i in range(kmesh):
    for j in range(kmesh):
        kpts.append([float(i)/float(kmesh),float(j)/float(kmesh)])
# solve the model on this mesh
evals=my_model.solve_all(kpts)
#flatten completely the matrix
evals=evals.flatten()

# plotting DOS
print('Plotting DOS...')

# now plot density of states
bins = 101 #want bins to be an odd number
dos, energy, patches = ax1.hist(evals, bins)
# put title
#ax1.set_title("DOS")
ax1.set_xlabel("Energy [eV]")
ax1.set_ylabel("Number of states")
print('Done.\n')


energy_bin = []
length = len(energy)
for i in range(0, length-1):
    energy_bin.append([(energy[i]+energy[i+1])/2])


ax1.plot(energy_bin, dos, color="r")
plt.show()

