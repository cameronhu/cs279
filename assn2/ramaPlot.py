import sys,os.path

if len(sys.argv) < 2:
  print("Usage: python ramaPlot.py <filename>.pdb")
  sys.exit(1)

pdbs = []

for f in sys.argv[1:]:
  if not os.path.isfile(f):
    print("%s is not a valid path to a pdb file"%f)
  else:
    pdbs.append(f)

if not pdbs:
  sys.exit(1)

from pyrosetta.rosetta import *
from pyrosetta.toolbox import *
from pyrosetta import init
from pyrosetta.io import pose_from_pdb
import matplotlib.pylab as pl
from numpy import mgrid,zeros
import matplotlib.colors
from math import ceil

init()


for i,pdbFile in enumerate(pdbs):
  p = pose_from_pdb(pdbFile)

  nResidues = p.total_residue()

  nDivs = 2*int(ceil(nResidues**(0.67)))
  bins = zeros((nDivs,nDivs))
  psis,phis = mgrid[slice(-180,181,360./nDivs),slice(-180,181,360./nDivs)]

  for res in range(1,p.total_residue()+1):
    phiBin = int(nDivs*(p.phi(res) + 180)/360)%nDivs
    psiBin = int(nDivs*(p.psi(res) + 180)/360)%nDivs
    bins[psiBin,phiBin] += 1
    # print (p.phi(res),p.psi(res))

  pl.figure(i+1)
  pl.pcolor(phis,psis,bins,cmap='YlOrBr')
  pl.colorbar()
  pl.clim(0,bins.max()/2)
  pl.axis([phis.min(),phis.max(),psis.min(),psis.max()])
  pl.xlabel("phi")
  pl.ylabel("psi")
  pl.title(pdbFile)

pl.show()

