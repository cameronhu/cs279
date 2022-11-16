# You are to implement the stochastic and Laplacian diffusers in the code below.
# Specifically, you will implement the update(self) functions
# in the StochasticDiffuser and LaplacianDiffuser classes.
# Once you implement these update functions, you can run the script from the terminal:
#
# To run StochasticDiffuser: python diffuser.py stoc
# To run LaplacianDiffuser : python diffuser.py lap
#
# Note: use 'stoc' for the stochastic and 'lap' for the Laplacian diffuser.
# Change the initial settings and parameters for the simulation below,
# in the main method. The units in the StochasticDiffuser are the number of particles.
# The units in the LaplacianDiffuser are arbitrary measures of concentration.

import numpy as np
from math import floor, ceil

from random import random # imports random generator from 0-1

import matplotlib.pyplot as plt
from matplotlib import animation, colors

import sys

class Diffuser(object):

  def __init__(self,
               nX = 25,
               nY = 25,
               diffusion = 1,
               ):
    self.nX = nX
    self.nY = nY

    # creates a 2D array of zeros to represent the particles (initialized with zero particles) at each x,y coordinate in the cell
    self.blocks = np.zeros((nY,nX))
    self.diffusion = diffusion

  def setBlock(self,block,val):
    # NOTE: block should be a tuple containing (y-coordinate, x-coordinate)
    if block[0] >= 0 and block[0] < self.nY and block[1] >= 0 and block[1] < self.nX:
      self.blocks[block[0], block[1]] = val
    else:
      raise ValueError("X, Y coordinates out of bounds.")

  def getBlock(self,block):
    # NOTE: block should be a tuple containing (y-coordinate, x-coordinate)
    if block[0] >= 0 and block[0] < self.nY and block[1] >= 0 and block[1] < self.nX:
        return self.blocks[block[0], block[1]]
    return None

  def getGrid(self):
    return self.blocks

  def simulateTimeStep(self,nIters=1):
    for step in range(nIters):
      self.blocks = self.update()

  def update(self):
    raise NotImplementedError("Only instantiate Stochastic or Laplacian Diffusers")


  # Returns the correct x index after moving step blocks
  # to the right (to the left for negative step). Accounts for
  # the "wrap-around" at the border.
  # Note: When calling nextX() in your update() function,
  #   you should use self.nextX(x, step)
  def nextX(self,x,step):
    return int((x+step)%self.nX)
  # Returns the correct y index after moving step blocks
  # up (down for negative step). Accounts for
  # the "wrap-around" at the border.
  # Note: When calling nextY() in your update() function,
  #   you should use self.nextY(y, step)
  def nextY(self,y,step):
    return int((y+step)%self.nY)


class StochasticDiffuser(Diffuser):

  def update(self):
    # First, we will initialize a new blocks array with the same dimensions.
    # All values in new blocks are set to 0, which you will update in your code below.
    newBlocks = np.zeros((self.nY, self.nX))

    # Next, we will iterate through the coordinates of our cell,
    # so that we can then iterate through each particle at that coordinate.

    for row in range(self.nY):
      for col in range(self.nX):
        numberOfPart = self.blocks[row, col]
        # The delta/step size depend on the diffusion constant which is set in
        # the main method
        delta = ceil(self.diffusion)

        # Now we will iterate through every particle in each block
        for part in range(int(numberOfPart)):
          # Introduce some randomness in step size
          stepSize = delta - (random() > 0.7)

          # TODO:  YOUR CODE HERE
          # Fill in this section.
          # The particle should choose left, up, right, or down with equal probability.
          # Hint: use the random() function, which will return a random float number between 0 and 1.
          #
          # After you've chosen a random direction,
          # "add" a particle to the newBlock coordinate that is "stepSize" away in the chosen direction.
          #
          # Note: To achieve the "wrap-around" effect it will be helpful to use
          # the self.nextX() and self.nextY() methods.
          up = (1, 0)
          down = (-1, 0)
          left = (0, -1)
          right = (0, 1)
          directions = [up, down, left, right]
          step_direct = directions



          # END YOUR CODE HERE

    return newBlocks


class LaplacianDiffuser(Diffuser):

  def update(self):
    # First, we will initialize the derivative matrices, with all values set to 0
    ddx = np.zeros((self.nY, self.nX))
    ddy = np.zeros((self.nY, self.nX))
    ddt = np.zeros((self.nY, self.nX))

    # Next, we will iterate through the coordinates of self.blocks
    # to compute the first discrete derivatives (change in particles with respect to x or y).
    for row in range(self.nY):
      for col in range(self.nX):

        # TODO: YOUR CODE HERE
        # Compute a discrete (first) derivative in x and y.
        # i.e., fill in ddx, ddy
        # Again, self.nextX() and self.nextY() will be helpful.

        pass # remove this line once you start implementing


        # END YOUR CODE HERE

    # Finally, we will once again iterate through the coordinates of self.blocks to compute the second discrete derivatives (ie the change in the first derivative with respect to x or y)
    for row in range(self.nY):
      for col in range(self.nX):
        # TODO:  YOUR CODE HERE
        # Compute the discrete second derivative in x and y.
        # i.e. the derivative of the derivatives you computed above
        # Make sure to balance your discrete derivatives so that the are not
        # using the same interval over x or y for both derivatives. (see assignment sheet for details)

        lap_x = None # fill in with your own code
        lap_y = None # fill in with your own code

        # END YOUR CODE HERE

        ddt[row, col] = self.diffusion*(lap_x + lap_y)

    # Update the changes in concentration
    for row in range(self.nY):
      for col in range(self.nX):
        self.blocks[row, col] += ddt[row, col]
    return self.blocks


def main():
  usage = "Usage: python diffuser.py [stoc | lap]"
  cm = "cool"

  iters = 1000 # number of iterations for diffusion
  length = 25 # size (x = y) of the plane

  if len(sys.argv) < 2:
    print(usage)
    exit(1)

  # EDIT DIFFUSION PARAMETERS below for Stochastic and Laplacian Diffusers
  if sys.argv[1] == "stoc":
    # play around with different values for diffusion
    # (it only makes sense to use integers here)
    diffusion = 1
    c = StochasticDiffuser(nX = length,
                           nY = length,
                           diffusion = diffusion,
                           )
  elif sys.argv[1] == "lap":
    # play around with different values for diffusion
    # (note what happens when diffusion > 0.25)
    diffusion = 0.24
    c = LaplacianDiffuser(nX = length,
                          nY = length,
                          diffusion = diffusion,
                          )
  else:
    print(usage)
    exit(1)


  """ Here are a variety of settings for the initial conditions of the
      diffusion simulations.  Feel free to try some of your own.
      Note: The units in StochasticDiffuser are number of particles.
            The units in LaplacianDiffuser are an aribitrary measure of concentration.
  """

  ##############################################################
  # EDIT INITIAL SIMULATION PARAMETERS BELOW.

  """SINGLE PARTICLE"""
  c.setBlock((length//2,length//2),1)

  """Point Mass"""
  # c.setBlock((length//2,length//2),1250)

  """1D Diffusion"""
  # for i in range(length):
    # c.setBlock((i,length//2),125)

  """1D Gradient Diffusion"""
  #for i in range(length):
  #  for j in range(length):
  #    c.setBlock((i,j), 2*abs(int(length/2-i)))

  # END PARAMETERS
  ##############################################################

  # This portion of the code runs the simulation through matplotlib
  for l in range(iters):
    plt.clf()
    grid = c.getGrid()
    plt.pcolor(np.array([i for i in range(length+1)]),np.array([j for j in range(length+1)]),grid,cmap = cm)
    plt.colorbar()
    # un-comment the line below to 'fix' the colorbar to the range [0,3] so it is not longer dynamic
    # 'Fixing' the colorbar will make visualization easier, but may cut off the higher values
    # plt.clim(0,3)
    plt.xlabel("Iteration %d"%l)
    plt.pause(0.1)
    c.simulateTimeStep(1)


# Boiler plate invokes the main function
if __name__ == "__main__":
  main()
