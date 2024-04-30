#######################################################################
# solve L[phi] - m^2 phi = rho using multigrid relaxation
#######################################################################

import numpy as np
from PIL import Image
from scipy.ndimage import convolve, maximum_filter

#######################################################################

# enable diagnostic output
verbose = True

class grid:
	# initialize grid hierarchy (with optional mask)
	def __init__(self, phi, dx, m2=0.0, rhs=0.0, mask=None, bc='constant', cval=0.0, level=0):
		self.phi = phi; self.rho = rhs; self.mask = mask
		self.level = level; self.bc = bc; self.cval = cval

		# inherit size attributes from initial conditions
		self.shape = phi.shape; self.size = phi.size
		self.pixels = phi.size - np.count_nonzero(mask)

		# log initialization progress to stdout
		if verbose: print(f'Initializing {self.shape} grid, {self.pixels} unmasked pixels')

		# second order Laplacian stencil (isotropic)
		laplacian = np.array([[1,4,1],[4,-20,4],[1,4,1]])/(6.0*dx**2)
		laplacian[1,1] -= m2; self.stencil = laplacian

		# set solver timestep to maximally stable one
		assert laplacian[1,1] < 0.0, "multigrid: mass term is unstable"
		self.dx = dx; self.dt = -1.0/laplacian[1,1]

		# initialize coarse subgrid (if large enough)
		if self.pixels > 256:
			cmask = None if (mask is None) else self.shrink(mask)
			self.coarse = grid(self.downgrade(phi), 2.0*dx, m2, mask=cmask, bc=bc, cval=cval, level=level+1)
		else:
			self.coarse = None

	# residual of stencil[phi] = rho (optionally masked)
	def residual(self):
		R = convolve(self.phi, self.stencil, mode=self.bc, cval=self.cval) - self.rho
		return R if self.mask is None else np.where(self.mask == 0.0, R, 0.0)
	
	# smooth the solution by a few diffusive steps
	def smooth(self, iterations):
		for i in range(iterations):
			self.phi += self.residual()*self.dt
	
	# shrink mask to coarser grid
	def shrink(self, mask):
		return maximum_filter(mask,5)[::2,::2]

	# downgrade data to coarser grid
	def downgrade(self, data, resample=Image.BILINEAR):
		h,w = data.shape; img = Image.fromarray(data)
		return np.array(img.resize((w>>1,h>>1), resample=resample))

	# interpolate data to finer grid
	def upgrade(self, data, resample=Image.BILINEAR):
		h,w = self.shape; img = Image.fromarray(data)
		return np.array(img.resize((w,h), resample=resample))

	# multigrid W-stroke
	def wstroke(self, fine=4, last=64):
		if self.coarse is None:
			self.smooth(last)
		else:
			self.smooth(fine)
			self.coarse.phi.fill(0.0)
			self.coarse.rho = self.downgrade(self.residual())
			self.coarse.wstroke(); self.coarse.wstroke()
			self.phi -= self.upgrade(self.coarse.phi)
			self.smooth(fine)

	# solve the problem
	def solve(self, iterations=16, feather=8):
		for i in range(iterations):
			self.wstroke(fine=4 if self.mask is None else 32)
			if verbose:
				residual = self.residual()
				norm = np.linalg.norm(residual)
				maxv = np.max(np.abs(residual))
				print(f'Iteration {i}: |residual| = {norm/self.pixels}, max = {maxv}')
		self.smooth(feather)
		return self.phi
