import numpy as np

class AbstractData:

	def __init__(self, data, nX, nY, nK):
		self._data = data
		self._nX = nX
		self._nY = nY
		self._nK = nK

	@property
	def getData(self):
		return self._data
	
	def level_determine(self, data):
		if len(data) == self._nK and len(data[0]) == self._nY and len(data[0][0]) == self._nX:
			lvl = "mp"
		elif len(data) == self._nK and len(data[0]) == self._nY and len(data[0][0]) == self._nX+1:
			lvl = "u"
		elif len(data) == self._nK and len(data[0]) == self._nY+1 and len(data[0][0]) == self._nX:
			lvl = "v"
		elif len(data) == self._nK+1 and len(data[0]) == self._nY and len(data[0][0]) == self._nX:
			lvl = "w"
	
		return lvl
	

class NormalData(AbstractData):

	def __init__(self, data, boundary_data,
	nX, nY, nK, wd):
		super().__init__(data, nX, nY, nK)
		self._boundary_data = boundary_data
		level = self.level_determine(data)
		self._level = level
		self._wd = wd
		
	@property
	def level(self):
		return self._level 

	def spatial(self, axis):#axis option (x,y)

		var_forward = np.concatenate((self._data.take(np.arange(1, self._data.shape[axis]), axis), self._boundary_data), axis)
		var_backward = np.concatenate((self._boundary_data, self._data.take(np.arange(0, self._data.shape[axis]-1), axis)), axis)

		Var = (var_forward - var_backward)/self._wd[axis]

		return NormalData(Var, self._boundary_data, self._nX, self._nY, self._nK, self._wd)

	def interp(self, axis):#x is supposed to be a tuple


		wd_forward = self._wd[axis]
		wd_backward = self._wd[axis]
		
		#what is bc_forward and bc_backward
		var_forward = np.concatenate((self._data, self._boundary_data), axis)
		var_backward = np.concatenate((self._boundary_data, self._data), axis)
		
		Var = 0.5*(wd_backward/(0.5*(wd_backward + wd_forward))*var_forward + wd_forward/(0.5*(wd_backward + wd_forward))*var_backward)
		return NormalData(Var, self._boundary_data, self._nX, self._nY, self._nK, self._wd)

class ConstData(AbstractData):
	def __init__(self, data, nX, nY, nK, wd):
		super().__init__(data, nX, nY, nK)
		level = self.level_determine(data)
		self._level = level
		self._wd = wd
	
	def spatial(self, axis):
		var_forward = self._data.take(np.arange(1, self._data.shape[axis]), axis)
		var_backward = self._data.take(np.arange(0, self._data.shape[axis]-1), axis)
		Var = (var_forward - var_backward)/self._wd[axis]
		return ConstData(Var, self._nX, self._nY, self._nK, self._wd)

	def interp(self, axis):
		var_forward = self._data.take(np.arange(1, self._data.shape[axis]), axis)
		var_backward = self._data.take(np.arange(0, self._data.shape[axis]-1), axis)
		Var = 0.5*(var_forward + var_backward)	
		return ConstData(Var, self._nX, self._nY, self._nK, self._wd)

