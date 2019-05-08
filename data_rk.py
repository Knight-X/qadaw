import numpy as np

class AbstractData:

	def __init__(self, data):
		self._data = data

	def getData(self):
		return self._data
	

class NormalData(AbstractData):

	def __init__(self, data, boundary_data,
	level, wd):
		super().__init__(data)
		self._boundary_data = boundary_data
		self._level = level
		self._wd = wd
	def level(self):
		return self._level 

	def spatial(self, axis):#axis option (x,y)

		var_forward = np.concatenate((self._data.take(np.arange(1, self._data.shape[axis]), axis), self._boundary_data), axis)
		var_backward = np.concatenate((self._boundary_data, self._data.take(np.arange(0, self._data.shape[axis]-1), axis)), axis)

		Var = (var_forward - var_backward)/self._wd[axis]

		return NormalData(Var, self._boundary_data, self._level, self._wd)

	def interp(self, axis):#x is supposed to be a tuple


		wd_forward = self._wd[axis]
		wd_backward = self._wd[axis]
		
		#what is bc_forward and bc_backward
		var_forward = np.concatenate((self._data, self._boundary_data), axis)
		var_backward = np.concatenate((self._boundary_data, self._data), axis)
		
		Var = 0.5*(wd_backward/(0.5*(wd_backward + wd_forward))*var_forward + wd_forward/(0.5*(wd_backward + wd_forward))*var_backward)
		return NormalData(Var, self._boundary_data, self._level, self._wd)

class ConstData(AbstractData):
	def __init__(self, data, level, wd):
		super().__init__(data)
		self._level = level
		self._wd = wd
	
	def spatial(self, axis):
		var_forward = self._data.take(np.arange(1, self._data.shape[axis]), axis)
		var_backward = self._data.take(np.arange(0, self._data.shape[axis]-1), axis)
		Var = (var_forward - var_backward)/self._wd[axis]
		return ConstData(Var, self._level, self._wd)

	def interp(self, axis):
		var_forward = self._data.take(np.arange(1, self._data.shape[axis]), axis)
		var_backward = self._data.take(np.arange(0, self._data.shape[axis]-1), axis)
		Var = 0.5*(var_forward + var_backward)	
		return ConstData(Var, self._level, self._wd)

