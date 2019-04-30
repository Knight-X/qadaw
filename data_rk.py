import numpy as np

class AbstractData:

	def __init__(self, data):
		self._data = data

	def getData(self):
		return self._data
	
	def setData(self, data):
		self._data = data


class NormalData(AbstractData):

	def __init__(self, data, boundary_data,
	level, wd):
		self._boundary_data = boundary_data
		self._level = level
		self._wd = wd


	def Spatial(self, axis):#axis option (x,y)

		var_forward = np.concatenate((self._data.take(np.arange(1, self._data.shape[axis]), axis), self._boundary_data), axis)/2
		var_backward = np.concatenate((self._boundary_data, self._data.take(np.arange(0, self._data.shape[axis]-1), axis)), axis)/2


		Var = (var_forward - var_backward)/self._wd[axis]

		return NormalData(Var, self._boundary_data, self._level, self._wd)

	def Interp(self, axis):#x is supposed to be a tuple


		wd_forward = np.concatenate((self._wd[axis], self._boundary_data), axis)
		wd_backward = np.concatenate((self._boundary_data, self._wd[axis]), axis)
		
		var_forward = np.concatenate((self._data, self._boundary_data), axis)
		var_backward = np.concatenate((self._boundary_data, self._data), axis)
		
		Var = 0.5*(wd_backward/(0.5*(wd_backward + wd_forward))*var_forward + wd_forward/(0.5*(wd_backward + wd_forward))*var_backward)
		return NormalData(Var, self._boundary_data, self._level, self._wd)

class ConstData(AbstractData):
	def __init__(self, data, level):
		self._data = data
		self._level = level
	
	def spatial(self):
		var_forward = self._data.take(np.arange(1, self._data.shape[axis]), axis)
		var_backward = self._data.take(np.arange(0, self._data.shape[axis]-1), axis)
		Var = (var_forward - var_backward)/wd[axis]
		return Var

	def interp(self):
		var_forward = self._data.take(np.arange(1, self._data.shape[axis]), axis)
		var_backward = self._data.take(np.arange(0, self._data.shape[axis]-1), axis)
		Var = 0.5*(var_forward + var_backward)	
		return Var