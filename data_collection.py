import numpy as np
import datetime
import os
from netCDF4 import Dataset

class DataCollection:
    def __init__(self, data):
        self._data = None
    def collect(self):
        return

    def readFile(self):
        raise NotImplementedError("not implemeted")

    def writeFile(self):
        raise NotImplementedError("not implemeted")

class AirDynamicData(DataCollection):
    def __init__(self, data):
        super().__init__(data)
        
    def air_collect(self, name, time):
        timezone = self._data["time"][:]
        for i in range(0, len(timezone)):
            if timezone[i] == time:
                return 
        return 

    def collect(self, name, time):
        return self.air_collect(name, time)

    def readFile(self, name):
        self._data = Dataset("%s.nc"%name, "r")

    def writeFile(self, name):
        if "%s.nc"%name not in os.listdir(os.getcwd()):
            f = Dataset("%s.nc"%name, "w")
        elif "%s.nc"%name in os.listdir(os.getcwd()):
            f = Dataset("%s.nc"%name, "a")
    



def Input(name, t):
	time =  Dataset("%s.nc"%name, "r")["time"][:]
	for i in range(0,len(time)):
		if t == time[i]:
			Var = Dataset("%s.nc"%name, "r+")[name]
	return Var

def Output(name, var, t):
    global nX, nY
    if "%s.nc"%name not in os.listdir(os.getcwd()):
        f = Dataset("%s.nc"%name, "w")
    elif "%s.nc"%name in os.listdir(os.getcwd()):
        f = Dataset("%s.nc"%name, "a")
    return f	

   