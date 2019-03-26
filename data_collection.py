import numpy as np
import datetime
import os
from netCDF4 import Dataset

class DataCollection:
    def __init__(self, data=None):
        self._data = data

    def collect(self):
        raise NotImplementedError("not implemeted")

    def readFile(self):
        raise NotImplementedError("not implemeted")

    def writeFile(self):
        raise NotImplementedError("not implemeted")

class AirDynamicData(DataCollection):
    def __init__(self, name):
        data = Dataset("%s.nc"%name, "r")
        super().__init__(data)
        self._time_map = self.prepare(data)
        
    def air_collect(self, name, time):
        return self._data[name][self._time_map[time]].filled()

    def collect(self, name, time):
        return self.air_collect(name, time)

    def readFile(self, name):
        self._data = Dataset("%s.nc"%name, "r")

    def writeFile(self, name):
        if "%s.nc"%name not in os.listdir(os.getcwd()):
            f = Dataset("%s.nc"%name, "w")
        elif "%s.nc"%name in os.listdir(os.getcwd()):
            f = Dataset("%s.nc"%name, "a")
    
    def prepare(self, data):
        timezone = data["time"][:]
        time_map = {}
        for i in range(0, len(timezone)):
            time_map[timezone[i]] = i
        return time_map

    
   