from netCDF4 import Dataset

class DataCollection:
    def __init__(self, data):
        self._data = data
    def collect(self):
        return

class AirDynamicData(DataCollection):
    def __init__(self, data):
        super().__init__(data)
        
    def air_collect(self):
        return 
    def collect(self):
        return self.air_collect()

def Input(name, t):
	time =  Dataset("%s.nc"%var, "r")["time"][:]
	for i in range(0,len(time)):
		if t == time[i]:
			Var = Dataset("%s.nc"%var, "r+")[var]
	return Var


def Output(name, var, t):
	global nX, nY
	if "%s.nc"%name not in os.listdir(os.getcwd()):
		f = Dataset("%s.nc"%name, "w")
		
	elif "%s.nc"%name in os.listdir(os.getcwd()):
		f = Dataset("%s.nc"%name, "a")