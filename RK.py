import numpy as np
from data_collection import AirDynamicData
import datetime
from data_rk import NormalData
#Initialization
data_collector 		= AirDynamicData("climate2014-01.nc")
BC_collecctor 		=  AirDynamicData("wrfbdy_d01")
#user define
t_start 			= datetime.datetime(2016, 8, 14, 0, 0, 0) 
t_end 				= datetime.datetime(2016, 8, 14, 6, 0, 0)
Step 				= datetime.timedelta(seconds = 3600*6) # coorelated to interval seconds in wrf namelist
wdX					= 11 #kilometer
wdY 				= 11 #kilometer
nX 					= 200
nY   				= 300
nK					= 31 #equals to the numbers of vertical isobaric layers from global model
num_bdy_layer 		= 1

wdK 				= None
wd 					= None
idxK 				= None,
idxW 				= None,
lat 				= None,
lon 				= None,

t_now 				= t_start
RK3_start 			= None
RK3_end				= None
Step_RK3 			= datetime.timedelta(seconds = 6) # 6*wdX sec is suggested
acoustic_start 		= None
acoustic_end 		= None
Step_acoustic		= None
n 					= None
ns 					= 8

p_sd		= None, 
p_dry_top   = None,
p0 			= None#pa
T0 			= None #K, mean sea level temp
A 			= None #K, difference between p0 and p0/e

g 			= 9.81
Rd 			= 287 #J/kg/K
Rv 			= 461.6 #J/kg/K
rEarth 		= 6.307*10^6 #meter
Beta 		= 0.1
e 			= None 
f 			= None 
sinalpha 	= None 
cosalpha 	= None

#input from wrfd01
mapFactorM 	= None,
mapFactorU 	= None,
mapFactorV 	= None,

p_refer 	= None, 
gp_refer 	= None, 
Mdry_refer 	= None, 
Ddry_refer 	= None,

p_pfr   	= None, 
gp_pfr   	= None, 
Mdry_pfr   	= None, 
Ddry_pfr   	= None,

u  			= None, 
v  			= None, 
w    		= None, 
o 			= None, 
tp  		= None,
p  			= None,
gp 			= None, 
Mdry 		= None, 
Ddry 		= None,
D  			= None, 
qv 			= None, 
qc 			= None, 
qr   		= None, 
qi 			= None, 
qm   		= None,

#coupled
U 			= None, 
V			= None, 
W			= None, 
O 			= None, 
Tp 			= None, 
Qv 			= None, 
Qc 			= None, 
Qr 			= None, 
Qi 			= None, 
Qm 			= None,


#reference state
def cp():
    global Rd
	c_p= 7*Rd/2
	return c_p

def cv():
	global Rd
	c_v = cp()−Rd
	return c_v

def gamma():
    gma = cp()/cv()
	return gma

def cs_sqr():
	global p,Ddry
	cs2 = gamma()*p*Ddry
	return cs2

def p_dry_surface():
	global A, p0, T0, Rd, gp_surface
	pdhs = p0*np.exp(-T0/A+((T0/A)^2-2*gp_surface/A/Rd)^0.5)
	return pdhs

def p_reference():
	global K
	pdr = K*Mdry_reference()+p_dry_top
	return pdr

def T_reference():
	global A, p0, T0
	Tr = T0+A*np.log(p_reference()/p0)
	return Tr

def tp_reference():
	global T0, A, p0, Rd
	tpdr = (T0+A*np.log(p_reference()/p0))*(p0/p_reference())^(Rd/cp())
	return tpdr

def Ddry_reference():
	global p0, Rd
	Ddr = Rd*tp_reference()/p0*(p_reference()/p0)^(-cv()/cp())
	return Ddr

def Mdry_reference():
	Mdr = p_dry_surface()-p_dry_top
	return Mdr

#perturbation state from reference
def Dd():
	global Rd, Rv, p0, qv, tp, p_pfr, p_refer
	dry_density = Rd/p0*tp*(1+Rv/Rd*qv)*(p_pfr + p_refer/p0)^(-cv()/cp())
	return dry_density

def Ddry_perturb_from_refer():
	Dd_pfr = Dd() - Ddry_reference()
	return Dd_pfr

#diagnostic eq
def tp_moist():
	global tp, qv
	tp_m = tp*(1+(Rv/Rd)*qv)
	return tp_m

def p_moist():
	global p0, Ddry
	pressure = p0*(Rd*tp_moist()/p0/Ddry)^gamma()
	return pressure 

#level identification
def level_determine(var):
	global nX, nY, nK

	if len(var) == nK and len(var[0]) == nY and len(var[0][0]) == nX:
			lvl = "mp"
	elif len(var) == nK and len(var[0]) == nY and len(var[0][0]) == nX+1:
			lvl = "u"
	elif len(var) == nK and len(var[0]) == nY+1 and len(var[0][0]) == nX:
			lvl = "v"
	elif len(var) == nK+1 and len(var[0]) == nY and len(var[0][0]) == nX:
			lvl = "w"
	
	return lvl


#import variable and turn it to coupled form
def coupled(var):# a is any variable
	global Mdry, mapFactorM, mapFactorU, mapFactorV

	if var.level() == "mp":
		M = Mdry, mF = 1
	elif var.level() == "u":
		M = Mdry.Inter(2), mf = mapFactorU
	elif var.level() == "v":
		M = Mdry.Inter(1), mf = mapFactorV
	elif var.level() == "w":
		M = Mdry, mf = mapFactorM
	Var = M.getData()*var.getData()/mF

	return Var

def uncoupled(Var):
	global Mdry
	if var.level() == "mp":
		M = Mdry, mF = 1
	elif var.level() == "u":
		M = Mdry.Inter(2), mf = mapFactorU
	elif var.level() == "v":
		M = Mdry.Inter(1), mf = mapFactorV
	elif var.level() == "w":
		M = Mdry, mf = mapFactorM
	var = *mf*Var.getData()/M.getData()

	return var

#Operator
def timeAvg(name, var_advanced, start, Step):
	global Beta, cplst

	var_prev = data_collector.collect(name, start - step)
	var_timeAvg = (1+Beta)/2*var_advanced+(1-Beta)/2*var_prev

	return var_timeAvg




def Grad_Interp(var, axis = 0)
	global wd
	
	
	wd_forward  = wd[axis]
	wd_backward = wd[axis]

	wd = 0.5*(wd_backward + wd_forward)
	Var = var*wd

	for i in range(Var.shape[axis], -1, -1):
		if i == 0 :
			Var[i] = Var[i] + upper_bc
		else:
			Var[i] = Var[i] + Var[i+1]

	if Var[0] != lower_bc:
		Var[0] = lower_bc
	
	return Var

def Vert_Integrat(var):
	global wd

	product = 0
	axis  = 0
	for i in range(0, var.shape[axis]):
		product += var[i]
	product = np.expand_dims(product, 0)
	
	return product

#tend Term
def Fcor(var):#axis is supposed to be "U", "V" or "W"
	global t_now, rEarth, e, f, sinalpha, cosalpha, mapFactorM, rotationAngle, U, V, W

	if var == "U":
		Vxy = V.interp(1).interp(2)
		Wxk = W.interp(0).interp(2)
		uncopuled_u = uncopule(U)
		uncopuled_v = uncopule(V)
		internal_data = uncopuled_u.interp(2).getData() * mapFactorM.spatial(1).getData() -
			uncopuled_v.interp(1).getData() * mapFactorM.spatial(2)
		internal = NormalData(internal_data)
		tend1 = Vxy*(
			f.interp(2).getData() + internal.interp(2).getData()
			)
		tend2 = Interp(e,2)*Wxk*Interp(cosplpha, 2)
		tend3 = uncouple(U)*Wxk/rEarth
        tend = tend1 - tend2 - tend3
	elif var == "V":
		Uxy = Interp(Interp(U,2),1)
		Uxy = U.interp(2).interp(1)
		Wyk = W.interp(0).interp(1)
		uncopuled_u = uncopule(U)
		uncopuled_v = uncopule(V)
		internal_data = uncopuled_u.interp(2).getData() * mapFactor.spatail(1).getData()
			- uncoupled_v.interp(1).getData() * mapFactorM.spatial(2).getData()
		internal = NormalData(internal_data)
		tend = - Uxy.getData() *(
			f.interp(1).getData() + internal.interp(1).getData()
			)
			+ e.interp(1).getData() *Wyk.getData() * sinalpha.interp(1).getData()
			- uncopuled_v.getData() *Wyk.getData() /rEarth
		internal = 
	elif var == "W":
		Uxk = U.interp(2).interp(0)
		Vyk = V.interp(1).interp(0)
		uncopuled_v = uncopule(V)
		uncopuled_u = uncopule(U)
		tend = e.getData() *(
			Uxk.getData() *cosplpha.getData()
			- Vyk.getData() *sinalpha.getData()
			)
			+ (
				uncouple_u.interp(2).interp(0).getData() * Uxk.getData()
				+ uncouple_v.interp(1).interp(0).getData() * Vyk.getData()
				)/rEarth

	else:
		tend = 0

	return tend

def Adv(var):
	global mapFactorM, mapFactorM, U, V, O, gp, acoustic_start, Step_acoustic, ns
	if var == "U":
		U_spatial_2 = U.spatial(2)
		U_spatial_1 = U.spatial(1)
		U_spatial_0 = U.spatial(1)
		V_Spatial = V.spatial(1)
		O_Spatial = O.spatial(0)
		tend = - mapFactorM*(U_spatial_2.getData() * uncoupled(U_spatial_2) + V_spatial.getData() *uncoupled(U_spatial_1)) - O_spatial.getData() *uncoupled(U_spatial_0)
	elif var == "V":
		V_spatial_2 = V.spatial(2)
		V_spatial_1 = V.spatial(1)
		V_spatial_0 = V.spatial(1)
		U_spatial = U.spatial(2)
		O_spatial = O.spatial(0)
		tend = - mapFactorM*(U_spatial.getData() * uncoupled(V_spatial_2) + V_spatial_1.getData() *uncoupled(V_spatial_1)) - mapFactorX/mapFactorY* O_spatial.getData() * uncoupled(V_spatial_0)
	elif var == "Mdry":
		U_spatial_2 = U.spatial(2)
		V_spatial_1 = V.spatial(1)
		O_spatial = O.spatial(0)
		tend = - mapFactorM*mapFactorM*(U_spatial_2.getData() + V_spatial_1.getData()) - mapFactorY* O_spatial.getData()
	elif var == "Tp":
		Tp_spatial_2 = Tp.spatial(2)
		Tp_spatial_1 = Tp.spatial(1)
		Tp_spatial_0 = Tp.spatial(1)
		Tp_coupled_2 = uncopuled(Tp_spatial_2)
		Tp_coupled_1 = uncopuled(Tp_spatial_1)
		Tp_coupled_0 = uncopuled(Tp_spatial_0)
		Tp_interp_2 = Tp_coupled_2.interp()
		Tp_interp_1 = Tp_coupled_1.interp()
		Tp_interp_0 = Tp_coupled_0.interp()
		U_spatial = U.spatial(2)
		V_spatial = V.spatial(1)
		O_spatial = O.spatial(0)
		tend = - mapFactorM*mapFactorM*(U_spatial.getData() *Tp_interp_2.getData() + V_spatial.getData()* Tp_interp_1.getData()) - mapFactorY*O_spatial.getData()*Tp_interp_0.getData()
	elif var == "W":
		U_spatial = U.spatial(2)
		V_spatial = V.spatial(1)
		O_spatial = O.spatial(0)
		W_spatial_2 = W.spatial(2)
		W_spatial_1 = W.spatial(1)
		W_spatial_0 = W.spatial(0)
		tend = - mapFactorM*mapFactorM/mapFactorM*(U_spatial.getData() * uncoupled(W_spatial_2) + V_spatial.getData()*uncoupled(W_spatial_1)) - O_spatial.getData() *uncoupled(W_spatial_0)
	elif var == "gp":
		Gp_spatial_2 = gp.spatial(2)
		gp_spatial_1 = gp.spatial(1)
		gp_spatial_0 = gp.spatial(0)
		U_spatial = U.spatial(2)
		V_spatial = V.spatial(1)
		O_spatial = O.spatial(0)
		tend = - (mapFactorM*mapFactorM*(U.getData()*gp_spatial_2.getData() + V.getData() * gp_spatial_1.getData() + mapFactorY*O_spatial.getData() * gp_spatial_0.getData()/Mdry
	else:
		tend = 0

	return tend

#Time integration
def ModeS(var, tend):
	global mapFactorM, U, V, O, gp, D, Ddry, Mdry, g, W, qc, qr, qv, p, p_refer, gp_refer, Mdry_refer, Ddry_refer, p_pfr, gp_pfr, Mdry_pfr, Ddry_pfr 

	if var == "U":
		ratio = D.getData() / Ddry.getData()
		interal = NormalData(ratio)
		first_part = Mdry.interp(2).getData() * D.interp(2).getData() * p_pfr.interp(2).getData()
			- Mdry.interp(2).getData() * Ddry_pfr.interp(2).getData() * p_refer.interp(2).getData()
		
		S = -(first_part) 
				- interal.interp(2).getData()
				*(Mdry.interp(2).getData() * gp_pfr.interp(0).spatial(2).getData()
					-p_pfr.interp(0).interp(2).spatail(0).getData() * gp.interp(0).spatial(2).getData() 
					+Mdry_pfr.interp(2).getData() * gp.interp(0).spatial(2).getData()
			+ tend

	elif var == "V":
		ratio = D.getData() / Ddry.getData()
		interal = NormalData(ratio)
		first_part = Mdry.interp(1).getData() * D.interp(1).getData() * p_pfr.spatial(1).getData()
			- Mdry.interp(1).getData() * Ddry_pfr.interp(1).getData() * p_refer.spatial(1).getData()
		S = -(first_part) 
				-internal.interp(1).getData()
				*(Mdry.interp(1).getData() * gp_pfr.interp(0).spatial(1).getData()
					-p_pfr.interp(0).interp(1).spatial(0).getData() * gp.interp(0).spatial(1).getData()
					+Mdry_pfr.interp(1).getData() * gp.interp(0).spatial(1).getData()
			+ tend

	elif var == "Mdry":
		S =  tend

	elif var == "Tp":
		S =  tend

	elif var == "W":
		ratio = D.getData() / Ddry.getData()
		interal = NormalData(ratio)
		S = g*(internal.interp(0).getData() *( p_pfr.spatial(0).getData() + 
		Mdry_refer.getData() * (qr.getData() +qc.getData() +qv.getData()))-Mdry_pfr.getData() )/mapFactorM.getData()
			+ tend

	elif var == "gp":
		S = mapFactorM.getData() *g* W.getData() /Mdry.getData()
			+ tend
	
	return S


def acoustic_time_integration(var_S, n, RK3_start, Step_acoustic):
	global acoustic_start, acoustic_end, Step_acoustic, RK3_start, mapFactorM, U, V, O, gp, D, Ddry, Mdry, g, W, qc, qr, qv, p, p_refer, gp_refer, Mdry_refer, Ddry_refer, p_pfr, gp_pfr, Mdry_pfr, Ddry_pfr 

	U_S = var_S[0], V_S = var_S[1], W_S = var_S[2],
	Mdry_S = var_S[3], Tp_S = var_S[4], gp_S = var_S[5]

	U_perturb 	 = coupled(data_collector.collect("U", RK3_start )) - U
	V_perturb 	 = coupled(data_collector.collect("V", RK3_start )) - V
	W_perturb 	 = coupled(data_collector.collect("W", RK3_start )) - W
	O_perturb 	 = coupled(data_collector.collect("O", RK3_start )) - O
	Tp_perturb 	 = coupled(data_collector.collect("T", RK3_start )) - Tp
	gp_perturb 	 = data_collector.collect("PH", RK3_start ) 		- gp_pfr
	Mdry_perturb = data_collector.collect("MU", RK3_start ) 		- Mdry_pfr
	Ddry_perturb = -(Spatial(gp_perturb, 0) + Ddry*Mdry_perturb)/Mdry 
	p_perturb_data 	 = cs_sqr()/Ddry*(Tp_perturb/Tp - Ddry_perturb/Ddry - Mdry_perturb/Mdry)

	Output("U_perturb", U_perturb, RK3_start )
	Output("V_perturb", V_perturb, RK3_start )
	Output("W_perturb", W_perturb, RK3_start )
	Output("O_perturb", O_perturb, RK3_start )
	Output("Mdry_perturb", Mdry_perturb, RK3_start )
	Output("Tp_perturb", Tp_perturb, RK3_start )
	Output("gp_perturb", gp_perturb, RK3_start )
	Output("Ddry_perturb", Ddry_perturb, RK3_start )
	Output("p_perturb", p_perturb_data, RK3_start)


	for j in range(0, n):
		acoustic_start = RK3_start + j*Step_acoustic
		p_perturb = NormalData(data_collector.collect("p_perturb", acoustic_start))
		ddry_perturb = NormalData(data_collector.collect("Ddry_perturb", acoustic_start))
		gp_perturb = NormalData(data_collector.collect("gp_perturb", acoustic_start))
		p_perturb = NormalData(data_collector.collect("p_perturb", acoustic_start))
		mdry_perturb = NormalData(data_collector.collect("Mdry_perturb", acoustic_start))
		U_perturb_data = U_perturb + Step_acoustic*
		( U_S - 
			((mapFactorM/mapFactorM)*(D/Ddry)*(
				Mdry*(Ddry*p_perturb.spatial(2).getData()+Spatial(p_refer,2)*ddry_perturb.getData()+gp_perturb.spatial(2).getData())
				+Spatial(gp, 2)*(p_perturb.spatial(0).getData()-mdry_perturb.getData()))
			)
		)
		U_perturb = NormalData(U_perturb_data)
		V_perturb_data = V_perturb + Step_acoustic*
		(V_S - 
			((mapFactorM/mapFactorM)*(D/Ddry)*(
				Mdry.getData() *(Ddry.getData() *p_perturb.spatial(1).getData() + p_refer.spatial(1).getData() * ddry_perturb.getData() + gp_perturb.spatial(1).getData())
				+gp.spatial(1).getData() * (p_perturb.spatial(0).getData() - mdry_perturb.getData()))
			)
		)
		V_perturb = NormalData(V_perturb_data)
		Mdry_perturb_difference = mapFactorM*mapFactorM*Vert_Integrat(Spatial(U_perturb, 2) + Spatial(V_perturb, 1))
		
		Mdry_perturb_data = Mdry_perturb.getData() + Step_acoustic*Mdry_perturb_difference 		  
		O_perturb_data = Grad_Interp(Mdry_S - Mdry_perturb_difference - mapFactorM*mapFactorM*(U_perturb.spatial(2).getData()  + V_perturb.spatial(1).getData()))/mapFactorM
		U_perturb_tp = NormalData(U_perturb.getData() * tp.getData())
		V_perturb_tp = NormalData(V_perturb.getData() * tp.getData())
		O_perturb_tp = NormalData(O_perturb_data * tp.getData())
		Tp_perturb = Tp_perturb + Step_acoustic*(Tp_S - (mapFactorM*mapFactorM*(U_perturb_tp.spatial(2).getData() + V_perturb_tp.spatial(1).getData())+mapFactorM* O_perturb_tp.spatial(0).getData()))
		 

		W_perturb = W_perturb + Step_acoustic*(W_S  + g/mapFactorM*(D/Ddry)*(p_perturb.spatial(0).getData() - timeAvg("Mdry_perturb", Mdry_perturb_data, acoustic start, Step_acoustic)))
		gp_perturb = gp_perturb + Step_acoustic*(gp_S - mapFactorM*( O_perturb_data * gp.spatial(0).getData() - g*timeAvg("W_perturb", W_perturb.getData(), acoustic start, Step_acoustic))/Mdry)
		 
		#diagnostic equation
		C = cs_sqr()/Mdry_perturb/Ddry^2
		Ddry_perturb_data = -(gp_perturb.spatial(0).getData() + Ddry.getData() * Mdry_perturb_data)/Mdry.getData()
		p_perturb = cs_sqr()/Ddry.getData() *(Tp_perturb/Tp - Ddry_perturb_data/Ddry.getData() - Mdry_perturb_data/Mdry.getData())

		acoustic_end = RK3_start + (j+1)*Step_acoustic

		Output("U_perturb", U_perturb, acoustic_end)
		Output("V_perturb", V_perturb, acoustic_end)
		Output("W_perturb", W_perturb, acoustic_end)
		Output("O_perturb", O_perturb, acoustic_end)
		Output("Mdry_perturb", Mdry_perturb, acoustic_end)
		Output("Tp_perturb", Tp_perturb, acoustic_end)
		Output("gp_perturb", gp_perturb, acoustic_end)
		Output("Ddry_perturb", Ddry_perturb, acoustic_end)
		Output("p_perturb", p_perturb, acoustic_end)

	U 		 = coupled(data_collector.collect("U", RK3_start )) - U_perturb		
	V 		 = coupled(data_collector.collect("V", RK3_start )) - V_perturb		
	W 		 = coupled(data_collector.collect("W", RK3_start )) - W_perturb		
	O 		 = coupled(data_collector.collect("O", RK3_start )) - O_perturb		
	Tp 		 = coupled(data_collector.collect("T", RK3_start )) - Tp_perturb		
	Mdry_pfr = data_collector.collect("PH", RK3_start ) 		- Mdry_perturb  	
	gp_pfr 	 = data_collector.collect("PH", RK3_start ) 		- gp_perturb 		
	Mdry 	 = Mdry_pfr + Mdry_refer
	gp 		 = gp_pfr 	+ gp_refer

	Output("U",	uncouple(U), acoustic_end)
	Output("V", uncouple(V), acoustic_end)
	Output("W", uncouple(W), acoustic_end)
	Output("O", uncouple(O), acoustic_end)
	Output("T", uncouple(Tp), acoustic_end)
	Output("MU", Mdry_pfr, acoustic_end)
	Output("PH", gp_pfr, acoustic_end)

def prepare():
	#confifuration
	idxK		= data_collector.collect("ZNU", t_now)
	idxW		= data_collector.collect("ZNW", t_now)
	lat 		= data_collector.collect("XLAT", t_now)
	lon 		= data_collector.collect("XLONG", t_now)
	mapFactorM  = data_collector.collect("MAPFAC_M", t_now)
	mapFactorU  = data_collector.collect("MAPFAC_U", t_now)
	mapFactorV  = data_collector.collect("MAPFAC_V", t_now)
	wdK 		= data_collector.collect("DN", t_now)
	wdW 		= data_collector.collect("DNW", t_now)
	wd 			= (wdK, wdY, wdK)
	#const
	e 			= ddata_collector.collect("E", t_now)
	f 			= data_collector.collect("F", t_now)
	sinalpha 	= data_collector.collect("SINALPHA", t_now)
	cosplpha 	= data_collector.collect("COSALPHA", t_now)
	#base state
	p_sd		= data_collector.collect("PSFC", t_now)
	p_dry_top 	= data_collector.collect("P_TOP", t_now)
	p0 			= data_collector.collect("P00", t_now)
	T0 			= data_collector.collect("T00", t_now)
	A  			= data_collector.collect("TLP", t_now)
	#reference state
	p_refer 	= data_collector.collect("PB", t_now)
	gp_refer 	= data_collector.collect("PHB", t_now)
	Mdry_refer 	= data_collector.collect("MUB", t_now)
	Ddry_refer 	= Ddry_reference()
	#perturbaton state from referenvce
	p_pfr 		= data_collector.collect("P", t_now)
	gp_pfr 		= data_collector.collect("PH", t_now)
	Mdry_pfr 	= data_collector.collect("MU", t_now)
	Ddry_pfr 	= Ddry_perturb_from_refer()
	# features participating in RK advance 
	u 			= data_collector.collect("U", t_now)
	v 			= data_collector.collect("V", t_now)
	w 			= data_collector.collect("W", t_now)
	o 			= 
	output("O", o, t_now)
	tp 			= data_collector.collect("T",t_now)
	qv 			= data_collector.collect("QVAPOR", t_now)
	qc 			= data_collector.collect("QCLOUD", t_now)
	qr 			= data_collector.collect("QRAIN", t_now)
	qi 			= data_collector.collect("QICE", t_now)
	qm 			= qv+qc+qr+qi
	output("Qm", qm, t_now)
	output("Qv", qv, t_now)
	output("Qc", qc, t_now)
	output("Qr", qr, t_now)
	output("Qi", qi, t_now)
	p 			= data_collector.collect("P_HYD", t_now) # = p_pfr + p_refer
	gp 			= gp_pfr 	+ 	gp_refer
	Mdry 		= Mdry_pfr 	+ 	Mdry_refer
	Ddry 		= Ddry_pfr 	+ 	Ddry_refer 
	D 			= Ddry/(1+qm)
	
 	# dry air mass coupled process
	U 			= coupled(u)
	V 			= coupled(v)
	W 			= coupled(w)
	O 			= coupled(o)
	Tp 			= coupled(tp)
	Qv 			= coupled(qv)
	Qc 			= coupled(qc) 
	Qr 			= coupled(qr) 
	Qi 			= coupled(qi)
	Qm 			= coupled(qm)

def BC(name, axis, pos, tend_opt):
	global t_now, num_bdy_layer
	opt  	= int(np.heaviside(tend_opt, 0))*"T"
	dxn  	= int(np.heaviside(axis-1 ,  0))*"X" + (1-int(np.heaviside(axis-1, 0)))*"Y"
	p  		= int(np.heaviside(pos,      0))*"E" + (1-int(np.heaviside(pos,    0)))*"S"
	vert  	= int(np.heaviside(1-axis,   0))
	if vert == 1:
		bc 		= data_collector.collect(name, ,t_now)
	elif:
		key  	= "%s_B%s%s%s"%(name, opt, dxn, p)
		bc 		= np.expand_dims(BC_collector.collect(key, t_now)[num_bdy_layer-1], axis)

	return bc

#Runge-Kutta Step	main-loop
def main():
	global t_now, t_start, t_end, RK3_start, RK3_end, Step_RK3, Step_acoustic, n, ns, mapFactorM, U, V, O, gp, D, Ddry, Mdry, g, W, qc, qr, qv, qm, p, p_refer, gp_refer, Mdry_refer, Ddry_refer, p_pfr, gp_pfr, Mdry_pfr, Ddry_pfr 

	for t in range(0, (t_end-t_start)/Step_RK3): 
		t_now = t_start + t*Step_RK3
			
		for i in range(3,0,-1):
			RK3_start = t_now

			if i == 3:
				n = 1, Step_acoustic = Step_RK3/3

				U_tend    = Adv("U")   +Fcor("U")
				V_tend    = Adv("V")   +Fcor("V")
				W_tend    = Adv("W")   +Fcor("W")
				Tp_tend   = Adv("Tp")  +Fcor("Tp")
				gp_tend   = Adv("gp")  +Fcor("gp")
				Mdry_tend = Adv("Mdry")+Fcor("Mdry")
				
			else:
				n = ns/i, Step_acoustic = Step_RK3/ns

			var_S = (
					ModeS("U", U_tend),
					ModeS("V", V_tend),
					ModeS("W", W_tend),
					ModeS("Mdry", Mdry_tend),
					ModeS("Tp", Tp_tend),
					ModeS("gp", gp_tend)
					)

			acoustic_time_integration(var_S, n, RK3_start, Step_acoustic)

			
			Qm_tend = - mapFactorM*mapFactorM*(Spatial(timeAvg("U", acoustic_start, Step_acoustic*n)*Interp(qm, 2), 2) + Spatial(timeAvg("V", acoustic_start, Step_acoustic*n)*Interp(qm, 1), 1)) - mapFactorM*Spatial(timeAvg("O", acoustic_start, Step_acoustic*n)*Interp(qm, 0), 0)
			Qm_advanced = coupled(data_collector.collect("Qm", RK3_start)) + Step_RK3*Qm_tend

			p_pfr 		= p_moist()
			Ddry_pfr 	= (Spatial(data_collector.collect("PH", acoustic_end), 0) + Mdry*data_collector.collect("MU", acoustic_end))/-Mdry_refer
			p 			= p_pfr  		+ p_refer
			Ddry 		= Ddry_pfr  	+ Ddry_refer

			RK3_end = RK3_start + n*Step_acoustic

			Output("Qm", uncouppled(Qm_advanced), RK3_end)
			Output("P", p_pfr, RK3_end)
			Output("P_HYD", p, RK3_end)
			Output("Ddry_pfr", Ddry_pfr, RK3_end)
			Output("Ddry", Ddry, RK3_end)

			