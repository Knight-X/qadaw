import numpy as np
from data_collection import AirDynamicData
import datetime
#Initialization
data_collector 		= AirDynamicData("climate2014-01.nc")
BC_collecctor 		=  AirDynamicData("wrfbdy_d01")

t_start 			= datetime.datetime() 
t_end 				= datetime.datetime()
t_now 				= t_start
RK3_start 			= None
RK3_end				= None
acoustic_start 		= None
acoustic_end 		= None
n 					= None
ns 					= None
Step 				= datetime.timedelta(seconds = ) # coorelated to interval seconds in wrf namelist
Step_RK3 			= datetime.timedelta(seconds = ) # 6*wdX sec is suggested
Step_acoustic		= None

wdX					= 1100
wdY 				= 1100
wdK 				= None
wd 					= None
nX 					= 200
nY   				= 300
nK					= 31 "equals to the numbers of vertical isobaric layers from global model"
num_bdy_layer 		= 1

idxK 				= None,
idxW 				= None,
lat 				= None,
lon 				= None,


p_dry_top = None,
p0 = None#pa
T0 = None #K, mean sea level temp
A = None #K, difference between p0 and p0/e
p_sd"grid"
gp_surface

g = 9.81
Rd = 287 #J/kg/K
Rv = 461.6 #J/kg/K
rEarth = 6.307*10^6 #meter
Beta = 0.1
e = None 
f = None 
sinalpha = None 
cosalpha = None
MW_wv = "molecular weight of water vapor"
MW_dry = "molecular weight of dry air"
epcilon = MW_wv/MW_dry

#input from wrfd01
mapFactorM = None,
mapFactorU = None,
mapFactorV = None,


p_refer = None, gp_refer = None, Mdry_refer = None, Ddry_refer = None,
p_pfr   = None, gp_pfr   = None, Mdry_pfr   = None, Ddry_pfr   = None,

u  = None, v  = None, w    = None, o  = None, tp   = None,
gp = None, p  = None, Mdry = None, D  = None, Ddry = None,
qv = None, qc = None, qr   = None, qi = None, qm   = None,

#coupled
cplst = ["U","V","W","O","Tp","Qm"]
U = None, V = None, W = None, O = None, Tp = None, Qv = None, Qc = None, Qr = None, Qi = None, Qm = None


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

def gp_reference_grad():
	grad = Ddry_reference()*Mdry_reference()
	return grad

def gp_reference():
	gpr = Grad_Interp(gp_reference_grad(), 0)

#perturbation state from reference
def Md():
	global p_sd
	dry_mass = p_sd - p_dry_surface()
	return dry_mass

def Dd():
	global Rd, Rv, p0, tp, p_pfr, p_refer
	dry_density = Rd/p0*tp*(1+Rv/Rd*qv)*(p_pfr + p_refer/p0)^(-cv()/cp())
	return dry_density

def Ddry_perturb_from_refer():
	Dd_pfr = Dd() - Ddry_reference()
	return Dd_pfr

def gp_perturb_from_refer_grad():
	grad = -(Md()*Ddry_perturb_from_refer()+Mdry_perturb_from_refer()*Ddry_reference())
	return grad

def gp_perturb_from_refer ():
	gp_pfr = Grad_Interp(gp_perturb_from_refer_grad(), 0)
	return gp_pfr 

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

	if level_determine(var) == "mp":
		M = Mdry, mF = 1
	elif level_determine(var) == "u":
		M = Interp(Mdry,2), mf = mapFactorU
	elif level_determine(var) == "v":
		M = Interp(Mdry,1), mf = mapFactorV
	elif level_determine(var) == "w":
		M = Mdry, mf = mapFactorM
	Var = M*var/mF

	return Var

def uncoupled(Var):
	global Mdry
	if level_determine(var) == "mp":
		M = Mdry, mF = 1
	elif level_determine(var) == "u":
		M = Interp(Mdry,2), mf = mapFactorU
	elif level_determine(var) == "v":
		M = Interp(Mdry,1), mf = mapFactorV
	elif level_determine(var) == "w":
		M = Mdry, mf = mapFactorM
	var = *mf*Var/M

	return var

#Operator
def timeAvg(name, var_advanced, start, Step):
	global Beta, cplst

	var_prev = data_collector.collect(name, start - step)
	var_timeAvg = (1+Beta)/2*var_advanced+(1-Beta)/2*var_prev

	return var_timeAvg


def Spatial(var, axis):#axis option (x,y)
	global wd, cplst
	

	if level_determine(var) == "mp":
		var_forward  = np.concatenate((var.take(np.arange(1, var.shape[axis]), axis), bc_forward   ), axis)
		var_backward = np.concatenate((bc_backward, var.take(np.arange(0, var.shape[axis]-1), axis)), axis)
	else:
		var_forward  = var.take(np.arange(1, var.shape[axis]  ), axis)
		var_backward = var.take(np.arange(0, var.shape[axis]-1), axis)

	Var = (var_forward - var_backward)/wd[axis]/2

	return Var

def Interp(name, var, axis):#x is supposed to be a tuple
	global wd, cplst


	if level_determine(var) == "mp":
		wd_forward  = wd[axis]
		wd_backward = wd[axis]
		
		var_forward  = np.concatenate((var, bc_forward),  axis)
		var_backward = np.concatenate((bc_backward, var), axis)
		
		Var = 0.5*(wd_backward/(0.5*(wd_backward + wd_forward))*var_forward + wd_forward/(0.5*(wd_backward + wd_forward))*var_backward)
	else:
		var_forward = var.take(np.arange(1, var.shape[axis]), axis)
		var_backward = var.take(np.arange(0, var.shape[axis]-1), axis)
		
		Var = 0.5*(var_forward + var_backward)	
	
	return Var

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
		Vxy = Interp(Interp(V,1),2)
		Wxk = Interp(Interp(W,0),2)
		tend = Vxy*(
			Interp(f,2)+Interp((Interp(uncouple(U),2)*Spatial(mapFactorM,1)-Interp(uncouple(V),1)*Spatial(mapFactorM,2)),2)
			)
			- Interp(e,2)*Wxk*Interp(cosplpha, 2)
			- uncouple(U)*Wxk/rEarth

	elif var == "V":
		Uxy = Interp(Interp(U,2),1)
		Wyk = Interp(Interp(W,0),1)
		tend = - Uxy*(
			Interp(f,1)+Interp((Interp(uncople(U),2)*Spatial(mapFactorM,1)-Interp(uncouple(V),1)*Spatial(mapFactorM,2)),1)
			)
			+ Interp(e,1)*Wyk*Interp(sinalpha, 1)
			- uncouple(V)*Wyk/rEarth

	elif var == "W":
		Uxk = Interp(Interp(U,2),0)
		Vyk = Interp(Interp(V,1),0)
		tend = e*(
			Uxk*cosplpha
			- Vyk*sinalpha
			)
			+ (
				Interp(Interp(uncouple(U),2),0)*Uxk
				+ Interp(Interp(uncouple(V),1),0)*Vyk
				)/rEarth

	else:
		tend = 0

	return tend

def Adv(var):
	global mapFactorM, mapFactorM, U, V, O, gp, acoustic_start, Step_acoustic, ns
	if var == "U":
		tend = - mapFactorM*(Spatial(U*uncoupled(U), 2) + Spatial(V*uncoupled(U), 1)) - Spatial(O*uncoupled(U), 0)
	elif var == "V":
		tend = - mapFactorM*(Spatial(U*uncoupled(V), 2) + Spatial(V*uncoupled(V), 1)) - mapFactorM/mapFactorM*Spatial(O*uncoupled(V), 0)
	elif var == "Mdry":
		tend = - mapFactorM*mapFactorM*(Spatial(U, 2) + Spatial(V, 1)) - mapFactorM*Spatial(O, 0)
	elif var == "Tp":
		tend = - mapFactorM*mapFactorM*(Spatial(U*Interp(uncoupled(Tp)), 2) + Spatial(V*Interp(uncoupled(Tp)), 1)) - mapFactorM*Spatial(O*uncoupled(Tp), 0)
	elif var == "W":
		tend = - mapFactorM*mapFactorM/mapFactorM*(Spatial(U*uncoupled(W), 2) + Spatial(V*uncoupled(W), 1)) - Spatial(O*uncoupled(W), 0)
	elif var == "gp":
		tend = - (mapFactorM*mapFactorM*(U*Spatial(gp, 2) + V*Spatial(gp, 1)) + mapFactorM*O*Spatial(gp, 0))/Mdry
	else:
		tend = 0

	return tend

#Time integration
def ModeS(var, tend):
	global mapFactorM, U, V, O, gp, D, Ddry, Mdry, g, W, qc, qr, qv, p, p_refer, gp_refer, Mdry_refer, Ddry_refer, p_pfr, gp_pfr, Mdry_pfr, Ddry_pfr 

	if var == "U":
		S = -(Interp(Mdry,2)*Interp(D,2)*Spatial(p_pfr,2)
				-Interp(Mdry,2)*Interp(Ddry_pfr,2)*Spatial(p_refer,2)) 
				-Interp(D/Ddry, 2)
				*(Interp(Mdry,2)*Spatial(Interp(gp_pfr,0),2)
					-Spatial(Interp(Interp(p_pfr,0),2),0)*Spatial(Interp(gp,0),2)
					+Interp(Mdry_pfr,2)*Spatial(Interp(gp,0),2))
			+ tend

	elif var == "V":
		S = -(Interp(Mdry,1)*Interp(D,1)*Spatial(p_pfr,1)
				-Interp(Mdry,1)*Interp(Ddry_pfr,1)*Spatial(p_refer,1)) 
				-Interp(D/Ddry, 1)
				*(Interp(Mdry,1)*Spatial(Interp(gp_pfr,0),1)
					-Spatial(Interp(Interp(p_pfr,0),1),0)*Spatial(Interp(gp,0),1)
					+Interp(Mdry_pfr,1)*Spatial(Interp(gp,0),1)) 
			+ tend

	elif var == "Mdry":
		S =  tend

	elif var == "Tp":
		S =  tend

	elif var == "W":
		S = g*(Interp(D/Ddry, 0)*(Spatial(p_pfr,0)+Mdry_refer*(qr+qc+qv))-Mdry_pfr)/mapFactorM 
			+ tend

	elif var == "gp":
		S = mapFactorM*g*W/Mdry
			+ tend
	
	return S


def acoustic_time_integration(var_S, n, RK3_start, Step_acoustic):
	global mapFactorM, U, V, O, gp, D, Ddry, Mdry, g, W, qc, qr, qv, p, p_refer, gp_refer, Mdry_refer, Ddry_refer, p_pfr, gp_pfr, Mdry_pfr, Ddry_pfr 

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
	p_perturb 	 = cs_sqr()/Ddry*(Tp_perturb/Tp - Ddry_perturb/Ddry - Mdry_perturb/Mdry)

	Output("U_perturb", U_perturb, RK3_start )
	Output("V_perturb", V_perturb, RK3_start )
	Output("W_perturb", W_perturb, RK3_start )
	Output("O_perturb", O_perturb, RK3_start )
	Output("Mdry_perturb", Mdry_perturb, RK3_start )
	Output("Tp_perturb", Tp_perturb, RK3_start )
	Output("gp_perturb", gp_perturb, RK3_start )
	Output("Ddry_perturb", Ddry_perturb, RK3_start )
	Output("p_perturb", p_perturb, RK3_start)

	for j in range(0, n):
		acoustic_start = RK3_start + j*Step_acoustic

		U_perturb = U_perturb + Step_acoustic*
		( U_S - 
			((mapFactorM/mapFactorM)*(D/Ddry)*(
				Mdry*(Ddry*Spatial(data_collector.collect("p_perturb", acoustic_start),2)+Spatial(p_refer,2)*data_collector.collect("Ddry_perturb", acoustic_start)+Spatial(data_collector.collect("gp_perturb", acoustic_start),2))
				+Spatial(gp, 2)*(Spatial(data_collector.collect("p_perturb", acoustic_start),0)-data_collector.collect("Mdry_perturb", acoustic_start)))
			)
		)

		V_perturb = V_perturb + Step_acoustic*
		(V_S - 
			((mapFactorM/mapFactorM)*(D/Ddry)*(
				Mdry*(Ddry*Spatial(data_collector.collect("p_perturb", acoustic_start),1)+Spatial(p_refer,1)*data_collector.collect("Ddry_perturb", acoustic_start)+Spatial(data_collector.collect("gp_perturb", acoustic_start),1))
				+Spatial(gp, 1)*(Spatial(data_collector.collect("p_perturb", acoustic_start),0)-data_collector.collect("Mdry_perturb", acoustic_start)))
			)
		)

		Mdry_perturb_difference = mapFactorM*mapFactorM*Vert_Integrat(Spatial(U_perturb, 2) + Spatial(V_perturb, 1))
		
		Mdry_perturb = Mdry_perturb + Step_acoustic*Mdry_perturb_difference 		  
		O_perturb = Grad_Interp(Mdry_S - Mdry_perturb_difference - mapFactorM*mapFactorM*(Spatial(U_perturb, 2) + Spatial(V_perturb, 1)))/mapFactorM
		Tp_perturb = Tp_perturb + Step_acoustic*(Tp_S - (mapFactorM*mapFactorM*(Spatial(U_perturb*tp,2)+Spatial(V_perturb*tp,1))+mapFactorM*Spatial(O_perturb*tp,0)))
		 

		W_perturb = W_perturb + Step_acoustic*(W_S  + g/mapFactorM*(D/Ddry)*(Spatial(p_perturb ,0)-timeAvg("Mdry_perturb", Mdry_perturb, acoustic start, Step_acoustic)))
		gp_perturb = gp_perturb + Step_acoustic*(gp_S - mapFactorM*( O_perturb * Spatial(gp, 0) - g*timeAvg("W_perturb", W_perturb, acoustic start, Step_acoustic))/Mdry)
		 
		#diagnostic equation
		C = cs_sqr()/Mdry_perturb/Ddry^2
		Ddry_perturb = -(Spatial(gp_perturb ,0) + Ddry*Mdry_perturb)/Mdry
		p_perturb = cs_sqr/Ddry*(Tp_perturb/Tp - Ddry_perturb/Ddry - Mdry_perturb/Mdry)

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
	Ddry_pfr 	= Dd()
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
	global mapFactorM, U, V, O, gp, D, Ddry, Mdry, g, W, qc, qr, qv, qm, p, p_refer, gp_refer, Mdry_refer, Ddry_refer, p_pfr, gp_pfr, Mdry_pfr, Ddry_pfr 

	for t in range(0, (t_end-t_start)/Step): 
		t_now = t_start + t*Step
			
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

			