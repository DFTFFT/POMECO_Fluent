#include  "udf.h"
#include "stdlib.h"

//zone index
#define bed_fluid_zone_ID 2
#define bed_solid_zone_ID 4

//domain index
#define mixture_domain_ID 1
#define liq_domain_ID 2			//obtained from GUI
#define vap_domain_ID 3			//obtained from GUI

//phase domain index
#define liq_phase_ID 0			//primary phase
#define vap_phase_ID 1			//secondary phase

//
#define Dp 0.98e-3				//particle diameter, m
#define g 9.81					//gravitational acceleration, m/s2


//for direct heat transfer between solid and fluid (liquid & vapor)
real Fd(real a);
real Area_dir(real por, real a);
real HTC_dir(cell_t c, Thread *t);

//for heat transfer involing interface
real Fi(real sl);
real Area_i(real por, real sl);
real HTC_boiling(cell_t c, Thread *t_liq, Thread *t_vap, real Tsat, real Ts);
real HTC_PB(cell_t c, Thread *t_liq, Thread *t_vap, real Tsat, real Ts);
real HTC_FB(cell_t c, Thread *t_liq, Thread *t_vap, real Tsat, real Ts);

//for bubbly flow
real Area_B(real D_B, real por, real sat);
void HTC_bubbly(cell_t c, Thread *t, real D_B, real *h_l, real *h_g);

//for droplet flow
real Area_D(real D_D, real por, real sat);
void HTC_droplet(cell_t c, Thread *t, real D_D, real *h_l, real *h_g);

//heat transfer power density between interface and fluid
void source_intf_fluid(cell_t c, Thread *t, real *Qil, real *Qig);

//fitting function of saturated temperature (K) with respect to fluid mixture pressure (Pa)
real Tsat_fit(real p);

//
real Decay_power();


//
DEFINE_ADJUST(init, domain)
{
	//save the solid temperature field to UDM-0 at the begining of every iteration
	Thread *t;
	cell_t c;

	t = Lookup_Thread(domain, bed_solid_zone_ID);

	begin_c_loop_int(c, t)
	{
		//Save the solid temperature field (UDS results) to UDM at the begining of every iteration
		C_UDMI(c, t, 0) = C_UDSI(c, t, 0);
	}
	end_c_loop_int(c, t)
}

/******************************************************/
// direct heat transfer (solid-liquid & solid-vapor)  //
/******************************************************/

DEFINE_SOURCE(source_sl, c, t, dS, eqn)
{
	//heat density exchanged directly between liquid and solid particles, W/m3
	Domain *domain;
	Thread *t_solid;

	real source;
    real Tl;
    real Ts;
    real por;
    real sat;											//liquid saturation (i.e. volume fraction), sat = 1 - alpha

    real Als;											//interfacial area density, 1/m
    real hls;											//heat transfer coefficent between liquid and solid particles, W/(m2-k)

    //
    Tl = C_T(c, t);										//liquid temperature
    por = C_POR(c, t);									//porosity
    sat = C_VOF(c, t);									//liquid saturation					
    
    Message("sat=%g\n", sat);

	domain = Get_Domain(mixture_domain_ID);
	t_solid = Lookup_Thread(domain, bed_solid_zone_ID);

    Ts = C_UDMI(c, t_solid, 0);
    
    //
    //Message("Liquid temperature = %g\n", Tl);
    //Message("Solid temperature = %g\n", Ts);
    
    //
    Als = Area_dir(por, sat);
    hls = HTC_dir(c, t);
    
    //Message("Als=%g\n", Als);
    //Message("hls=%g\n", hls);
    //

    source = -hls*Als*(Tl - Ts);

    //Message("source from liquid=%g\n", source);

    //source for solid zone
    C_UDMI(c, t, 1) = -source;

    dS[eqn] = 0.0;

    return source;
}


DEFINE_SOURCE(source_sg, c, t, dS, eqn)
{
	//heat density exchanged directly between vapor and solid particles, W/m3
	Domain *domain;
	Thread *t_solid;

	real source;
    real Tg;
    real Ts;
    real por;
    real alpha;											//vapor saturation (i.e. volume fraction)

    real Ags;											//interfacial area density, 1/m
    real hgs;											//heat transfer coefficent between vapor and solid particles, W/(m2-k)

    //
    Tg = C_T(c, t);										//vapor temperature
    por = C_POR(c, t);									//porosity
    alpha = C_VOF(c, t);								//vapor saturation					
    
    //Message("alpha=%g\n", alpha);

	domain = Get_Domain(mixture_domain_ID);
	t_solid = Lookup_Thread(domain, bed_solid_zone_ID);

    Ts = C_UDMI(c, t_solid, 0);
    
    //
    //Message("Vapor temperature = %g\n", Tg);
    //Message("Solid temperature = %g\n", Ts);
    
    //
    Ags = Area_dir(por, alpha);
    hgs = HTC_dir(c, t);
    
    //Message("Ags=%g\n", Ags);
    //Message("hgs=%g\n", hgs);
    //

    source = -hgs*Ags*(Tg - Ts);

    //Message("source from vapor=%g\n", source);

    //source for solid zone
    C_UDMI(c, t, 2) = -source;

    dS[eqn] = 0.0;

    return source;
}


DEFINE_SOURCE(source_solid, c, t, dS, eqn)
{
	//heat source density for solid particles bed, W/m3
	//consists of four parts:
	//		1. transferred directly to liquid when liquid phase is dominant;
	//		2. transferred directly to vapor when vapor phase is dominant;
	//		3. extracted by two-phase interface when fluid is saturated;
	//		4. external decay heat power.

	Domain *domain;
	Thread *t_mix;					//thread for mixture
    Thread *t_liq;					//thread for liquid phase
	Thread *t_vap;					//thread for vapor phase

    real source;					//return source result, W/m3
    real Quds;						//total heat source density for UDS, W/m3 
    real Qsl, Qsg;					//solid <-> liquid, solid <-> vapor, W/m3 
    real Qsi;						//solid <-> interface, W/m3 
    real Qdecay;					//internal decay heat power, W/m3

    real P_op;						//operating pressure, Pa
    real Pmix;						//fluid mixture pressure (absolute), Pa
    real por;						//porosity of fluid zone
    real sat;						//volume fraction of liquid
    real alpha;						//volume fraction of vapor; alpha = 1 - sat
    real Cp_s;						//specific heat capacity of solid, J/(kg-K)
    real Ts;						//solid particles temperature, K
    real Tsat;						//saturation temperature of fluid corresponding to Pmix, K
    real Asi;						//interfacial area density between solid and interface, 1/m
    real hsi;						//HTC between solid and interface

    //Obtain thread pointers for fluid zone of each phase (mixture, liquid and vapor)
    domain = Get_Domain(mixture_domain_ID);
	t_mix = Lookup_Thread(domain, bed_fluid_zone_ID);
	t_liq = THREAD_SUB_THREAD(t_mix, liq_phase_ID); 
	t_vap = THREAD_SUB_THREAD(t_mix, vap_phase_ID);

	//solid particles properties
	Ts = C_UDMI(c, t, 0);
    Cp_s = C_CP(c, t);

    //fluid mixture properties
 	P_op = RP_Get_Real("operating-pressure");
	Pmix = C_P(c, t_mix) + P_op;    
	Tsat = Tsat_fit(Pmix);		//saturation temperature of fluid computed by fitting function
								//getsatvalues_NIST(index, Pmix, y) doesn't work

	//liquid properties
	por = C_POR(c, t_liq);
	sat = C_VOF(c, t_liq);

	//vapor properties
	alpha = C_VOF(c, t_vap);

    //source_1: direct heat exchange between liquid and solid
	Qsl = C_UDMI(c, t_liq, 1);


	//source_2: direct heat exchange between vapor and solid
	Qsg = C_UDMI(c, t_vap, 2);

	//source_3: heat exchange power between solid and interface (only exists when alpha > 0)
	if(alpha > 0.0)
	{
		Asi = Area_i(por, sat);
		hsi = HTC_boiling(c, t_liq, t_vap, Tsat, Ts);
		Qsi = -hsi*Asi*(Ts - Tsat);
	}
	else
		Qsi = 0.0;

	//source_4: decay heat power
	Qdecay = Decay_power();
	//Message("Qdecay=%g\n", Qdecay);

	//total heat source density for UDS
    Quds = (Qdecay + (Qsl + Qsg +Qsi)/(1.0 - por))/Cp_s;
    C_UDMI(c, t, 3) = Quds;
	//Message("Quds=%g\n", Quds);

    source = Quds;

    dS[eqn] = 0.0;

    return source;
}


/*************************************************************************/
//  Heat trasnfer between solid and fluid (liquid & vapor) via interface //
/*************************************************************************/

DEFINE_SOURCE(source_il, c, t, dS, eqn)
{
	//heat density power exchanged between interface and bulk liquid, W/m3
	real source;
	real Qil; 						//heat density transferred between interface and bulk liquid, W/m3
	real Qig;						//heat density transferred between interface and bulk vapor, W/m3

	source_intf_fluid(c, t, &Qil, &Qig);

	source = -Qil;

	dS[eqn] = 0.0;

	return source;
}


DEFINE_SOURCE(source_ig, c, t, dS, eqn)
{
	//heat density power exchanged between interface and bulk vapor, W/m3
	real source;
	real Qil;						//heat density transferred between interface and bulk liquid, W/m3
	real Qig;						//heat density transferred between interface and bulk vapor, W/m3

	source_intf_fluid(c, t, &Qil, &Qig);

	source = -Qig;

	dS[eqn] = 0.0;

	return source;
}



/*******************************/
//  Functions implementations  //
/*******************************/

real Fd(real a)
{
	const real THD = 0.7;		//Threshold value
	real f;

	if(a<THD)
		f = 0.0;
	else
		f = (a - THD)/(1.0-THD);

	//Message("f=%g\n", f);

	return f;
}


real Fi(real sl)
{
	const real THD = 0.3;					//Threshold value
	real f;

	if(sl<THD)
		f = sl/THD;
	else
		f = 1.0;

	//Message("f=%g\n", f);

	return f;
}


real Area_dir(real por, real a)
{
	real area;

	area = 6.0*(1-por)/Dp*Fd(a);
	
	return area;
}


real Area_i(real por, real sl)
{
	real area;

	area = 6.0*(1-por)/Dp*Fi(sl);

	//Message("area=%g\n", area);
	
	return area;
}


real Area_B(real D_B, real por, real sat)
{
	//for bubbly flow
	real area;

	area = 6.0*por*(1 - sat)/D_B;

	return area;
}


real Area_D(real D_D, real por, real sat)
{
	//for droplet flow
	real area;

	area = 6.0*por*sat/D_D;

	return area;
}


real HTC_dir(cell_t c, Thread *t)
{
	//heat transfer coefficient for direct heat exchange between solid and fluid (liquid and vapor)
	//when one fluid phase is continuous/dominant while another is dispersed phase
	//Ranz-Marshall correlation

	real h;						//heat transfer coefficient between fluid and solid, W/(m2-K)
	real Nu;					//Nusselt number
	real Re;					//Reynolds number
	real k;						//Thermal conductivity, W/(m-K)
	real rho;					//density, kg/m3
	real u;						//u velocity, m/s
	real v;						//v velocity, m/s
	real vel;					//velocity (magnitude), m/s
	real mu_L;					//laminar viscosity, kg/(m-s)
	real mu_T;					//turbulent viscosity, kg/(m-s)
	real mu;					//dynamic viscosity, kg/(m-s

	//fluid properties
	rho = C_R(c, t);
	u = C_U(c, t);
	v = C_V(c ,t);
	mu_L = C_MU_L(c, t);		//laminar dynamic viscosity
	mu_T = C_MU_T(c, t);		//turbulent dynamic viscosity
	k = C_K_L(c, t);

	mu = mu_L;					//use laminar viscosity? 
	vel = sqrt(u*u + v*v);		//magnitude of velocity vector
	Re = fabs(rho*vel*Dp)/mu;

	Nu = 2.0 + 0.6*sqrt(Re);

	h = Nu*k/Dp;

	//
	/*
	Message("rho=%g\n", rho);
	Message("u=%g\n", u);
	Message("v=%g\n", v);
	Message("vel=%g\n", vel);
	Message("k=%g\n", k);
	Message("mu_L=%g\n", mu_L);
	Message("Re=%g\n", Re);
	Message("Nu=%g\n", Nu);
	Message("h=%g\n", h);
	Message("\n");
	*/

	return h;
}


real HTC_boiling(cell_t c, Thread *t_liq, Thread *t_vap, real Tsat, real Ts)
{
	//heat transfer coefficient for heat transfer between solid and interface
	//when the fluid is two-phase
	//physical mechanism: boiling heat transfer
	//		pure pool boiling: Rohsenow correlation
	//		transition region: interpolation
	//		pure film boiling: Lienhard correlation

	real h;						//boiling heat transfer coefficient, W/(m2-K)
	real Tmin_FB;				//minimum film boiling temperature, K
	real Tmax_PB;				//maximum pool boiling temperature, K

	real h_PB;					//HTC for pool boiling condition, W/(m2-K)
	real h_FB;					//HTC for pool film condition, W/(m2-K)

	real W;						//weight factor (depends on Ts)
								//W=0: pool boiling
								//W=1: film boiling

	//
	Tmin_FB = Tsat + 17.0;
	Tmax_PB = Tsat + 100.0;

	if(Ts < Tmin_FB)
	{
		//pure pool boiling, Rohsenow correlation
		W = 0;
	}
	else if(Ts > Tmax_PB)
	{
		//pure film boiling, Lienhard correlation
		W = 1;
	}
	else
	{
		//transition region
		W = (Ts - Tmin_FB)/(Tmax_PB - Tmin_FB);
	}

	h_PB = HTC_PB(c, t_liq, t_vap, Tsat, Ts);
	h_FB = HTC_FB(c, t_liq, t_vap, Tsat, Ts);

	h = (1.0 - W)*h_PB + W*h_FB;

	return h;
}


real HTC_PB(cell_t c, Thread *t_liq, Thread *t_vap, real Tsat, real Ts)
{
	//pool boiling, Rohsenow correlation
	const real Csf = 0.013;
	const real m = 0.33;
	const real n = 1.7;

	real h;									//HTC, W/(m2-K)
	real Cp_l;								//specific heat capacity of liquid, J/(kg-K)
	real miu_l;								//dynamic viscosity of liquid, kg/(m-s)
	real k_l;								//conductivity of liquid, W/(m-K)
	real rho_l;								//saturated liquid density, kg/m3
	real rho_g;								//saturated vapor density, kg/m3
	real Hl;								//saturated liquid enthalpy, J/kg
	real Hg;								//saturated vapor enthalpy, J/kg
	real Hfg;								//evaporation/latent heat, J/kg
	real sigma;								//surface tension, N/m
	real Pr;								//Prandtl number
	real dT;								//temperature difference (Ts - Tsat)
	real dRho;								//density difference (rho_sat_l - rho_sat_g)
	real buoy;								//buoyance term

	Cp_l = C_CP(c, t_liq);
	miu_l = C_MU_L(c, t_liq);				//laminar dynamic viscosity
	k_l = C_K_L(c, t_liq);

	//For the sake of simplicity, following saturated properties are set to be constant obtained from other application
	//saturation condition P = 1.01325e5 Pa
	rho_l = 958.3727;
	rho_g = 0.597623;
	Hl =  418.99e3;
	Hg = 2675.53e3;
	sigma = 0.0589;							//liquid water of 100C

	Hfg = Hg - Hl;
	Pr = Cp_l*miu_l/k_l;
	dT = Ts - Tsat;
	dRho = rho_l - rho_g; 
	buoy = sigma/(g*dRho);

	h = pow(Csf, -3.0)*pow(Cp_l, 3.0)*miu_l*(dT*dT)/(Hfg*Hfg*pow(buoy, 0.5)*pow(Pr, n*3.0));

	return h;
}


real HTC_FB(cell_t c, Thread *t_liq, Thread *t_vap, real Tsat, real Ts)
{
	//film boiling, Lienhard correlation
	const real C = 0.67;
	const int n = 0.25;

	real h;									//HTC, W/(m2-K)
	real miu_g;								//dynamic viscosity of vapor, kg/(m-s)
	real k_g;								//conductivity of vapor, W/(m-K)
	real rho_l;								//saturated liquid density, kg/m3
	real rho_g;								//saturated vapor density, kg/m3
	real Hl;								//saturated liquid enthalpy, J/kg
	real Hg;								//saturated vapor enthalpy, J/kg
	real Hfg;								//evaporation/latent heat, J/kg
	real Hfg_m;								//modified latent heat, J/kg
	real Ja;								//Jacobs number
	real Cp_g;								//specific heat capacity of vapor, J/(kg-K)
	real Pr_g;								//Prandtl number for vapor
	real dT;								//temperature difference (Ts - Tsat)
	real dRho;								//density difference (rho_sat_l - rho_sat_g)
	real Nu;								//Nusselt number

	real x;									//interim variable

	miu_g = C_MU_L(c, t_vap);				//laminar dynamic viscosity
	k_g = C_K_L(c, t_vap);
	Cp_g = C_CP(c, t_vap);

	//For the sake of simplicity, following saturated properties are set to be constant obtained from other application
	//saturation condition P = 1.01325e5 Pa
	rho_l = 958.3727;
	rho_g = 0.597623;
	Hl =  418.99e3;
	Hg = 2675.53e3;

	Hfg = Hg - Hl;
	dT = Ts - Tsat;
	dRho = rho_l - rho_g;

	//modified latent heat
	Pr_g = Cp_g*miu_g/k_g;
	Ja = Cp_g*dT/Hfg;
	Hfg_m = Hfg*(1 + (0.968 - 0.163/Pr_g)*Ja);

	//
	x = rho_g*dRho*g*Hfg_m*pow(Dp, 3.0)/(miu_g*k_g*dT);

	Nu = C*pow(x, n);

	h = Nu*k_g/Dp;

	return h;
}


void HTC_bubbly(cell_t c, Thread *t, real D_B, real *h_l, real *h_g)
{
	//HTC for heat transfer between interface and fluid under bubbly flow condition
	Thread *t_mix;
	Thread *t_liq;
	Thread *t_vap;

	real k_l;								//conductivity of liquid, W/(m-K)
	real k_g;								//conductivity of vapor, W/(m-K)
	real rho_l;								//liquid density, kg/m3
	real Cp_l;								//specific heat capacity of liquid, J(kg-K)
	real miu_l;								//dynamic viscosity of liquid (laminar), kg/(m-s)
	real Re_rel_l;							//Reynolds number of liquid with respect to relative velocity
	real Pr_l;								//Prandtl number of liquid
	real Nu_l;								//Nusselt number of liquid
	real Nu_g;								//Nusselt number of vapor

	real ul, ug;							//u-component velocity of liquid and vapor, m/s
	real vl, vg;							//v-component velocity of liquid and vapor, m/s
	real urel;								//relative velocity in x direction, m/s
	real vrel;								//relative velocity in y direction, m/s
	real vel_rel;							//magnitude of relative velocity vectorbetween phases, m/s
											//vel_rel = ((ul - ug)^2 + (vl - vg)^2 + (wl - wg)^2)^0.5

	//
	t_mix = THREAD_SUPER_THREAD(t);
	t_liq = THREAD_SUB_THREAD(t_mix, liq_phase_ID); 
	t_vap = THREAD_SUB_THREAD(t_mix, vap_phase_ID);

	//liquid properties
	rho_l = C_R(c, t_liq);
	Cp_l = C_CP(c, t_liq);
	k_l = C_K_L(c, t_liq);
	miu_l = C_MU_L(c, t_liq);				//laminar dynamic viscosity

	//vapor properties
	k_g = C_K_L(c, t_vap);

	//for liquid
	ul = C_U(c, t_liq);
	ug = C_U(c, t_vap);
	vl = C_V(c, t_liq);
	vg = C_V(c, t_vap);
	urel = ul - ug;
	vrel = vl - vg;
	vel_rel = sqrt(urel*urel + vrel*vrel);
	Re_rel_l = rho_l*vel_rel*D_B/miu_l;
	Pr_l = Cp_l*miu_l/k_l;
	Nu_l = 2.0 + 0.6*sqrt(Re_rel_l)*pow(Pr_l, 0.33333);

	//for vapor
	Nu_g = 10.0;

	//
	*h_l = Nu_l*k_l/D_B;
	*h_g = Nu_g*k_g/D_B;
}


void HTC_droplet(cell_t c, Thread *t, real D_D, real *h_l, real *h_g)
{
	//HTC for heat transfer between interface and fluid under droplet flow condition
	Thread *t_mix;
	Thread *t_liq;
	Thread *t_vap;

	real k_l;								//conductivity of liquid, W/(m-K)
	real k_g;								//conductivity of vapor, W/(m-K)
	real rho_g;								//vapor density, kg/m3
	real Cp_g;								//specific heat capacity of vapor, J(kg-K)
	real miu_g;								//dynamic viscosity of vapor (laminar), kg/(m-s)
	real Re_rel_g;							//Reynolds number of vapor with respect to relative velocity
	real Pr_g;								//Prandtl number of vapor
	real Nu_l;								//Nusselt number of liquid
	real Nu_g;								//Nusselt number of vapor

	real ul, ug;							//u-component velocity of liquid and vapor, m/s
	real vl, vg;							//v-component velocity of liquid and vapor, m/s
	real urel;								//relative velocity in x direction, m/s
	real vrel;								//relative velocity in y direction, m/s
	real vel_rel;							//magnitude of relative velocity vectorbetween phases, m/s
											//vel_rel = ((ul - ug)^2 + (vl - vg)^2 + (wl - wg)^2)^0.5


	//
	t_mix = THREAD_SUPER_THREAD(t);
	t_liq = THREAD_SUB_THREAD(t_mix, liq_phase_ID); 
	t_vap = THREAD_SUB_THREAD(t_mix, vap_phase_ID);

	//liquid properties
	k_l = C_K_L(c, t_liq);

	//vapor properties
	rho_g = C_R(c, t_vap);
	Cp_g = C_CP(c, t_vap);
	k_g = C_K_L(c, t_vap);
	miu_g = C_MU_L(c, t_vap);				//laminar dynamic viscosity

	//for liquid
	Nu_l = 10.0;

	//for vapor
	ul = C_U(c, t_liq);
	ug = C_U(c, t_vap);
	vl = C_V(c, t_liq);
	vg = C_V(c, t_vap);
	urel = ul - ug;
	vrel = vl - vg;
	vel_rel = sqrt(urel*urel + vrel*vrel);
	Re_rel_g = rho_g*vel_rel*D_D/miu_g;
	Pr_g = Cp_g*miu_g/k_g;
	Nu_g = 2.0 + 0.738*sqrt(Re_rel_g)*pow(Pr_g, 0.33333);

	//
	*h_l = Nu_l*k_l/D_D;
	*h_g = Nu_g*k_g/D_D;
	
}


void source_intf_fluid(cell_t c, Thread *t, real *Qil, real *Qig)
{
	//heat transfer power density exchanged between interface and bulk fluid, W/m3

	const real slim_B = 0.7;					//lower boundary value of liquid saturation for bubbly flow
	const real slim_D = 0.3;					//upper boundary value of liquid saturation for droplet flow

	Thread *t_mix;								//thread for mixture
	Thread *t_liq;								//thread for liquid phase
	Thread *t_vap;								//thread for vapor phase
									
	real Qil_B;									//heat density transferred between interface and bulk liquid (bubbly flow), W/m3
	real Qil_D;									//heat density transferred between interface and bulk liquid (droplet flow), W/m3
	real Qig_B;									//heat density transferred between interface and bulk vapor (bubbly flow), W/m3
	real Qig_D;									//heat density transferred between interface and bulk vapor (droplet flow), W/m3	
	real Tl;									//liquid temperature, K
	real Tg;									//vapor temperature, K
	real Tsat;									//fluid saturation temperature, K
	real area_B;								//interfacial area density for bubbly flow regime, 1/m
	real area_D;								//interfacial area density for droplet flow regime, 1/m	
	real hl_B;									//HTC of liquid for bubbly flow regime, W/(m2-K)
	real hl_D;									//HTC of liquid for droplet flow regime, W/(m2-K)
	real hg_B;									//HTC of vapor for bubbly flow regime, W/(m2-K)
	real hg_D;									//HTC of vapor for droplet flow regime, W/(m2-K)
	real D_B;									//bubble diameter, m
	real D_D;									//droplet diameter, m

	real Z_B;									//weight factor for bubbly flow
	real Z_D;									//weight factor for droplet flow
	real sat;									//volume fraction (saturation) of liquid
	real P_op;									//operating pressure, Pa
	real Pmix;									//fluid mixture pressure (absolute), Pa
	real por;									//porosity

	//
	t_mix = THREAD_SUPER_THREAD(t);
	t_liq = THREAD_SUB_THREAD(t_mix, liq_phase_ID); 
	t_vap = THREAD_SUB_THREAD(t_mix, vap_phase_ID);

	sat = C_VOF(c, t_liq);

	if(sat > 0.0)
	{
		//only exists when it is two-phase flow condition, in which heat is transferred betweeen solid and fluid via interface
		//flow pattern:
		//		- Bubbly flow
		//		- Transition region
		//		- Droplet flow

		Tl = C_T(c, t_liq);
		Tg = C_T(c, t_vap);
		por = C_POR(c, t_liq);
		P_op = RP_Get_Real("operating-pressure");
		Pmix = C_P(c, t_mix) + P_op;    
		Tsat = Tsat_fit(Pmix);					//saturation temperature of fluid computed by fitting function
												//getsatvalues_NIST(index, Pmix, y) doesn't work

		if(sat > slim_B)
		{
			//bubbly flow regime
			D_B = 0.125*Dp;
			area_B = Area_B(D_B, por, sat);
			HTC_bubbly(c, t, D_B, &hl_B, &hg_B);

			*Qil = area_B*hl_B*(Tl - Tsat);
			*Qig = area_B*hg_B*(Tg - Tsat);

		}
		else if(sat < slim_D)
		{
			//droplet flow regime
			D_D = 0.125*Dp;
			area_D = Area_D(D_D, por, sat);
			HTC_droplet(c, t, D_D, &hl_D, &hg_D);

			*Qil = area_D*hl_D*(Tl - Tsat);
			*Qig = area_D*hg_D*(Tg - Tsat);
		}
		else
		{
			//transition region		
			//weight factor
			Z_B = (sat - slim_D)/(slim_B - slim_D);
			Z_D = 1.0 - Z_B;

			//bubbly flow
			D_B = 0.125*Dp;
			area_B = Area_B(D_B, por, slim_B)*Z_B;			//Note: the limit value is used instead!
			HTC_bubbly(c, t, D_B, &hl_B, &hg_B);
			Qil_B = area_B*hl_B*(Tl - Tsat);
			Qig_B = area_B*hg_B*(Tg - Tsat);

			//droplet flow
			D_D = 0.125*Dp;
			area_D = Area_D(D_D, por, slim_D)*Z_D;			//Note: the limit value is used instead!
			HTC_droplet(c, t, D_D, &hl_D, &hg_D);
			Qil_D = area_D*hl_D*(Tl - Tsat);
			Qig_D = area_D*hg_D*(Tg - Tsat);

			//total heat power
			*Qil = Qil_B + Qil_D;
			*Qig = Qig_B + Qig_D;
		}
	}
	else
	{
		//Fluid is subcooling single liquid phase. No interface exists.
		*Qil = 0.0;
		*Qig = 0.0;
	}
}


real Decay_power()
{
	//The decay heat power is defined as the volumetric power density imposed only on the solid particles, W/m3
	//Receive the value from a user-defined parameter in the workbench
	real Ss;

	Ss = RP_Get_Input_Parameter("real-1");
	//Message("Ss=%g\n", Ss);

	return Ss;
}


real Tsat_fit(real p)
{
	//calculate the saturation tempearture of water (K) based on absolute pressure (Pa)
	//use linear polynomial function

	real Tsat;

	Tsat = 353.0 + 0.0002022*p;

	return Tsat;
}
