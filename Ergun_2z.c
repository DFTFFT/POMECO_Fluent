#include "udf.h"
#include "mem.h"

//zone index
#define bed_large_fluid_zone_ID 5
#define bed_large_solid_zone_ID 37
#define bed_small_fluid_zone_ID 6
#define bed_small_solid_zone_ID 24
#define pool_zone_ID 4

//
#define EPS 1.0e-6						//error tolerance
#define Dps 1.5e-3						//diameter of small particle 
#define Dpl 3.0e-3						//diameter of large particle 
#define Dw 1.0e-2						//diameter for water



DEFINE_PROFILE(vis_coef, t, i)
{
	//viscous resistance coefficient for various zones
	//Reed model is used

	cell_t c;
	real epsilon;						//porosity
	real Kr;							//relative permeability
	real a;								//phase fraction: for vapor is alpha; for liquid is saturation (sat = 1 - alpha)
	real dia;							//diameter of particles in this zone
	real D;								//Fluent coefficient for permeability

	int zone_id = THREAD_ID(t);			//obtain the index of current zone to determine the particle size

	//Determine the particle size in a particular zone
	switch(zone_id)
	{
		case bed_small_fluid_zone_ID:
			dia = Dps;
			break;

		case bed_large_fluid_zone_ID:
			dia = Dpl;
			break;

		case pool_zone_ID:
			dia = Dw;
			break;

		default:
			Error("Error in the particle diameter.\n");
			break;
	}

	//Message("Zone ID = %d. Diameter = %g\n", zone_id, dia);

	//
	begin_c_loop(c, t)
	{
		epsilon = C_POR(c, t);
		a = C_VOF(c, t);				//phase fraction; the phase is determined by the zone that UDF is hooked
		Kr = pow(a, 3.0);				//Reed model

		//avoid being divided by zero
		if(Kr < EPS)						
			Kr = EPS;

		D = 150.0*pow(1.0 - epsilon, 2.0)/(dia*dia*pow(epsilon, 3.0))/Kr;
		
		/*
		Message("\n");
		Message("porosity = %f\n", epsilon);
		Message("a = %f\n", a);
		Message("Kr = %f\n", Kr);
		Message("D = %f\n", D);
		Message("\n");
		*/
		
		F_PROFILE(c, t, i) = D;
	}
	end_c_loop(c, t)
}



DEFINE_PROFILE(iner_coef, t, i)
{
	//inertial resistance coefficient for various zones
	//Reed model is used

	cell_t c;
	real epsilon;						//porosity
	real eta_r;							//relative passability
	real a;								//phase fraction: for vapor is alpha; for liquid is saturation (sat = 1 - alpha)
	real dia;							//diameter of particles in this zone
	real C;								//Fluent coefficient for passability for liquid phase


	int zone_id = THREAD_ID(t);			//obtain the index of current zone to determine the particle size

	//Determine the particle size in a particular zone
	switch(zone_id)
	{
		case bed_small_fluid_zone_ID:
			dia = Dps;
			break;

		case bed_large_fluid_zone_ID:
			dia = Dpl;
			break;

		case pool_zone_ID:
			dia = Dw;
			break;

		default:
			Error("Error in the particle diameter.\n");
			break;
	}

	//Message("Zone ID = %d. Diameter = %g\n", zone_id, dia);

	//
	begin_c_loop(c, t)
	{
		epsilon = C_POR(c, t);
		a = C_VOF(c, t);				//phase fraction; the phase is determined by the zone that UDF is hooked
		eta_r = pow(a, 5.0);			//Reed model

		//avoid being divided by zero
		if(eta_r < EPS)
			eta_r = EPS;

		C = 2*1.75*(1.0 - epsilon)/(dia*pow(epsilon, 3.0))/eta_r;

		F_PROFILE(c, t, i) = C;
	}
	end_c_loop(c, t)
}
