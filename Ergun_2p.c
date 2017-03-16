#include "udf.h"
#include "mem.h"

#define EPS 1.0e-6
#define Dp 0.95e-3
#define Dw 1.0e-2						//diameter for water

DEFINE_PROFILE(vis_bed, t, i)
{
	//viscous resistance coefficient for bed fluid zone
	cell_t c;
	real epsilon;						//porosity
	real Kr;							//relative permeability
	real a;								//phase fraction: for vapor is alpha; for liquid is saturation = 1 - alpha
	real D;								//Fluent coefficient for permeability

	begin_c_loop(c, t)
	{
		epsilon = C_POR(c, t);
		a = C_VOF(c, t);				//phase fraction; the phase is determined by the zone that UDF is hooked
		Kr = pow(a, 3.0);				//Reed model
		if(Kr<EPS)
			Kr = EPS;

		D = 150.0*pow(1-epsilon, 2.0)/(pow(Dp, 2.0)*pow(epsilon, 3.0))/Kr;
		
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

DEFINE_PROFILE(vis_pool, t, i)
{
	//viscous resistance coefficient for pool zone
	cell_t c;
	real epsilon;						//porosity
	real Kr;							//relative permeability
	real a;								//phase fraction: for vapor is alpha; for liquid is saturation = 1 - alpha
	real D;								//Fluent coefficient for permeability

	begin_c_loop(c, t)
	{
		epsilon = C_POR(c, t);
		a = C_VOF(c, t);				//phase fraction; the phase is determined by the zone that UDF is hooked
		Kr = pow(a, 3.0);				//Reed model
		if(Kr<EPS)
			Kr = EPS;

		D = 150.0*pow(1-epsilon, 2.0)/(pow(Dw, 2.0)*pow(epsilon, 3.0))/Kr;		//Water diameter is used here

		F_PROFILE(c, t, i) = D;
	}
	end_c_loop(c, t)
}

DEFINE_PROFILE(iner_bed, t, i)
{
	//inertial resistance coefficient for bed zone
	cell_t c;
	real epsilon;						//porosity
	real eta_r;							//relative passability
	real a;								//phase fraction: for vapor is alpha; for liquid is saturation = 1 - alpha
	real C;								//Fluent coefficient for passability for liquid phase

	begin_c_loop(c, t)
	{
		epsilon = C_POR(c, t);
		a = C_VOF(c, t);				//phase fraction; the phase is determined by the zone that UDF is hooked
		eta_r = pow(a, 5.0);			//Reed model
		if(eta_r<EPS)
			eta_r = EPS;

		C = 2*1.75*(1-epsilon)/(Dp*pow(epsilon, 3.0))/eta_r;

		F_PROFILE(c, t, i) = C;
	}
	end_c_loop(c, t)
}

DEFINE_PROFILE(iner_pool, t, i)
{
	//inertial resistance coefficient for pool zone
	cell_t c;
	real epsilon;						//porosity
	real eta_r;							//relative passability
	real a;								//phase fraction: for vapor is alpha; for liquid is saturation = 1 - alpha
	real C;								//Fluent coefficient for passability for liquid phase

	begin_c_loop(c, t)
	{
		epsilon = C_POR(c, t);
		a = C_VOF(c, t);				//phase fraction; the phase is determined by the zone that UDF is hooked
		eta_r = pow(a, 5.0);			//Reed model
		if(eta_r<EPS)
			eta_r = EPS;

		C = 2*1.75*(1-epsilon)/(Dw*pow(epsilon, 3.0))/eta_r;			//Water diameter is used here

		F_PROFILE(c, t, i) = C;
	}
	end_c_loop(c, t)
}

