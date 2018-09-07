
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef MAIN_H_
#define MAIN_H_
/*** MACROS ***/
//Network Parameters///////////////////////////////////////////////////////////////////////////////////////////
#define RAND_INIT 0 // make it zero to facilitate debugging
#define SIMTIME 1// in ms, for when no input file is provided
//IO network size is IO_NETWORK_DIM1*IO_NETWORK_DIM2 times the TIME_MUX_FACTOR

//Change the Time Mux Factor to change the number of simulated cells in the network
#define TIME_MUX_FACTOR 1
#define HW_CELLS 1

#define IO_NETWORK_SIZE HW_CELLS*TIME_MUX_FACTOR
//Connection Matrix size in ints.
#define CONN_MATRIX_SIZE IO_NETWORK_SIZE*IO_NETWORK_SIZE

//Maximum Number of cells and tiTIME_MUX_FACTORme-multiplexing
#define MAX_N_SIZE HW_CELLS*TIME_MUX_FACTOR
#define MAX_TIME_MUX TIME_MUX_FACTOR
#define CONN_MATRIX_MAX MAX_N_SIZE*MAX_N_SIZE
#define CONN_MATRIX_HW  CONN_MATRIX_MAX/8
#define IAPP_MAX_CHARS 1 //2 integer, the dot, 2 decimals and the delimiter


// Cell properties///////////////////////////////////////////////////////////////////////////////////////////////
#define DELTA 0.05
//Conductance for neighbors' coupling
#define CONDUCTANCE 0.04
// Capacitance
#define C_M 1
// Somatic conductances (mS/cm2)
#define G_NA_S 150      // Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
#define G_KDR_S 9.0    // K delayed rectifier gate conductance (alternative value: 18)
#define G_K_S 5      // Voltage-dependent (fast) potassium
#define G_LS 0.016  // Leak conductance (0.015)
// Dendritic conductances (mS/cm2)
#define G_K_CA 35       // Potassium gate conductance (35)
#define G_CAH 4.5     // High-threshold Ca gate conductance (4.5)
#define G_LD 0.016   // Dendrite leak conductance (0.015)
#define G_H 0.125    // H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
// Axon hillock conductances (mS/cm2)
#define G_NA_A 240      // Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
#define G_NA_R 0      // Na (resurgent) gate conductance
#define G_K_A 20      // K voltage-dependent
#define G_LA 0.016  // Leak conductance
// Cell morphology
#define P1 0.25        // Cell surface ratio soma/dendrite (0.2)
#define P2 0.15      // Cell surface ratio axon(hillock)/soma (0.1)
#define G_INT 0.13       // Cell internal conductance (0.13)
// Reversal potentials
#define V_NA 55       // Na reversal potential (55)
#define V_K -75       // K reversal potential
#define V_CA 120       // Ca reversal potential (120)
#define V_H -43       // H current reversal potential
#define V_L 10       // leak current




//typedef double mod_prec;
typedef float mod_prec;

typedef struct returnState{
	mod_prec axonOut[IO_NETWORK_SIZE];
	//mod_prec dendOut[IO_NETWORK_SIZE];
	//mod_prec somaOut[IO_NETWORK_SIZE];
}returnState;

//Structs for internal calculations and general inputs
typedef struct Dend{
	mod_prec V_dend;
	mod_prec Hcurrent_q;
	mod_prec Calcium_r;
	mod_prec Potassium_s;
	mod_prec I_CaH;
	mod_prec Ca2Plus;
}Dend;

typedef struct Soma{
	mod_prec g_CaL;
	mod_prec V_soma;
	mod_prec Sodium_m;
	mod_prec Sodium_h;
	mod_prec Calcium_k;
	mod_prec Calcium_l;
	mod_prec Potassium_n;
	mod_prec Potassium_p;
	mod_prec Potassium_x_s;
}Soma;

typedef struct Axon{
	mod_prec V_axon;
	mod_prec Sodium_m_a;
	mod_prec Sodium_h_a;
	mod_prec Potassium_x_a;
}Axon;

typedef struct cellState{
	 Dend dend;
	 Soma soma;
	 Axon axon;
}cellState;

typedef struct cellCompParams{
	mod_prec iAppIn;
	mod_prec neighVdend[IO_NETWORK_SIZE];
	cellState newCellState;
}cellCompParams;

typedef struct channelParams{
	mod_prec v;
	mod_prec prevComp1, prevComp2;
	mod_prec newComp1, newComp2;
}channelParams;

typedef struct dendCurrVoltPrms{
	mod_prec iApp;
	mod_prec iC;
	mod_prec vDend;
	mod_prec vSoma;
	mod_prec q, r, s;
	mod_prec newVDend;
	mod_prec newI_CaH;
}dendCurrVoltPrms;

typedef struct somaCurrVoltPrms{
	mod_prec g_CaL;
	mod_prec vSoma;
	mod_prec vDend;
	mod_prec vAxon;
	mod_prec k, l, m, h, n, x_s;
	mod_prec newVSoma;
}somaCurrVoltPrms;

typedef struct axonCurrVoltPrms{
	mod_prec vSoma;
	mod_prec vAxon;
	mod_prec m_a, h_a, x_a;
	mod_prec newVAxon;
}axonCurrVoltPrms;

/*** FUNCTION PROTOTYPES ***/
//void ComputeNetwork(bool,bool ,cellState* , mod_prec*, int , int, mod_prec*,int,mod_prec* );
void ComputeNetwork(cellState* , mod_prec*, int,mod_prec*,mod_prec* );
void ComputeOneCellTimeMux0(int,  mod_prec*,mod_prec* , int, int, mod_prec*);


//void ComputeOneCell( int,int, mod_prec, mod_prec* );
void ComputeOneCell0( int,int, mod_prec ,mod_prec*, int, mod_prec*);

Dend CompDend(struct Dend, mod_prec , mod_prec,mod_prec*, int , mod_prec*,int);
mod_prec DendHCurr(struct channelParams, mod_prec );
mod_prec DendCaCurr(struct channelParams, mod_prec );
mod_prec DendKCurr(struct channelParams );
mod_prec DendCal(struct channelParams );
dendCurrVoltPrms DendCurrVolt(struct dendCurrVoltPrms );
mod_prec IcNeighbors(mod_prec*, mod_prec, int, mod_prec*,int);

Soma CompSoma(struct Soma, mod_prec ,mod_prec );
channelParams SomaCalcium(struct channelParams );
channelParams SomaSodium(struct channelParams );
channelParams SomaPotassium(struct channelParams );
mod_prec SomaPotassiumX(struct channelParams );
mod_prec SomaCurrVolt(struct somaCurrVoltPrms );

Axon CompAxon(struct Axon , mod_prec );
channelParams AxonSodium(channelParams );
mod_prec AxonPotassium(channelParams );
mod_prec AxonCurrVolt(axonCurrVoltPrms );



//mod_prec Open_exp(mod_prec );
cellState InitState();
int ReadFileLine(char *, int, FILE *, mod_prec *);
//inline mod_prec min(mod_prec a, mod_prec b);
bool Check_bit(int , int );



#endif /* MAIN_H_ */
