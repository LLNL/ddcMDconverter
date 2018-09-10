
martini_input="""
simulate SIMULATE
{
    type = MD;
    system=system;
   integrator=nglf;
   deltaloop=500000;
    maxloop =1000000;
   dt = 20;
   printrate=100;
   snapshotrate=0;
   checkpointrate=50000;
   nLoopDigits=12;
   gidFormat=hex;
    printinfo=printinfo;
    heap=heap;
    ddc = ddc;
   xanalysis = writeCharmm;
}
energyInfo ENERGYINFO{}
heap HEAP { size = 1000 ;}

//ddc DDC { }
//ddc DDC { lx=2;ly=2;lz=1; ddcRule = rule; loadbalance=balance;}
ddc DDC { lx=2;ly=2;lz=1; }

//balance LOADBALANCE {
//   type=bisection;
//}

rule DDCRULE { type = MARTINI;}

printinfo PRINTINFO
{
  PRESSURE=bar;
  VOLUME = Ang^3;
  TEMPERATURE = K;
  ENERGY = kJ/mol;
  TIME = ns;
  printStress=1;
}

martini  POTENTIAL
{
   type = MARTINI;
   excludePotentialTerm=0;
   use_transform=0;
   cutoff=11.0 Angstrom;
   rcoulomb=11.0 Angstrom; epsilon_r=15; epsilon_rf=-1;
   function=lennardjones;
   parmfile=martini.data;
}

//nglf INTEGRATOR {type = NGLF; }
nglf INTEGRATOR {type = NGLFCONSTRAINT; T=310K; P0 = 1.0 bar; beta = 3.0e-4/bar; tauBarostat = 1.0 ps;}

system SYSTEM
{
   type = NORMAL;
   potential = martini;
   neighbor=nbr;
   groups= group free ;
   random = lcg64;
   box = box;
   collection=collection;
   moleculeClass = moleculeClass;
}


thermal TRANSFORM { type=THERMALIZE; temperature=0.003eV/kB;}
vcm TRANSFORM { type=SETVELOCITY; vcm=0 0 0;}

box BOX { type=ORTHORHOMBIC; pbc=7; }

nbr NEIGHBOR { type = NORMAL; deltaR=4.0000;minBoxSide=6; }

group GROUP { type = LANGEVIN; Teq=310K; tau=1ps; useDefault=0;}
//group   GROUP { type = FREE;  }
free   GROUP { type = FREE;  }

lcg64 RANDOM {type = LCG64;randomizeSeed=1;}

martini ANALYSIS { type =  PAIRCORRELATION; eval_rate=100; delta_r =0.05; length=40; outputrate=1; }
vcm ANALYSIS { type =  vcmWrite; outputrate=1; }
writeCharmm  ANALYSIS { type = subsetWrite; outputrate=1000; format=binaryCharmm; }

"""