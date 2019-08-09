
martini_input="""
simulate SIMULATE
{
    type = MD;
    system=system;
   integrator=nglf;
   deltaloop=500;
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

ddc DDC { updateRate=20; }

printinfo PRINTINFO
{
  PRESSURE=bar;
  VOLUME = Ang^3;
  TEMPERATURE = K;
  ENERGY = kJ/mol;
  TIME = ns;
  printStress=0;
  printMolecularPressure=1;
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
   use_gpu=1;
}

restraint POTENTIAL { type=RESTRAINT; use_gpu=1; parmfile=restraint.data; }

nglf INTEGRATOR { type=NGLFCONSTRAINTGPULANGEVIN; T=310K; P0 = 1.0 bar; beta = 3.0e-4/bar; tauBarostat = 1.0 ps;}

system SYSTEM
{
   type = NORMAL;
   potential = martini restraint ;
   neighbor=nbr;
   groups= group free ;
   random = lcg64;
   box = box;
   collection=collection;
   moleculeClass = moleculeClass;
   nConstraints=5866;
}


box BOX { type=ORTHORHOMBIC; pbc=7; }

nbr NEIGHBOR { type = NORMAL; deltaR=4.0000; minBoxSide=6; }

group GROUP { type = LANGEVIN; Teq=310K; tau=1ps; useDefault=0;}
free GROUP { type = LANGEVIN; Teq=310K; tau=1ps; useDefault=0;}

lcg64 RANDOM {type = LCG64;randomizeSeed=1;}

//martini ANALYSIS { type =  PAIRCORRELATION; eval_rate=100; delta_r =0.05; length=40; outputrate=1; }
//vcm ANALYSIS { type =  vcmWrite; outputrate=1; }
writeCharmm  ANALYSIS { type = subsetWrite; outputrate=1000; format=binaryCharmm; }

"""