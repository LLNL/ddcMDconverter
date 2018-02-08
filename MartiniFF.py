beads = {  # Name: Mass
    "P5": 72, "P4": 72, "P3": 72, "P2": 72, "P1": 72,  # Polar
    "Nda": 72, "Nd": 72, "Na": 72, "N0": 72,  # Intermediate
    "C5": 72, "C4": 72, "C3": 72, "C2": 72, "C1": 72,  # Apolar
    "Qda": 72, "Qd": 72, "Qa": 72, "Q0": 72,  # Charged
    "SP5": 45, "SP4": 45, "SP3": 45, "SP2": 45, "SP1": 45,  # Ring type, polar
    "SNda": 45, "SNd": 45, "SNa": 45, "SN0": 45,  # Ring type, intermediate
    "SC5": 45, "SC4": 45, "SC3": 45, "SC2": 45, "SC1": 45,  # Ring type, apolar
    "SQda": 45, "SQd": 45, "SQa": 45, "SQ0": 45,  # Ring type, charged
    "AC1": 72, "AC2": 72,  # Amino acid side chains (Q-C interactions)
}


epsilon = {
    "0": 0.00,  # no interaction
    "A": 5.60,  # supra attractive
    "B": 5.00,  # attractive
    "C": 4.50,  # almost attractive
    "D": 4.00,  # semi attractive
    "E": 3.50,  # intermediate
    "F": 3.10,  # almost intermediate
    "G": 2.70,  # semi repulsive
    "H": 2.30,  # almost repulsive
    "I": 2.00,  # repulsive
    "Z": 0.25,  # with sigma=1 yields c6=c12=1
}

sigma = {
    "0": 0.00,  # no interaction
    "1": 0.43,  # ring bead type
    "2": 0.47,  # default bead type
    "3": 0.57,  # supra attractive bead type
    "4": 0.62,  # super repulsive bead type
    "5": 1.00,  # unity: c6/12 parity
    "6": 0.32,  # for tiny atom type
}

scale = {
    "a": 1.00,  # no scaling
    "b": 0.95,  # 95%: interaction with polarizable water
    "c": 0.90,  # 90%: S* - C1 interactions with BMW water
    "d": 0.75,  # 75%: ring bead types (S*) and uncharged interactions with BMW for epsilon >=4.5
    "e": 0.71,  # 71%: uncharged interactions with BMW for epsilon<4.5
    "f": 0.68,  # 68% for tiny atom type
}

# BEAD TYPES

# Default bead types; 4:1 mapping, 72 AMU
plain = {
    "P5":   72,  "P4":   72,  "P3":   72,  "P2":   72,  "P1":   72, # Polar
    "Nda":  72,  "Nd":   72,  "Na":   72,  "N0":   72,              # Intermediate
    "C5":   72,  "C4":   72,  "C3":   72,  "C2":   72,  "C1":   72, # Apolar
    "Qda":  72,  "Qd":   72,  "Qa":   72,  "Q0":   72,              # Charged
    "AC1":  72,  "AC2":  72,                                        # Amino acid specific (Q-C interactions)
    }

# Ring bead types; mapping 2:1 or 3:1, 45 AMU
small = {
    "SP5":  45,  "SP4":  45,  "SP3":  45,  "SP2":  45,  "SP1":  45, # Ring type, polar
    "SNda": 45,  "SNd":  45,  "SNa":  45,  "SN0":  45,              # Ring type, intermediate
    "SC5":  45,  "SC4":  45,  "SC3":  45,  "SC2":  45,  "SC1":  45, # Ring type, apolar
    "SQda": 45,  "SQd":  45,  "SQa":  45,  "SQ0":  45,              # Ring type, charged
    }

tiny = {
    "TP5":  45,  "TP4":  45,  "TP3":  45,  "TP2":  45,  "TP1":  45, # Ring type, polar
    "TNda": 45,  "TNd":  45,  "TNa":  45,  "TN0":  45,              # Ring type, intermediate
    "TC5":  45,  "TC4":  45,  "TG3":  45,  "TG2":  45,  "TC1":  45, # modified
    "TQda": 45,  "TQd":  45,  "TQa":  45,  "TQ0":  45,              # Ring type, charged
    }

# Virtual sites, plain type; mapping 4:1, no mass
vsite = {
    "vP5":   0,  "vP4":   0,  "vP3":   0,  "vP2":   0,  "vP1":   0, # Polar
    "vNda":  0,  "vNd":   0,  "vNa":   0,  "vN0":   0,              # Intermediate
    "vC5":   0,  "vC4":   0,  "vC3":   0,  "vC2":   0,  "vC1":   0, # Apolar
    "vQda":  0,  "vQd":   0,  "vQa":   0,  "vQ0":   0,              # Charged
    "vAC1":  0,  "vAC2":  0,                                        # Amino acid specific (Q-C interactions)
    }

# Virtual sites, small type; mapping 2:1 or 3:1, no mass
svste = {
    "vSP5":  0,  "vSP4":  0,  "vSP3":  0,  "vSP2":  0,  "vSP1":  0, # Ring type, polar
    "vSNda": 0,  "vSNd":  0,  "vSNa":  0,  "vSN0":  0,              # Ring type, intermediate
    "vSC5":  0,  "vSC4":  0,  "vSC3":  0,  "vSC2":  0,  "vSC1":  0, # Ring type, apolar
    "vSQda": 0,  "vSQd":  0,  "vSQa":  0,  "vSQ0":  0,              # Ring type, charged
    }

# Other: Special purpose types
other = {
    "BP4":  72, # Big water particle (antifreeze).
    "D":     0, # Dummy particle type
    }

classes = ("plain","small","tiny","vsite","svste","other")

table_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1  AC2  AC1
  Qda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
   Qd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
   Qa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
   Q0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
   P5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
   P4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
   P3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2  Ga2  Ha2
   P2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2  Ga2  Ha2
   P1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
  Nda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
   Nd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
   Na  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
   N0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
   C5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
   C4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
   C3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
   C2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
   C1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
  AC2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
  AC1  Ia2  Ia2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
"""

table_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
 SQda  Ad1  Ad1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  SQd  Ad1  Bd1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Dd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  SQa  Ad1  Ad1  Bd1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  SQ0  Cd1  Cd1  Cd1  Ed1  Bd1  Ad1  Bd1  Cd1  Dd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  SP5  Ad1  Ad1  Ad1  Bd1  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
  SP4  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Cd1  Cd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
  SP3  Ad1  Ad1  Ad1  Bd1  Ad1  Bd1  Bd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Fd1  Gd1  Hd1
  SP2  Bd1  Bd1  Bd1  Cd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Hd1
  SP1  Bd1  Bd1  Bd1  Dd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Ed1  Fd1  Gd1
 SNda  Bd1  Bd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
  SNd  Bd1  Dd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Dd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
  SNa  Bd1  Bd1  Dd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
  SN0  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Dd1  Dd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Gd1
  SC5  Fd1  Fd1  Fd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
  SC4  Gd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Ed1  Ed1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
  SC3  Hd1  Hd1  Hd1  Hd1  Gd1  Gd1  Fd1  Fd1  Ed1  Gd1  Gd1  Gd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1
  SC2  Ia4  Ia4  Ia4  Ia4  Hd1  Hd1  Gd1  Gd1  Fd1  Gd1  Gd1  Gd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1
  SC1  Ia4  Ia4  Ia4  Ia4  Id1  Id1  Hd1  Hd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Fd1  Ed1  Ed1  Ed1
"""

table_tiny = """
      TQda  TQd  TQa  TQ0  TP5  TP4  TP3  TP2  TP1 TNda  TNd  TNa  TN0  TC5  TC4  TG3  TG2  TC1
 TQda  Af1  Af1  Af1  Cf1  Af1  Af1  Af1  Bf1  Bf1  Bf1  Bf1  Bf1  Ef1  Ff1  Gf1  Hf1  If4  If4
  TQd  Af1  Bf1  Af1  Cf1  Af1  Af1  Af1  Bf1  Bf1  Bf1  Df1  Bf1  Ef1  Ff1  Gf1  Hf1  If4  If4
  TQa  Af1  Af1  Bf1  Cf1  Af1  Af1  Af1  Bf1  Bf1  Bf1  Bf1  Df1  Ef1  Ff1  Gf1  Hf1  If4  If4
  TQ0  Cf1  Cf1  Cf1  Ef1  Bf1  Af1  Bf1  Cf1  Df1  Df1  Df1  Df1  Ef1  Ff1  Gf1  Hf1  If4  If4
  TP5  Af1  Af1  Af1  Bf1  Af1  Af1  Af1  Af1  Af1  Bf1  Bf1  Bf1  Ef1  Ff1  Gf1  Gf1  Hf1  If1
  TP4  Af1  Af1  Af1  Af1  Af1  Bf1  Bf1  Cf1  Cf1  Df1  Df1  Df1  Ef1  Ff1  Gf1  Gf1  Hf1  If1
  TP3  Af1  Af1  Af1  Bf1  Af1  Bf1  Bf1  Cf1  Cf1  Cf1  Cf1  Cf1  Ef1  Ef1  Ff1  Ff1  Gf1  Hf1
  TP2  Bf1  Bf1  Bf1  Cf1  Af1  Cf1  Cf1  Cf1  Cf1  Cf1  Cf1  Cf1  Df1  Ef1  Ef1  Ff1  Gf1  Hf1
  TP1  Bf1  Bf1  Bf1  Df1  Af1  Cf1  Cf1  Cf1  Cf1  Cf1  Cf1  Cf1  Df1  Ef1  Ef1  Ef1  Ff1  Gf1
 TNda  Bf1  Bf1  Bf1  Df1  Bf1  Df1  Cf1  Cf1  Cf1  Cf1  Cf1  Cf1  Ef1  Ef1  Ff1  Gf1  Gf1  Gf1
  TNd  Bf1  Df1  Bf1  Df1  Bf1  Df1  Cf1  Cf1  Cf1  Cf1  Df1  Cf1  Ef1  Ef1  Ff1  Gf1  Gf1  Gf1
  TNa  Bf1  Bf1  Df1  Df1  Bf1  Df1  Cf1  Cf1  Cf1  Cf1  Cf1  Da6  Ea6  Ef1  Ff1  Ca6  Ca6  Gf1
  TN0  Ef1  Ef1  Ef1  Ef1  Ef1  Ef1  Ef1  Df1  Df1  Ef1  Ef1  Ea6  Ea6  Ef1  Ef1  Da6  Da6  Gf1
  TC5  Ff1  Ff1  Ff1  Ff1  Ff1  Ff1  Ef1  Ef1  Ef1  Ef1  Ef1  Ef1  Ef1  Ef1  Ef1  Ef1  Ff1  Ff1
  TC4  Gf1  Gf1  Gf1  Gf1  Gf1  Gf1  Ff1  Ef1  Ef1  Ff1  Ff1  Ff1  Ef1  Ef1  Ef1  Ef1  Ff1  Ff1
  TG3  Hf1  Hf1  Hf1  Hf1  Gf1  Gf1  Ff1  Ff1  Ef1  Ca6  Gf1  Ca6  Da6  Ef1  Ef1  Ga6  Ga6  Ef1
  TG2  If4  If4  If4  If4  Hf1  Hf1  Gf1  Gf1  Ff1  Ca6  Gf1  Ca6  Da6  Ff1  Ff1  Ga6  Ha6  Ef1
  TC1  If4  If4  If4  If4  If1  If1  Hf1  Hf1  Gf1  Gf1  Gf1  Gf1  Gf1  Ff1  Ff1  Ef1  Ef1  Ef1
"""

table_small_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1  AC2  AC1
 SQda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  SQd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  SQa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  SQ0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  SP5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  SP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  SP3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2  Ga2  Ha2
  SP2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2  Ga2  Ha2
  SP1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
 SNda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  SNd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  SNa  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  SN0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
  SC5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  SC4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  SC3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
  SC2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
  SC1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
"""

table_tiny_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1  AC2  AC1
 TQda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  TQd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  TQa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  TQ0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  TP5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  TP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  TP3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2  Ga2  Ha2
  TP2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2  Ga2  Ha2
  TP1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
 TNda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  TNd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  TNa  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  TN0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
  TC5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  TC4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  TG3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
  TG2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
  TC1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
"""

table_tiny_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
 TQda  Ad1  Ad1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  TQd  Ad1  Bd1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Dd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  TQa  Ad1  Ad1  Bd1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  TQ0  Cd1  Cd1  Cd1  Ed1  Bd1  Ad1  Bd1  Cd1  Dd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
  TP5  Ad1  Ad1  Ad1  Bd1  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
  TP4  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Cd1  Cd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
  TP3  Ad1  Ad1  Ad1  Bd1  Ad1  Bd1  Bd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Fd1  Gd1  Hd1
  TP2  Bd1  Bd1  Bd1  Cd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Hd1
  TP1  Bd1  Bd1  Bd1  Dd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Ed1  Fd1  Gd1
 TNda  Bd1  Bd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
  TNd  Bd1  Dd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Dd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
  TNa  Bd1  Bd1  Dd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
  TN0  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Dd1  Dd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Gd1
  TC5  Fd1  Fd1  Fd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
  TC4  Gd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Ed1  Ed1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
  TG3  Hd1  Hd1  Hd1  Hd1  Gd1  Gd1  Fd1  Fd1  Ed1  Gd1  Gd1  Gd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1
  TG2  Ia4  Ia4  Ia4  Ia4  Hd1  Hd1  Gd1  Gd1  Fd1  Gd1  Gd1  Gd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1
  TC1  Ia4  Ia4  Ia4  Ia4  Id1  Id1  Hd1  Hd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Fd1  Ed1  Ed1  Ed1
"""

table_vsite = """
      vQda  vQd  vQa  vQ0  vP5  vP4  vP3  vP2  vP1 vNda  vNd  vNa  vN0  vC5  vC4  vC3  vC2  vC1 vAC2 vAC1
 vQda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vQd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vQa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vQ0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vP5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  vP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  vP3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2  Ga2  Ha2
  vP2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2  Ga2  Ha2
  vP1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
 vNda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  vNd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  vNa  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  vN0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
  vC5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  vC4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  vC3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
  vC2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
  vC1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
 vAC2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
 vAC1  Ia2  Ia2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
"""

table_vsite_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1  AC2  AC1
 vQda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vQd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vQa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vQ0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
  vP5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  vP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
  vP3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2  Ga2  Ha2
  vP2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2  Ga2  Ha2
  vP1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
 vNda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  vNd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  vNa  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
  vN0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
  vC5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  vC4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
  vC3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
  vC2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
  vC1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
 vAC2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
 vAC1  Ia2  Ia2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
"""

table_vsite_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
 vQda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4
  vQd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4
  vQa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4
  vQ0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4
  vP5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2
  vP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2
  vP3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2
  vP2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2
  vP1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2
 vNda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2
  vNd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2
  vNa  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2
  vN0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2
  vC5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2
  vC4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2
  vC3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
  vC2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2
  vC1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2
 vAC2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2
 vAC1  Ia2  Ia2  Ia2  Ia2  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2
"""

table_svste = """
      vSQda vSQd vSQa vSQ0 vSP5 vSP4 vSP3 vSP2 vSP1 vSNda vSNd vSNa vSN0 vSC5 vSC4 vSC3 vSC2 vSC1
vSQda  Ad1  Ad1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSQd  Ad1  Bd1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Dd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSQa  Ad1  Ad1  Bd1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSQ0  Cd1  Cd1  Cd1  Ed1  Bd1  Ad1  Bd1  Cd1  Dd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSP5  Ad1  Ad1  Ad1  Bd1  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
 vSP4  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Cd1  Cd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
 vSP3  Ad1  Ad1  Ad1  Bd1  Ad1  Bd1  Bd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Fd1  Gd1  Hd1
 vSP2  Bd1  Bd1  Bd1  Cd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Hd1
 vSP1  Bd1  Bd1  Bd1  Dd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Ed1  Fd1  Gd1
vSNda  Bd1  Bd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
 vSNd  Bd1  Dd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Dd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
 vSNa  Bd1  Bd1  Dd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
 vSN0  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Dd1  Dd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Gd1
 vSC5  Fd1  Fd1  Fd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
 vSC4  Gd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Ed1  Ed1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
 vSC3  Hd1  Hd1  Hd1  Hd1  Gd1  Gd1  Fd1  Fd1  Ed1  Gd1  Gd1  Gd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1
 vSC2  Ia4  Ia4  Ia4  Ia4  Hd1  Hd1  Gd1  Gd1  Fd1  Gd1  Gd1  Gd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1
 vSC1  Ia4  Ia4  Ia4  Ia4  Id1  Id1  Hd1  Hd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Fd1  Ed1  Ed1  Ed1
"""

table_svste_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1  AC2  AC1
vSQda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSQd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSQa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSQ0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSP5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
 vSP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
 vSP3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2  Ga2  Ha2
 vSP2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2  Ga2  Ha2
 vSP1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
vSNda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
 vSNd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
 vSNa  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
 vSN0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
 vSC5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
 vSC4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
 vSC3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
 vSC2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
 vSC1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
"""

table_svste_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
vSQda  Ad1  Ad1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSQd  Ad1  Bd1  Ad1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Dd1  Bd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSQa  Ad1  Ad1  Bd1  Cd1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Bd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSQ0  Cd1  Cd1  Cd1  Ed1  Bd1  Ad1  Bd1  Cd1  Dd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Hd1  Ia4  Ia4
 vSP5  Ad1  Ad1  Ad1  Bd1  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Bd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
 vSP4  Ad1  Ad1  Ad1  Ad1  Ad1  Bd1  Bd1  Cd1  Cd1  Dd1  Dd1  Dd1  Ed1  Fd1  Gd1  Gd1  Hd1  Id1
 vSP3  Ad1  Ad1  Ad1  Bd1  Ad1  Bd1  Bd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Fd1  Gd1  Hd1
 vSP2  Bd1  Bd1  Bd1  Cd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Hd1
 vSP1  Bd1  Bd1  Bd1  Dd1  Ad1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Ed1  Fd1  Gd1
vSNda  Bd1  Bd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
 vSNd  Bd1  Dd1  Bd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Dd1  Cd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
 vSNa  Bd1  Bd1  Dd1  Dd1  Bd1  Dd1  Cd1  Cd1  Cd1  Cd1  Cd1  Dd1  Ed1  Ed1  Fd1  Gd1  Gd1  Gd1
 vSN0  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Dd1  Dd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Gd1
 vSC5  Fd1  Fd1  Fd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
 vSC4  Gd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Ed1  Ed1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1  Ed1  Fd1  Fd1
 vSC3  Hd1  Hd1  Hd1  Hd1  Gd1  Gd1  Fd1  Fd1  Ed1  Gd1  Gd1  Gd1  Ed1  Ed1  Ed1  Ed1  Ed1  Ed1
 vSC2  Ia4  Ia4  Ia4  Ia4  Hd1  Hd1  Gd1  Gd1  Fd1  Gd1  Gd1  Gd1  Fd1  Fd1  Fd1  Ed1  Ed1  Ed1
 vSC1  Ia4  Ia4  Ia4  Ia4  Id1  Id1  Hd1  Hd1  Gd1  Gd1  Gd1  Gd1  Gd1  Fd1  Fd1  Ed1  Ed1  Ed1
"""

table_svste_vsite = """
      vQda  vQd  vQa  vQ0  vP5  vP4  vP3  vP2  vP1 vNda  vNd  vNa  vN0  vC5  vC4  vC3  vC2  vC1 vAC2 vAC1
vSQda  Aa2  Aa2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSQd  Aa2  Ba2  Aa2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Da2  Ba2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSQa  Aa2  Aa2  Ba2  Ca2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ba2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSQ0  Ca2  Ca2  Ca2  Ea2  Ba2  Aa2  Ba2  Ca2  Da2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ha2  Ia4  Ia4  Ia2  Ia2
 vSP5  Aa2  Aa2  Aa2  Ba2  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ba2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
 vSP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
 vSP3  Aa2  Aa2  Aa2  Ba2  Aa2  Ba2  Ba2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Fa2  Ga2  Ha2  Ga2  Ha2
 vSP2  Ba2  Ba2  Ba2  Ca2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ha2  Ga2  Ha2
 vSP1  Ba2  Ba2  Ba2  Da2  Aa2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
vSNda  Ba2  Ba2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
 vSNd  Ba2  Da2  Ba2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Da2  Ca2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
 vSNa  Ba2  Ba2  Da2  Da2  Ba2  Da2  Ca2  Ca2  Ca2  Ca2  Ca2  Da2  Ea2  Ea2  Fa2  Ga2  Ga2  Ga2  Ga2  Ga2
 vSN0  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Da2  Da2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Ga2  Fa2  Ga2
 vSC5  Fa2  Fa2  Fa2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
 vSC4  Ga2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Ea2  Ea2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Fa2  Fa2  Fa2  Fa2
 vSC3  Ha2  Ha2  Ha2  Ha2  Ga2  Ga2  Fa2  Fa2  Ea2  Ga2  Ga2  Ga2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2  Ea2
 vSC2  Ia4  Ia4  Ia4  Ia4  Ha2  Ha2  Ga2  Ga2  Fa2  Ga2  Ga2  Ga2  Fa2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
 vSC1  Ia4  Ia4  Ia4  Ia4  Ia2  Ia2  Ha2  Ha2  Ga2  Ga2  Ga2  Ga2  Ga2  Fa2  Fa2  Ea2  Ea2  Ea2  Ea2  Ea2
"""

table_other = """
       BP4
  BP4  Ba2
"""

table_other_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1  AC2  AC1
  BP4  Aa2  Aa2  Aa2  Aa2  Aa2  Aa3  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
"""

table_other_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
  BP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2
"""

table_other_vsite = """
      vQda  vQd  vQa  vQ0  vP5  vP4  vP3  vP2  vP1 vNda  vNd  vNa  vN0  vC5  vC4  vC3  vC2  vC1 vAC2 vAC1
  BP4  Aa2  Aa2  Aa2  Aa2  Aa2  Aa3  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2  Ha2  Ia2
"""

table_other_svste = """
      vSQda vSQd vSQa vSQ0 vSP5 vSP4 vSP3 vSP2 vSP1 vSNda vSNd vSNa vSN0 vSC5 vSC4 vSC3 vSC2 vSC1
  BP4  Aa2  Aa2  Aa2  Aa2  Aa2  Ba2  Ba2  Ca2  Ca2  Da2  Da2  Da2  Ea2  Fa2  Ga2  Ga2  Ha2  Ia2
"""

class MartiniFF:
    def __init__(self):
        # Gather all atom types
        # all, mass = zip(*[i for j in MartiniFF.classes for i in eval(j).items()])
        # virtual = MartiniFF.vsite.keys() + MartiniFF.svste.keys()

        # rla = range(len(all))
        # cmb = [(all[i], all[j]) for i in rla for j in rla[i:]]

        tables = []
        for i in classes:
            tables.append(globals().get("table_%s" % i))
            for j in classes:
                tables.append(globals().get("table_%s_%s" % (i, j)))

        tables = [[j.split() for j in i.split("\n") if j not in ("", "\n")] for i in tables if i]

        self.pairs = dict([j for i in tables for j in self.table2pairs(i)])

    def table2pairs(self, x):
        labels2 = x.pop(0)
        return [((i[0], j), k) for i in x for j, k in zip(labels2, i[1:])]

    def sigeps2c(self, eps=None, scl="a", sig="2"):
        """
        Convert a string encoding for epsilon, sigma and scaling to C6 and C12 parameters
        """
        if eps == None:
            return None, None
        if scl in sigma.keys():
            sig, scl = scl, "2"
        return (epsilon[eps] * scale[scl], sigma[sig])


    def getLJparm(self, atomTypeI, atomTypeJ):
        pairs=self.pairs
        return self.sigeps2c(*pairs.get((atomTypeI, atomTypeJ), pairs.get((atomTypeJ, atomTypeI), "")))