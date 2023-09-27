#!/usr/bin/env python

import sys,math,random,re
import phaseHeight

version   = "20130516.10.TAW"
previous  = "20130416.16.TAW"


# Modify insane to take in arbitary lipid definition strings and use them as a template for lipids
# Also take in lipid name
# Edits: by Helgi I. Ingolfsson (all edits are marked with: # HII edit - lipid definition )

# PROTOLIPID (diacylglycerol), 18 beads
#
# 1-3-4-6-7--9-10-11-12-13-14
#  \| |/  |
#   2 5   8-15-16-17-18-19-20
#

lipidsx = {}
lipidsy = {}
lipidsz = {}
lipidsa = {}
#
## Diacyl glycerols
moltype = "lipid"
lipidsx[moltype] = (    0, .5,  0,  0, .5,  0,  0, .5,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipidsy[moltype] = (    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipidsz[moltype] = (   10,  9,  9,  8,  8,  7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipidsa.update({      # 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
## Phospholipids
    "DPPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DHPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    "DLPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DMPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DSPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
# POPC with 4 bead oleoyl
    "POPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DOPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
# POPC with 5 bead oleoyl
    "POP5": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
    "DOP5": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "DAPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  D1B D2B D3B D4B C5B  - "),
    "OIPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    "DIPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - "),
    "PAPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  C1B C2B C3B C4B  -   - "),
    "DEPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A D4A C5A C6A C1B C2B C3B D4B C5B C6B"),
    "DPPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DHPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    "DLPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DMPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DSPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "POPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "POPY": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DOPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    "PPCS": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - "),
    "DOPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "POPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
    "DOPS": (moltype, " -   -   -  CN0  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "POPS": (moltype, " -   -   -  CNO  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "PAPS": (moltype, " -   -   -  CNO  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  C1B C2B C3B C4B  -   - "),
    "OAPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  C1B D2B C3B C4B  -   - "),
    "DIPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - "),
    "PAPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  C1B C2B C3B C4B  -   - "),
    "DAPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  D1B D2B D3B D4B C5B  - "),
    "DAPY": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  D1B D2B D3B D4B C5B  - "),
    "DPSM": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A  -   -   -  C1B C2B C3B C4B  -   - "),
    "PBSM": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A C4A  -   -  C1B C2B C3B C4B C5B   - "),
    "DBSM": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A  -   -   -  C1B C2B C3B C4B C5B   - "),
    "DPSX": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A  -   -   -  C1B C2B C3B C4B  -   - "),
    "POPX": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DOPX": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
# PG for thylakoid membrane of T. vulcanus
     "CPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B  -   - "),
# PG for thylakoid membrane of spinach (PPT with a trans-unsaturated bond at sn1 and a triple-unsaturated bond at sn2,
# and PPG  with a transunsaturated bond at sn1 and a palmitoyl tail at sn2)
     "PPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  - "),
     "PPT": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A D3A D4A  -   -  D1B C2B C3B C4B  - "),
## Glycolipids
    "DSMG": (moltype, " -   -   -  C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "DSDG": (moltype, "C61 C41 C11 C62 C42 C12 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "DSSQ": (moltype, " -   -   S6 C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
     "CYF": (moltype, " -   -   -   -   -   -  F1   -  F2  F3  F4   -   -   -   -   -   -   -   - "),
     "FAR": (moltype, " -   -   -   -   -   -  C1   -  C2  C3   -   -   -   -   -   -   -   -   - "),
})

# HII fix for PI templates and new templates PI(s) with diffrent tails, PO-PIP1(3) and POPIP2(4,5)
#Prototopology for phosphatidylinositol type lipids 5,6,7 are potentail phosphates (PIP1,PIP2 and PIP3)
# 1,2,3 - is the inositol and 4 is the phosphate that links to the tail part.
#  5
#   \
#  6-2-1-4-8--10-11-12-13-14-15
#    |/    |
#  7-3     9--16-17-18-19-20-21
moltype = "INOSITOLLIPIDS"
lipidsx[moltype] = (   .5,  .5,   0,   0,   1, .5,  0,  0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipidsy[moltype] = (    0,   0,   0,   0,   0,  0,  0,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipidsz[moltype] = (    8,   9,   9,   7,  10, 10, 10,  6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipidsa.update({      # 1     2    3    4    5   6   7   8    9    10    11    12    13    14   15    16    17    18    19   20
    "DPPI": (moltype, " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "DOPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  D3A  C4A  C5A   -   C1B  C2B  D3B  C4B  C5B   - "),
    "PI"  : (moltype, " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
    "PI34": (moltype, " C1   C2   C3    CP PO1 PO2   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
   "POPI" : (moltype, " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
   "PIPI" : (moltype, " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  D2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
   "PAPI" : (moltype, " C1   C2   C3    CP   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
   "PUPI" : (moltype, " C1   C2   C3    CP   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  D5A  C6A  C1B  C2B  C3B  C4B   -    - "),
   "POP1" : (moltype, " C1   C2   C3    CP  P1   -   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
   "POP2" : (moltype, " C1   C2   C3    CP  P1  P2   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
   "POP3" : (moltype, " C1   C2   C3    CP  P1  P2  P3  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
   "PAP6" : (moltype, " C1   C2   C3   PO4  P4  P5   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
})


#Prototopology for longer and branched glycosil and ceramide based glycolipids
#
#     17-15-14-16
#         |/
#        13
#         |
# 12-10-9-7-6-4-3-1--18--26-27-28-29-30-31
#  |/   |/  |/  |/    |
#  11   8   5   2    19--20-21-22-23-24-25
moltype = "GLYCOLIPIDS"
lipidsx[moltype] = (    0,  .5,   0,   0,  .5,  0,  0, .5,  0,    0,   .5,    0,    0,    0,   0,    0,    0,    0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipidsy[moltype] = (    0,   0,   0,   0,   0,  0,  0,  0,  0,    0,    0,    0,   .5,    1,   1,    1,    1,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipidsz[moltype] = (    7,   8,   8,   9,  10, 10, 11, 12, 12,   13,   14,   14,   11,   10,  11,    9,   12,    6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipidsa.update({      # 1     2    3    4    5   6   7   8   9    10    11    12    13    14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31
    "GM1" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17  GM18  GM19 GM20 GM21 GM22 GM23 GM24   -  GM25 GM26 GM27 GM28   -    - "),
    "DGDG": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "MGDG": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "SQDG": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "CER" : (moltype, "  -    -    -    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "GCER": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
#lipids for thylakoid membrane of cyanobacteria: oleoyl tail at sn1 and palmiotyl chain at sn2. SQDG no double bonds
   "CDGDG": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  D3B  C4B   -    - "),
   "CMGDG": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  D3B  C4B   -    - "),
   "CSQDG": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
   "CSQDB": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  D3B  C4B   -    - "),
#lipids for thylakoid membrane of spinach: for the *T both chains are triple unsaturated and the *G have a triple unsaturated chain at sn1 and a pamitoyl chain at sn2.
   "PDGDG": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  D3B  D4B   -    - "),
   "PDGDT": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  D2A  D3A  D4A   -    -   C1B  D2B  D3B  D4B   -    - "),
   "PMGDG": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  D3B  D4B   -    - "),
   "PMGDT": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  D2A  D3A  D4A   -    -   C1B  D2B  D3B  D4B   -    - "),
   "PSQDG": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  D3B  D4B   -    - "),
# HII add GM1(s) with dififrent tails  -  2013.11.10 + extended matrix by x2 columns
   "DPG1" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
   "DXG1" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
   "PNG1" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
   "XNG1" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
   "DPG3" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
   "DXG3" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
   "PNG3" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
   "XNG3" : (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
})


# Prototopology for mycolic acid(s)
#
#  1--2--3--4--5--6--7--8
#                       |
# 16-15-14-13-12-11-10--9
# |
# 17-18-19-20-21-22-23-24
#                     /
# 32-31-30-29-28-27-25-26
#

moltype = "MYCOLIC ACIDS"
lipidsx[moltype] = (      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    1,    1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1)
lipidsy[moltype] = (      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,  0,   0,   0)
lipidsz[moltype] = (      7,   6,   5,   4,   3,   2,   1,   0,   0,   1,   2,   3,   4,   5,   6,    7,    7,    6,    5,   4,   3,   2,   1,   0,   1,   0,   2,   3,   4,  5,   6,   7)
lipidsa.update({        # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31   32
    "AMA":   (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    "AMA.w": (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    "KMA":   (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    "MMA":   (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
})


# Sterols
moltype = "sterol"
lipidsx[moltype] = (     0,  0,  0,  0,  0, 0,   0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0)
lipidsy[moltype] = (     0,  0,  0,  0,  0, 0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipidsz[moltype] = (     0,  0,  0,  0,  0, 0, 5.3,4.5,3.9,3.3, 3 ,2.6,1.4,  0,  0,  0,  0,  0)
lipidsa.update({
    "CHOL": (moltype, " -   -   -   -   -   -  ROH  R1  R2  R3  R4  R5  C1  C2  -   -   -   - "),
})

# Lists for automatic charge determination 
# 2020.06.11 Helgi - changed PAP6 from -5 to -4 
charges = {"ARG":1, "LYS":1, "ASP":-1, "GLU":-1, "DOPG":-1, "POPG":-1, "DOPS":-1, "PAPS":-1,  "POPS":-1, "PAP6":-4, "DSSQ":-1}

a,  b  = math.sqrt(2)/20, math.sqrt(2)/60
ct, st = math.cos(math.pi*109.47/180), math.sin(math.pi*109.47/180) # Tetrahedral

# Get a set of coordinates for a solvent particle with a given name
# Dictionary of solvents; First only those with multiple atoms
solventParticles = {
    "PW":       (("W",(-0.07,0,0)),                          # Polarizable water
                 ("WP",(0.07,0,0)),
                 ("WM",(0.07,0,0))),
    "BMW":      (("C",(0,0,0)),
                 ("Q1",(0.12,0,0)),
                 ("Q2",(-0.06,math.cos(math.pi/6)*0.12,0))), # BMW water
    "SPC":      (("OW",(0,0,0)),                             # SPC
                 ("HW1",(0.01,0,0)),
                 ("HW2",(0.01*ct,0.01*st,0))),
    "SPCM":     (("OW",(0,0,0)),                             # Multiscale/Martini SPC
                 ("HW1",(0.01,0,0)),
                 ("HW2",(0.01*ct,0.01*st,0)),
                 ("vW",(0,0,0))),
    "FG4W":     (("OW1",(a,a,a)),                            # Bundled water
                 ("HW11",(a,a-b,a-b)),
                 ("HW12",(a,a+b,a+b)),
                 ("OW2",(a,-a,-a)),
                 ("HW21",(a-b,-a,-a+b)),
                 ("HW22",(a+b,-a,-a-b)),
                 ("OW3",(-a,-a,a)),
                 ("HW31",(-a,-a+b,a-b)),
                 ("HW32",(-a,-a-b,a+b)),
                 ("OW4",(-a,a,-a)),
                 ("HW41",(-a+b,a,-a+b)),
                 ("HW42",(-a-b,a,-a-b))),
    "FG4W-MS":  (("OW1",(a,a,a)),                            # Bundled water, multiscaled
                 ("HW11",(a,a-b,a-b)),
                 ("HW12",(a,a+b,a+b)),
                 ("OW2",(a,-a,-a)),
                 ("HW21",(a-b,-a,-a+b)),
                 ("HW22",(a+b,-a,-a-b)),
                 ("OW3",(-a,-a,a)),
                 ("HW31",(-a,-a+b,a-b)),
                 ("HW32",(-a,-a-b,a+b)),
                 ("OW4",(-a,a,-a)),
                 ("HW41",(-a+b,a,-a+b)),
                 ("HW42",(-a-b,a,-a-b)),
                 ("VZ",(0,0,0))),
    "OCO":      (("PC",(-0.11, 0,   0)),
                 ("C" ,( 0.05, 0.16,0))),
    "GLUC":     (("B1",(-0.11, 0,   0)),
                 ("B2",( 0.05, 0.16,0)),
                 ("B3",( 0.05,-0.16,0))),
    "FRUC":     (("B1",(-0.11, 0,   0)),
                 ("B2",( 0.05, 0.16,0)),
                 ("B3",( 0.05,-0.16,0))),
    "SUCR":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "MALT":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "CELL":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "KOJI":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "SOPH":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "NIGE":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "LAMI":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    "TREH":     (("B1",(-0.25, 0.25,0)),
                 ("B2",(-0.25, 0,   0)),
                 ("B3",(-0.25,-0.25,0)),
                 ("B4",( 0.25, 0,   0)),
                 ("B5",( 0.25, 0.25,0)),
                 ("B6",( 0.25,-0.25,0))),
    }

# Update the solvents dictionary with single atom ones
for s in ["WN","W","NA","CL","Mg","K","BUT"]:
    solventParticles[s] = ((s,(0,0,0)),)

# Apolar amino acids nd stuff for orienting proteins in membrane
apolar = "ALA CYS PHE ILE LEU MET VAL TRP PLM CLR".split()

# 2020.06.19 Helgi - added code by Tomas Oppelstrup to fix asym rounding in sub-patches
def subsample(nsamples,propensities):
    ## Given a total number of lipids 'nsamples', and a list of propensities
    ## 'propensities' sample 'nsamples' lipids from the types to match the
    ## propensities as closely as possible, give the discretization

    proptot = sum(propensities)
    raw_counts = [ (nsamples+0.0)*p/proptot for p in propensities ]
    int_counts = [ int(x) for x in raw_counts ]
    int_sum = 0
    cdfsum = 0.0
    cdf = []
    for i in range(len(raw_counts)):
        cdfsum += raw_counts[i] - int_counts[i]
        cdf.append(cdfsum)
        int_sum += int_counts[i]

    if cdfsum == 0:
        return int_counts        
    cdf = [ x/cdfsum for x in cdf ]
    while int_sum < nsamples:
        j = 0
        r = random.random()
        while r > cdf[j]:
            j = j+1
        int_counts[j] += 1
        int_sum +=1

    return int_counts

## PRIVATE PARTS FROM THIS POINT ON ##

S = str
F = float
I = int
R = random.random

def vector(v):
    if type(v) == str and "," in v:
        return [float(i) for i in v.split(",")]
    return float(v)

def vvadd(a,b):
    if type(b) in (int,float):
        return [i+b for i in a]
    return [i+j for i,j in zip(a,b)]

def vvsub(a,b):
    if type(b) in (int,float):
        return [i-b for i in a]
    return [i-j for i,j in zip(a,b)]

def isPDBAtom(l):
    return l.startswith("ATOM") or l.startswith("HETATM")

def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    ## ===>   atom name,   res name,     res id, chain,       x,            y,             z
    return (S(a[12:16]),S(a[17:20]),I(a[22:26]),a[21],F(a[30:38])/10,F(a[38:46])/10,F(a[46:54])/10)

d2r = 3.14159265358979323846264338327950288/180
def pdbBoxRead(a):
    # Convert a PDB CRYST1 entry to a lattice definition.
    # Convert from Angstrom to nanometer
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    ca, cb, cg, sg         = math.cos(d2r*aa), math.cos(d2r*ab), math.cos(d2r*ac) , math.sin(d2r*ac)
    wx, wy                 = 0.1*fc*cb, 0.1*fc*(ca-cb*cg)/sg
    wz                     = math.sqrt(0.01*fc*fc - wx*wx - wy*wy)
    return [0.1*fa, 0, 0, 0.1*fb*cg, 0.1*fb*sg, 0, wx, wy, wz]

def groAtom(a):
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===>   atom name,   res name,     res id, chain,       x,          y,          z
    return (S(a[10:15]), S(a[5:10]),   I(a[:5]), " ", F(a[20:28]),F(a[28:36]),F(a[36:44]))

def groBoxRead(a):
    b = [F(i) for i in a.split()] + 6*[0] # Padding for rectangular boxes
    return b[0],b[3],b[4],b[5],b[1],b[6],b[7],b[8],b[2]

def readBox(a):
    x = [ float(i) for i in a.split(",") ] + 6*[0]
    if len(x) == 12: # PDB format
        return pdbBoxRead("CRYST1 "+" ".join([str(i) for i in x]))
    else:            # GRO format
        return x[0],x[3],x[4],x[5],x[1],x[6],x[7],x[8],x[2]

class Structure:
    def __init__(self,filename=None):
        self.title   = ""
        self.atoms   = []
        self.coord   = []
        self.rest    = []
        self.box     = []
        self._center = None

        if filename:
            lines = open(filename).readlines()
            # Try extracting PDB atom/hetatm definitions
            self.rest   = []
            self.atoms  = [pdbAtom(i) for i in lines if isPDBAtom(i) or self.rest.append(i)]
            if self.atoms:
                # This must be a PDB file
                self.title = "THIS IS INSANE!\n"
                for i in self.rest:
                    if i.startswith("TITLE"):
                        self.title = i
                self.box   = [0,0,0,0,0,0,0,0,0]
                for i in self.rest:
                    if i.startswith("CRYST1"):
                        self.box = pdbBoxRead(i)
            else:
                # This should be a GRO file
                self.atoms = [groAtom(i) for i in lines[2:-1]]
                self.rest  = [lines[0],lines[1],lines[-1]]
                self.box   = groBoxRead(lines[-1])
                self.title = lines[0]
            self.coord = [i[4:7] for i in self.atoms]
            self.center()


    def __bool__(self):
        return bool(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __iadd__(self,s):
        for i in range(len(self)):
            self.coord[i] = vvadd(self.coord[i],s)
        return self

    def center(self,other=None):
        if not self._center:
            self._center = [ sum(i)/len(i) for i in zip(*self.coord)]
        if other:
            s = vvsub(other,self._center)
            for i in range(len(self)):
                self.coord[i] = vvadd(self.coord[i],s)
            self._center = other
            return s # return the shift
        return self._center

    def diam(self):
        if self._center != (0,0,0):
            self.center((0,0,0))
        return 2*math.sqrt(max([i*i+j*j+k*k for i,j,k in self.coord]))

    def diamxy(self):
        if self._center != (0,0,0):
            self.center((0,0,0))
        return 2*math.sqrt(max([i*i+j*j for i,j,k in self.coord]))

    def fun(self,fn):
        return [fn(i) for i in zip(*self.coord)]

# Mean of deviations from initial value
def meand(v):
    return sum([i-v[0] for i in v])/len(v)

# Sum of squares/crossproducts of deviations
def ssd(u,v):
    return sum([(i-u[0])*(j-v[0]) for i,j in zip(u,v)])/(len(u)-1)

# Parse a string for a lipid as given on the command line (LIPID[:NUMBER])
def parse_mol(x):
    l = x.split(":")
    return l[0], len(l) == 1 and 1 or float(l[1])

def parse_patch(fname):
    fo = open(fname, "r")
    patch = fo.readline().strip("\n").split(" ")
    lip_list = []
    for tmp in fo:
      lip_list.append(tmp.strip("\n").split(" ") )
    fo.close()
    return int(patch[0]),int(patch[1]),lip_list

def parse_raslist(fname):
    fo = open(fname, "r")
    ras_list = []
    for tmp in fo:
      ras_list.append(tmp.strip("\n").split(" ") )
    fo.close()
    return ras_list

## MIJN EIGEN ROUTINE ##

# Quite short piece of code for diagonalizing symmetric 3x3 matrices :)

# Analytic solution for third order polynomial
def solve_p3( a, b, c ):
    Q,R,a3 = (3*b-a**2)/9.0, (-27*c+a*(9*b-2*a**2))/54.0, a/3.0
    if Q**3 + R**2:
        t,R13 = math.acos(R/math.sqrt(-Q**3))/3, 2*math.sqrt(-Q)
        u,v,w = math.cos(t), math.sin(t+math.pi/6), math.cos(t+math.pi/3)
        return R13*u-a3, -R13*v-a3, -R13*w-a3
    else:
        R13   = math.sqrt3(R)
        return 2*R13-a3, -R13-a3, -R13-a3

# Normalization of 3-vector
def normalize(a):
    f = 1.0/math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
    return f*a[0],f*a[1],f*a[2]

# Eigenvectors for a symmetric 3x3 matrix:
# For symmetric matrix A the eigenvector v with root r satisfies
#   v.Aw = Av.w = rv.w = v.rw
#   v.(A-rI)w = v.Aw - v.rw = 0 for all w
# This means that for any two vectors p,q the eigenvector v follows from:
#   (A-rI)p x (A-rI)q
# The input is var(x),var(y),var(z),cov(x,y),cov(x,z),cov(y,z)
# The routine has been checked and yields proper eigenvalues/-vectors
def mijn_eigen_sym_3x3(a,d,f,b,c,e):
    a,d,f,b,c,e=1,d/a,f/a,b/a,c/a,e/a
    b2, c2, e2, df = b*b, c*c, e*e, d*f
    roots = list(solve_p3(-a-d-f, df-b2-c2-e2+a*(f+d), a*e2+d*c2+f*b2-a*df-2*b*c*e))
    roots.sort(reverse=True)
    ux, uy, uz = b*e-c*d, b*c-a*e, a*d-b*b
    u = (ux+roots[0]*c,uy+roots[0]*e,uz+roots[0]*(roots[0]-a-d))
    v = (ux+roots[1]*c,uy+roots[1]*e,uz+roots[1]*(roots[1]-a-d))
    w = u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0] # Cross product
    return normalize(u),normalize(v),normalize(w),roots

# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __bool__(self):
        return self.value != None
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]

tm   = []
lipL = []
lipU = []
solv = []
rasl = []


# HII edit - lipid definition, for extra lipid definitaions
usrmols  = []
usrheads = []
usrlinks = []
usrtails = []
usrLipHeadMapp = { # Define supported lipid head beads. One letter name mapped to atom name
    "C":  ('NC3'), # NC3 = Choline
    "E":  ('NH3'), # NH3 = Ethanolamine
    "G":  ('GL0'), # GL0 = Glycerol
    "S":  ('CNO'), # CNO = Serine
    "P":  ('PO4'),  # PO4 = Phosphate
    "O":  ('PO4')  # PO4 = Phosphate acid
    }
usrIndexToLetter = "A B C D E F G H I J K L M N".split() # For naming lipid tail beads

# Description
desc = ""

# Option list
options = [
#   option           type number default description
# HII edit - lipid definition (last options are for additional lipid specification)
    """
Input/output related options
""",
    ("-f",      Option(tm.append,   1,        None, "Input GRO or PDB file 1: Protein")),
    ("-ras",    Option(rasl.append, 4,      None, "Input GRO or PDB file 4: Protein lu x y")),
    ("-rasf",   Option(str,         1,      None, "Input GRO or PDB file 1: Protein file")),
    ("-o",      Option(str,         1,        None, "Output GRO file: Membrane with Protein")),
    ("-p",      Option(str,         1,        None, "Optional rudimentary topology file")),
    """
Periodic boundary conditions
If -d is given, set up PBC according to -pbc such that no periodic
images are closer than the value given.  This will make the numbers
provided for lipids be interpreted as relative numbers. If -d is
omitted, those numbers are interpreted as absolute numbers, and the
PBC are set to fit the given number of lipids in.
""",
    ("-pbc",    Option(str,         1, "hexagonal", "PBC type: hexagonal, rectangular, square, cubic or optimal")),
    ("-d",      Option(float,       1,         2.5, "Distance between periodic images (nm)")),
    ("-dz",     Option(float,       1,           0, "Z distance between periodic images (nm)")),
    ("-x",      Option(vector,      1,           0, "X dimension or first lattice vector of system (nm)")),
    ("-y",      Option(vector,      1,           0, "Y dimension or first lattice vector of system (nm)")),
    ("-z",      Option(vector,      1,           0, "Z dimension or first lattice vector of system (nm)")),
    ("-box",    Option(readBox,     1,        None, "Box in GRO (3 or 9 floats) or PDB (6 floats) format, comma separated")),
    ("-n",      Option(str,         1,        None, "Index file --- TO BE IMPLEMENTED")),
    """
Membrane/lipid related options.
The options -l and -u can be given multiple times. Option -u can be
used to set the lipid type and abundance for the upper leaflet. Option
-l sets the type and abundance for the lower leaflet if option -u is
also given, or for both leaflets if option -u is not given. The
meaning of the number depends on whether option -d is used to set up
PBC
""",
    ("-l",      Option(lipL.append, 1,   None, "Lipid type and relative abundance (NAME[:#])")),
    ("-u",      Option(lipU.append, 1,   None, "Lipid type and relative abundance (NAME[:#])")),
    ("-lf",     Option(str,         1,   None, "Lipid type and relative abundance description file")),
    ("-uf",     Option(str,         1,   None, "Lipid type and relative abundance description file")),

    ("-a",      Option(float,       1,        0.60, "Area per lipid (nm*nm)")),
    ("-asym",   Option(int,         1,        None, "Membrane asymmetry (number of lipids)")),
    ("-hole",   Option(float,       1,        None, "Make a hole in the membrane with specified radius")),
    ("-rand",   Option(float,       1,         0.1, "Random kick size (maximum atom displacement)")),
    """
Protein related options.
""",
    ("-center", Option(bool,        0,        None, "Center the protein on z")),
    ("-orient", Option(bool,        0,        None, "Orient protein in membrane")),
    ("-rotate", Option(str,         0,        None, "Rotate protein (random|princ|angle(float))")),
    ("-od",     Option(float,       1,         1.0, "Grid spacing for determining orientation")),
    ("-op",     Option(float,       1,         4.0, "Hydrophobic ratio power for determining orientation")),
    ("-fudge",  Option(float,       1,         0.1, "Fudge factor for allowing lipid-protein overlap")),
    ("-ring",   Option(bool,        0,        None, "Put lipids inside the protein")),
    ("-dm",     Option(float,       1,        None, "Set distance between protein and membrane")),
    """
Solvent related options.
""",
    ("-sol",    Option(solv.append, 1,        None, "Solvent type and relative abundance (NAME[:#])")),
    ("-sold",   Option(float,       1,         0.5, "Solvent diameter")),
    ("-solr",   Option(float,       1,         0.1, "Solvent random kick")),
    ("-excl",   Option(float,       1,         1.5, "Exclusion range (nm) for solvent addition relative to membrane center")),
    """
Salt related options.
""",
    ("-salt",   Option(str,         1,        None, "Salt concentration")),
    ("-charge", Option(str,         1,      "auto", "Charge of system. Set to auto to infer from residue names")),
    """
Define additional lipid types (same format as in lipid-martini-itp-v01.py)
""",
    ("-alname",  Option(usrmols.append,         1,        None, "Additional lipid name, x4 letter")),
    ("-alhead",  Option(usrheads.append,        1,        None, "Additional lipid head specification string")),
    ("-allink",  Option(usrlinks.append,        1,        None, "Additional lipid linker specification string")),
    ("-altail",  Option(usrtails.append,        1,        None, "Additional lipid tail specification string")),
    ]

args = sys.argv[1:]

if '-h' in args or '--help' in args:
    print("\n",__file__)
    print(desc or "\nSomeone ought to write a description for this script...\n")
    for thing in options:
        print(type(thing) != str and "%10s  %s"%(thing[0],thing[1].description) or thing)
    print()
    sys.exit()


# Convert the option list to a dictionary, discarding all comments
options = dict([i for i in options if not type(i) == str])


# Process the command line
while args:
    ar = args.pop(0)
    options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])


absoluteNumbers = not options["-d"]

# HII edit - lipid definition
# Add specified lipid definition to insane lipid library
for name, head, link, tail in zip(usrmols,usrheads,usrlinks,usrtails):
    moltype = "usr_"+name
    lipidsx[moltype] = []
    lipidsy[moltype] = []
    lipidsz[moltype] = []
    headArray = (head).split()
    linkArray = (link).split()
    tailsArray = (tail).split()
    lipidDefString = ""

    if len(tailsArray) != len(linkArray):
        print("Error, Number of tails has to equal number of linkers")
        sys.exit()

    # Find longest tail
    maxTail = 0
    for cTail in tailsArray:
       if len(cTail) > maxTail:
           maxTail = len(cTail)
    cBeadZ = maxTail + len(headArray) # longest tail + linker (always x1) + lengths of all heads - 1 (as it starts on 0)

    # Add head beads
    for cHead in headArray:
        lipidsx[moltype].append(0)
        lipidsy[moltype].append(0)
        lipidsz[moltype].append(cBeadZ)
        cBeadZ -= 1
        lipidDefString += usrLipHeadMapp[cHead] + " "

    # Add linkers
    for i,cLinker in enumerate(linkArray):
        lipidsx[moltype].append(max(i-0.5,0))
        lipidsy[moltype].append(0)
        lipidsz[moltype].append(cBeadZ)
        # HII add 2013.10.11 AM linker support
        if cLinker == 'G':
            lipidDefString += "GL" + str(i+1) + " "
        elif cLinker == 'A':
            lipidDefString += "AM" + str(i+1) + " "
        else:
            print("Error, linker type not supported")
            sys.exit()

    # Add tails
    for i,cTail in enumerate(tailsArray):
        cBeadZ = maxTail - 1

        for j,cTailBead in enumerate(cTail):
            lipidsx[moltype].append(i)
            lipidsy[moltype].append(0)
            lipidsz[moltype].append(cBeadZ)
            cBeadZ -= 1
            lipidDefString += cTailBead + str(j+1) + usrIndexToLetter[i] + " "

    lipidsa[name] = (moltype,lipidDefString)
# End user lipid definition

# HII edit - lipid definition, had to move this one below the user lipid definitions to scale them to.
# First all X/Y coordinates of templates are centered and scaled (magic numbers!)
for i in list(lipidsx.keys()):
    cx = (min(lipidsx[i])+max(lipidsx[i]))/2
    lipidsx[i] = [0.25*(j-cx) for j in lipidsx[i]]
    cy = (min(lipidsy[i])+max(lipidsy[i]))/2
    lipidsy[i] = [0.25*(j-cy) for j in lipidsy[i]]

# Periodic boundary conditions

# option -box overrides everything
if options["-box"]:
    options["-x"].value = options["-box"].value[:3]
    options["-y"].value = options["-box"].value[3:6]
    options["-z"].value = options["-box"].value[6:]

# options -x, -y, -z take precedence over automatic determination
pbcSetX = 0
if type(options["-x"].value) in (list,tuple):
    pbcSetX = options["-x"].value
elif options["-x"].value:
    pbcSetX = [options["-x"].value,0,0]

pbcSetY = 0
if type(options["-y"].value) in (list,tuple):
    pbcSetY = options["-y"].value
elif options["-y"].value:
    pbcSetY = [0,options["-y"].value,0]

pbcSetZ = 0
if type(options["-z"].value) in (list,tuple):
    pbcSetZ = options["-z"].value
elif options["-z"].value:
    pbcSetZ = [0,0,options["-z"].value]

lipd  = math.sqrt(options["-a"].value)

tm    = [ Structure(i) for i in tm ]

if rasl:
  rasl  = [ [Structure(i),lu,x,y] for i,lu,x,y in zip(rasl[::4],rasl[1::4],rasl[2::4],rasl[3::4])]
elif options["-rasf"]:
  rasl = parse_raslist(option["-rasf"].value)


################
## I. PROTEIN ##
################


protein  = Structure()
prot     = []
shift    = [0] # Shift in x direction per protein

## A. NO PROTEIN ---
if not tm:

    resi = 0
    # Set the box -- If there is a hole, add its radius to the distance
    pbcx = pbcy = pbcz = options["-d"].value + (options["-hole"] and options["-hole"].value or 0)
    if "hexagonal".startswith(options["-pbc"].value):
        # Hexagonal prism -- y derived from x directly
        pbcy = math.sqrt(3)*pbcx/2
        pbcz = options["-dz"].value or options["-z"].value or options["-d"].value
    elif "optimal".startswith(options["-pbc"].value):
        # Rhombic dodecahedron with hexagonal XY plane
        pbcy = math.sqrt(3)*pbcx/2
        pbcz = math.sqrt(6)*options["-d"].value/3
    if "rectangular".startswith(options["-pbc"].value):
        pbcz = options["-dz"].value or options["-z"].value or options["-d"].value

    # Possibly override
    pbcx = pbcSetX and pbcSetX[0] or pbcx
    pbcy = pbcSetY and pbcSetY[1] or pbcy
    pbcz = pbcSetZ and pbcSetZ[2] or pbcz


## B. PROTEIN ---
else:

    for prot in tm:

        ## a. NO MEMBRANE --
        if ((not lipL) or (not options["-lf"])):

            # A protein, but don't add lipids... Just solvate the protein
            # Maybe align along principal axes and then build a cell according to PBC

            # Set PBC starting from diameter and adding distance
            if "cubic".startswith(options["-pbc"].value):
                pbcx = pbcy = pbcz = prot.diam()+options["-d"].value
            elif "rectangular".startswith(options["-pbc"].value):
                pbcx, pbcy, pbcz = vvadd(vvsub(prot.fun(max),prot.fun(min)),options["-d"].value)
            else:
                # Rhombic dodecahedron
                pbcx = pbcy = prot.diam()+options["-d"].value
                pbcz = math.sqrt(2)*pbcx/2

            # Possibly override
            pbcx = pbcSetX and pbcSetX[0] or pbcx
            pbcy = pbcSetY and pbcSetY[1] or pbcy
            pbcz = pbcSetZ and pbcSetZ[2] or pbcz

            # Center coordinates in rectangular brick -- Add solvent next
            if len(tm) == 1:
                prot.center((0.5*pbcx, 0.5*pbcy, 0.5*pbcz))

            # Do not set an exclusion range for solvent
            # HII remove this for now
            #options["-excl"].value = -1


        ## b. PROTEIN AND MEMBRANE --
        else:

            # Have to build a membrane around the protein.
            # So first put the protein in properly.


            # Center the protein and store the shift
            shift = prot.center((0,0,0))


            ## 1. Orient with respect to membrane
            # Orient the protein according to the TM region, if requested
            # This doesn't actually work very well...
            if options["-orient"]:

                # Grid spacing (nm)
                d  = options["-od"].value
                pw = options["-op"].value

                # Determine grid size
                mx,my,mz = prot.fun(min)
                rx,ry,rz = prot.fun(lambda x: max(x)-min(x)+1e-8)

                # Number of grid cells
                nx,ny,nz = int(rx/d+0.5),int(ry/d+0.5),int(rz/d+0.5)

                # Initialize grids
                atom     = [[[0 for i in range(nz+2)] for j in range(ny+2)] for k in range(nx+2)]
                phobic   = [[[0 for i in range(nz+2)] for j in range(ny+2)] for k in range(nx+2)]
                surface  = []
                for i, (ix, iy, iz) in zip(prot.atoms,prot.coord):
                    if i[1] != "DUM":
                        jx,jy,jz = int(nx*(ix-mx)/rx), int(ny*(iy-my)/ry), int(nz*(iz-mz)/rz)
                        atom[jx][jy][jz]   += 1
                        phobic[jx][jy][jz] += (i[1] in apolar)

                # Determine average density
                occupd = sum([bool(k) for i in atom for j in i for k in j])
                avdens = float(sum([sum(j) for i in atom for j in i]))/occupd

                #cgofile  = open('density.cgo',"w")
                #cgofile.write('[\n')
                for i in range(nx):
                    for j in range(ny):
                        for k in range(nz):
                            if atom[i][j][k] > 0.1*avdens:
                                # Check the neighbouring cells; If one of them is not occupied, count cell as surface
                                if not (atom[i-1][j][k] and atom[i+1][j][k] and
                                        atom[i][j-1][k] and atom[i][j+1][k] and
                                        atom[i][j][k-1] and atom[i][j][k+1]):
                                    sx,sy,sz = mx+rx*(i+0.5)/nx, my+ry*(j+0.5)/ny, mz+rz*(k+0.5)/nz
                                    sw       = (2.0*phobic[i][j][k]/atom[i][j][k])**pw
                                    surface.append((sx,sy,sz,sw))
                                    #cgofile.write("    7.0, %f, %f, %f, %f,\n"%(10*sx,10*sy,10*sz,0.25*sw))
                #cgofile.write(']\n')
                #cgofile.close()

                sx, sy, sz, w = list(zip(*surface))
                W             = 1.0/sum(w)

                # Weighted center of apolar region; has to go to (0,0,0)
                sxm,sym,szm   = [sum(p)*W for p in zip(*[(m*i,m*j,m*k) for m,i,j,k in zip(w,sx,sy,sz)])]

                # Place apolar center at origin
                prot.center((-sxm,-sym,-szm))
                sx, sy, sz    = list(zip(*[(i-sxm,j-sym,k-szm) for i,j,k in zip(sx,sy,sz)]))

                # Determine weighted deviations from centers
                dx,dy,dz      = list(zip(*[(m*i,m*j,m*k) for m,i,j,k in zip(w,sx,sy,sz)]))

                # Covariance matrix for surface
                xx,yy,zz,xy,yz,zx = [sum(p)*W for p in zip(*[(i*i,j*j,k*k,i*j,j*k,k*i) for i,j,k in zip(dx,dy,dz)])]

                # PCA: u,v,w are a rotation matrix
                (ux,uy,uz),(vx,vy,vz),(wx,wy,wz),r = mijn_eigen_sym_3x3(xx,yy,zz,xy,zx,yz)

                # Rotate the coordinates
                prot.coord = [(ux*i+uy*j+uz*k,vx*i+vy*j+vz*k,wx*i+wy*j+wz*k) for i,j,k in prot.coord]


            ## 4. Orient the protein in the xy-plane
            ## i. According to principal axes and unit cell
            if options["-rotate"].value == "princ":

                x, y, z = list(zip(*prot.coord))

                # The rotation matrix in the plane equals the transpose
                # of the matrix of eigenvectors from the 2x2 covariance
                # matrix of the positions.
                # For numerical stability we do
                # d_i     = x_i - x_0
                # mean(x) = x_0 + sum(d_i)/N =
                # var(x)  = sum((d_i - mean(d))**2)/(N-1)
                xy        = ssd(x,y)
                if xy != 0:
                    xx     = ssd(x,x)
                    yy     = ssd(y,y)

                    # The eigenvalues are the roots of the 2nd order
                    # characteristic polynomial, with the coefficients
                    # equal to the trace and the determinant of the
                    # matrix.
                    t,  d  = xx+yy, xx*yy - xy*xy
                    # The two eigenvectors form a 2D rotation matrix
                    # R = ((cos,sin),(-sin,cos)), which means that
                    # the second eigenvector follows directly from
                    # the first. We thus only need to determine one.
                    l1     = t/2 + math.sqrt(0.25*t*t-d)

                    ux, uy = l1-yy, xy
                    lu     = math.sqrt(ux*ux+uy*uy)

                    ux    /=  lu
                    uy    /=  lu

                    # Finally we rotate the system in the plane by
                    # matrix multiplication with the transpose of
                    # the matrix of eigenvectors
                    prot.coord = [(ux*i+uy*j,ux*j-uy*i,k) for i,j,k in zip(x,y,z)]

            ## ii. Randomly
            elif options["-rotate"].value == "random":
                ux   = math.cos(R()*2*math.pi)
                uy   = math.sqrt(1-ux*ux)
                prot.coord = [(ux*i+uy*j,ux*j-uy*i,k) for i,j,k in prot.coord]

            ## iii. Specifically
            elif options["-rotate"]:
                ux   = math.cos(float(options["-rotate"].value)*math.pi/180.)
                uy   = math.sin(float(options["-rotate"].value)*math.pi/180.)
                prot.coord = [(ux*i+uy*j,ux*j-uy*i,k) for i,j,k in prot.coord]



            ## 5. Determine the minimum and maximum x and y of the protein
            pmin, pmax = prot.fun(min), prot.fun(max)
            prng       = (pmax[0]-pmin[0],pmax[1]-pmin[1],pmax[2]-pmin[2])
            center     = (0.5*(pmin[0]+pmax[0]),0.5*(pmin[1]+pmax[1]))

            # Set the z-dimension
            pbcz = pbcSetZ and pbcSetZ[2] or min( options["-dz"].value+prng[2], options["-z"].value )

            # BLABLABLA...

            # At this point we should shift the subsequent proteins such that they end up
            # at the specified distance, in case we have a number of them to do
            # y-shift is always -ycenter
            # x-shift is -xmin+distance+xmax(current)
            xshft, yshft = shift[-1]-pmin[0]+options["-d"].value, -center[1]
            shift.append(shift[-1]+pmax[0]+options["-d"].value)

            ## 6. Set box (brick) dimensions
            pbcx = options["-d"].value + prng[0]
            if "square".startswith(options["-pbc"].value):
                pbcy = pbcx
            elif "rectangular".startswith(options["-pbc"].value):
                pbcy = options["-d"].value + prng[1]
            else:
                # This goes for a hexagonal cell as well as for the optimal arrangement
                # The latter is hexagonal in the membrane plane anyway...
                pbcy  = math.cos(math.pi/6)*pbcx

            ## 6. Adjust PBC for hole
            # If we need to add a hole, we have to scale the system
            # The scaling depends on the type of PBC
            if options["-hole"]:
                if ("square".startswith(options["-pbc"].value) or
                    "rectangular".startswith(options["-pbc"].value)):
                    scale = 1+options["-hole"].value/min(pbcx,pbcy)
                else:
                    area  = options["-hole"].value**2/math.cos(math.pi/6)
                    scale = 1+area/(pbcx*pbcy)
                pbcx, pbcy = scale*pbcx, scale*pbcy

            pbcx = pbcSetX and pbcSetX[0] or pbcx
            pbcy = pbcSetY and pbcSetY[1] or pbcy

            ## 2. Shift of protein relative to the membrane center
            if options["-dm"]:
                if options["-dm"].value < 0:
                    zshift = options["-dm"].value - max(zip(*prot.coord)[2])
                else:
                    zshift = options["-dm"].value - min(zip(*prot.coord)[2])
            elif not options["-center"]:
                zshift = -shift[2]
            else:
                zshift = 0

            # Now we center the system in the rectangular
            # brick corresponding to the unit cell
            # If -center is given, also center z in plane
            prot += (0.5*pbcx, 0.5*pbcy, zshift)


        # And we collect the atoms
        protein.atoms.extend(prot.atoms)
        protein.coord.extend(prot.coord)

    # Extract the parts of the protein that are in either leaflet
    prot_up,prot_lo = [],[]
    for ix,iy,iz in protein.coord:
        if   iz > 0 and iz <  2.4:
            prot_up.append((ix,iy))
        elif iz < 0 and iz > -2.4:
            prot_lo.append((ix,iy))

    # Current residue ID is set to that of the last atom
    resi = protein.atoms[-1][2]

atid      = len(protein)+1
molecules = []

rx, ry, rz = pbcx+1e-8, pbcy+1e-8, pbcz+1e-8

#################
## 2. MEMBRANE ##
#################

membrane = Structure()

if lipL or options["-lf"]:
    # Lipids are added on grid positions, using the prototypes defined above.
    # If a grid position is already occupied by protein, the position is untagged.


    #print "lipd=",lipd
    # Number of lipids in x and y if there were no solute
    # lipd = sqrt( area of lipid )
    lipids_x = int(pbcx/lipd+0.5)
    lipdx    = pbcx/lipids_x
    rlipx    = list(range(lipids_x))
    lipids_y = int(pbcy/lipd+0.5)
    lipdy    = pbcy/lipids_y
    rlipy    = list(range(lipids_y))

    print("; X: %.3f (%d lipids) Y: %.3f (%d lipids)"%(pbcx,lipids_x,pbcy,lipids_y), file=sys.stderr)

    # Set up grids to check where to place the lipids
    grid_up = [[0 for j in rlipy] for i in rlipx]
    grid_lo = [[0 for j in rlipy] for i in rlipx]

    lipidZbd = [[[ 0 for i in rlipy] for j in rlipx ] for k in range(-1,2) ]
    for i in rlipx:
      for j in rlipy:
        lipidZbd[1][i][j]=-100
        lipidZbd[-1][i][j]=100

    # If there is a protein, mark the corresponding cells as occupied
    if protein:
        for i in prot_up:
            grid_up[ int(lipids_x*i[0]/rx)%lipids_x ][ int(lipids_y*i[1]/ry)%lipids_y ] += 1
        for i in prot_lo:
            grid_lo[ int(lipids_x*i[0]/rx)%lipids_x ][ int(lipids_y*i[1]/ry)%lipids_y ] += 1

    # Determine which cells to consider occupied, given the fudge factor
    # The array is changed to boolean type here
    maxd    = float(max([max(i) for i in grid_up+grid_lo]))
    if  maxd == 0:
        if protein:
            print("; The protein seems not to be inside the membrane.", file=sys.stderr)
            print("; Run with -orient to put it in.", file=sys.stderr)
        maxd = 1


    fudge   = options["-fudge"].value
    grid_up = [[(j/maxd) <= fudge for j in i] for i in grid_up]
    grid_lo = [[(j/maxd) <= fudge for j in i] for i in grid_lo]


    # If we don't want lipids inside of the protein
    # we also mark everything from the center up to the first cell filled
    if not options["-ring"]:

        # Upper leaflet
        marked = [(i,j) for i in rlipx for j in rlipy if not grid_up[i][j]]
        if marked:
            # Find the center
            cx,cy  = [float(sum(i))/len(marked) for i in zip(*marked)]
            for i,j in marked:
                md = int(abs(i-cx)+abs(j-cy)) # Manhattan length/distance
                for f in range(md):
                    ii = int(cx+f*(i-cx)/md)
                    jj = int(cy+f*(j-cy)/md)
                    grid_up[ii][jj] = False

        # Lower leaflet
        marked = [(i,j) for i in rlipx for j in rlipy if not grid_lo[i][j]]
        if marked:
            # Find the center
            cx,cy  = [float(sum(i))/len(marked) for i in zip(*marked)]
            for i,j in marked:
                md = int(abs(i-cx)+abs(j-cy)) # Manhattan length
                for f in range(md):
                    ii = int(cx+f*(i-cx)/md)
                    jj = int(cy+f*(j-cy)/md)
                    grid_lo[ii][jj] = False


    # If we need to add a hole, we simply flag the corresponding cells
    # as occupied. The position of the hole depends on the type of PBC,
    # to ensure an optimal arrangement of holes around the protein. If
    # there is no protein, the hole is just put in the center.
    if options["-hole"]:
        if protein:
            if ("square".startswith(options["-pbc"].value) or
                "rectangular".startswith(options["-pbc"].value)):
                hx,hy = (0,0)
            else:
                hx,hy = (0,int(lipids_y*math.cos(math.pi/6)/9+0.5))
        else:
            hx,hy = (int(0.5*lipids_x), int(0.5*lipids_y))
        hr = int(options["-hole"].value/min(lipdx,lipdy)+0.5)
        ys = int(lipids_x*box[1][0]/box[0][0]+0.5)
        print("; Making a hole with radius %f nm centered at grid cell (%d,%d)"%(options["-hole"].value,hx, hy), hr, file=sys.stderr)
        hr -= 1
        for ii in range(hx-hr-1,hx+hr+1):
            for jj in range(hx-hr-1,hx+hr+1):
                xi, yj = ii, jj
                if (ii-hx)**2+(jj-hy)**2 < hr**2:
                    if jj < 0:
                        xi += ys
                        yj += lipids_y
                    if jj >= lipids_y:
                        xi -= ys
                        yj -= lipids_y
                    if xi < 0:
                        xi += lipids_x
                    if xi >= lipids_x:
                        xi -= lipids_x
                    grid_lo[xi][yj] = False
                    grid_up[xi][yj] = False


    #inject farnesyl where ras is placed
    far_mxy=[]
    if rasl:
        for ras in rasl:
          xx=float(ras[2])
          yy=float(ras[3])
          ii=int(lipids_x*xx/rx)%lipids_x
          jj=int(lipids_y*yy/ry)%lipids_y
          far_mxy.append((ii*pbcx/lipids_x,jj*pbcy/lipids_y,ii,jj))
          if ras[1] == 'l':
              grid_lo[ii][jj] = False
          elif ras[1] == 'u':
              grid_up[ii][jj] = False
          #print >>sys.stderr, "; location ",ras[1],"  ",ii,"  ",jj,"  deleted for ras\n"


    # Set the XY coordinates
    # To randomize the lipids we add a random number which is used for sorting
    random.seed()

    patch_x = 1
    patch_y = 1

    #todo handle two differnet patch sizes
    if options["-lf"]:
       #patch_x,patch_y, lipU_file = parse_patch(options["-lf"].value)
       patch_x,patch_y, lipL_file = parse_patch(options["-lf"].value)
    if options["-uf"]:
       #patch_x,patcy_y, lipL_file = parse_patch(options["-uf"].value)
       patch_x,patcy_y, lipU_file = parse_patch(options["-uf"].value)

    asymV  = options["-asym"].value or 0
    Npatch = patch_x * patch_y

    #  Add code to change -asym - now for all subpaces
    if ( asymV != 0):
        asymAbs = int(abs(asymV))
        asymAvg = int(asymAbs / Npatch)
        asymRem = asymAbs % Npatch
        asymArray = [ (asymAvg + (i<asymRem)) for i in range(Npatch) ]
        for i in range(Npatch):
            tmp=int(random.random()*float(Npatch))%Npatch
            Vtmp = asymArray[i];
            asymArray[i]= asymArray[tmp]
            asymArray[tmp]=Vtmp
        if ( asymV < 0 ):
            asymArray=[ (-1)*val for val in asymArray ]
    else:
        asymArray= [0] * Npatch
    #print 'asymArray=',asymArray

    #when lipdis_? is not multiple of patch_?, first (lipids_? % patch_?) gets one more lipids
    sub_lipids_x = [(int(lipids_x / patch_x) + (i< lipids_x%patch_x)) for i in range(patch_x)]
    sub_lipids_y = [(int(lipids_y / patch_y) + (i< lipids_y%patch_y)) for i in range(patch_y)]

    #print "lipids ",lipids_x,lipids_y
    #print "sub_lipids ",sub_lipids_x,sub_lipids_y
    offset_y = 0
    for jj,len_y in zip(range(patch_y),sub_lipids_y):
    #offset_x = 0
    #for ii,len_x in zip(range(patch_x),sub_lipids_x):

        #offset_y = 0
        #for jj,len_y in zip(range(patch_y),sub_lipids_y):
        offset_x = 0
        for ii,len_x in zip(range(patch_x),sub_lipids_x):

            # print "start ii,jj ",ii,jj,"================="

            upper, lower = [], []
            for i in range(offset_x, offset_x + len_x):
                for j in range(offset_y, offset_y + len_y):
                    if grid_up[i][j]:  #ck: False if there is a hole
                        upper.append((random.random(),i*pbcx/lipids_x,j*pbcy/lipids_y,i,j))  #ck: assign a coordinate to each grid point. lipids_x points from 0 to 1
                    if grid_lo[i][j]:
                        lower.append((random.random(),i*pbcx/lipids_x,j*pbcy/lipids_y,i,j))

            offset_x = offset_x + len_x
            #offset_y = offset_y + len_y
            # Sort on the random number
            upper.sort()
            lower.sort()

            print("; %d lipids in upper leaflet, %d lipids in lower leaflet"%(len(upper),len(lower)), file=sys.stderr)

            # print "ii,jj,upper=",ii,jj,upper
            # print "ii,jj,lower=",ii,jj,lower

            # Extract coordinates, taking asymmetry in account #ck:what is it doing? what if asym is not integer?
            # asym  = options["-asym"].value or 0
            # upper = [i[1:] for i in upper[max(0, asym):]]
            # lower = [i[1:] for i in lower[max(0,-asym):]]
            #upper = [ i[1:] for i in upper[max(0,asymArray[ii*patch_x+jj]):] ]
            #lower = [ i[1:] for i in lower[max(0,-asymArray[ii*patch_x+jj]):] ]
            upper = [ i[1:] for i in upper[max(0,asymArray[jj*patch_y+ii]):] ]
            lower = [ i[1:] for i in lower[max(0,-asymArray[jj*patch_y+ii]):] ]

            # Types of lipids, relative numbers, fractions and numbers

            #if(lipU_file):
            if(options["-uf"].value):
                 lipU = lipU_file.pop(0)
            #if(lipL_file):
            if(options["-lf"].value):
                 lipL = lipL_file.pop(0)

            lipU = lipU or lipL #ck: what's this?

            # Upper leaflet (+1)
            lipU, numU = list(zip(*[ parse_mol(i) for i in lipU ]))
            #print("HII-out")
            #print(lipU)
            #print(numU)
            totU       = float(sum(numU))  #ck: float berfore sum?
            #print(totU)
            num_up     = [int(len(upper)*i/totU) for i in numU]
            #print(num_up)
            num_up     = [(i + (j<(len(upper)-sum(num_up)))) for i,j in zip(num_up,list(range(len(num_up))))]
            #print(num_up)
            lip_up     = [l for i,l in zip(num_up,lipU) for j in range(i)]  #ck ?
            #print(lip_up)
            leaf_up    = ( 1,list(zip(lip_up,upper)))
            #print(leaf_up)

#            print "number of up=",sum(num_up),"   number of upper=",len(upper),"  number of lip_up=",len(lip_up)

            # Lower leaflet (-1)
            lipL, numL = list(zip(*[ parse_mol(i) for i in lipL ]))
            totL       = float(sum(numL))

            # 2020.06.19 Helgi - added code by Tomas Oppelstrup to fix asym rounding in sub-patches
            #num_lo     = [int(len(lower)*i/totL) for i in numL]
            #num_lo     = [(i + (j<(len(lower)-sum(num_lo)))) for i,j in zip(num_lo,list(range(len(num_lo))))]
            num_lo      = subsample(len(lower),numL)

            lip_lo     = [l for i,l in zip(num_lo,lipL) for j in range(i)]
            leaf_lo    = (-1,list(zip(lip_lo,lower)))

            molecules  += list(zip(lipU,num_up)) + list(zip(lipL,num_lo))

            kick       = options["-rand"].value

            # Build the membrane
            for leaflet,leaf_lip in [leaf_up,leaf_lo]:
                for lipid, pos in leaf_lip:
                    # Increase the residue number by one
                    resi += 1
                    # Set the random rotation for this lipid
                    rangle   = random.random()*math.pi
                    rcos     = math.cos(rangle)
                    rsin     = math.sin(rangle)
                    # Fetch the atom list with x,y,z coordinates
                    atoms    = list(zip(lipidsa[lipid][1].split(),lipidsx[lipidsa[lipid][0]],lipidsy[lipidsa[lipid][0]],lipidsz[lipidsa[lipid][0]]))
                    # Only keep atoms appropriate for the lipid
                    at,ax,ay,az = list(zip(*[i for i in atoms if i[0] != "-"]))
                    # The z-coordinates are spaced at 0.3 nm,
                    # starting with the first bead at 0.15 nm
                    h = phaseHeight.height(pos[0]+lipdx/2,pos[1]+lipdy/2,pbcx,pbcy); #default is 0.15
                    az       = [ h + leaflet*(0.15 + 0.3*(i-min(az))) for i in az ]

                    if (leaflet==1):
                      lipidZbd[1][pos[2]][pos[3]] = max(az)
                    else:
                      lipidZbd[-1][pos[2]][pos[3]] = min(az)

                    xx       = list(zip( ax,ay ))
                    nx       = [rcos*i-rsin*j+pos[0]+lipdx/2+random.random()*kick for i,j in xx]
                    ny       = [rsin*i+rcos*j+pos[1]+lipdy/2+random.random()*kick for i,j in xx]
                    #print >>sys.stderr, nx
                    # Add the atoms to the list
                    for i in range(len(at)):
                        atom  = "%5d%-5s%5s%5d"%(resi,lipid,at[i],atid)
                        membrane.coord.append((nx[i],ny[i],az[i]))
                        membrane.atoms.append((at[i],lipid,resi,0,0,0))
                        atid += 1
            # print "i,j,mem=",i,j,membrane.coord
        #end of jj loop over patch
        #offset_x = offset_x + len_x
        offset_y = offset_y + len_y
    #end of ii loop over patch

    # -ras RAS.gro side x y RAS.gro x y ...
    # later a file can be fed with similar format
    # -rasf RAS.conf
    # RAS.conf
    # RAS.gro side x y
    # RAS.gro side x y
    # ....
    if rasl:
        #m = re.compile("FAR")
        m = re.compile('CYF +F\d')
        memz = [ i[2] for i in membrane.coord ]
        memz_min = min(memz)
        memz_max = max(memz)
        leaflet_dict = {'l':-1, 'u':1}
        for ras in rasl:

           molecules  += [("Protein",1)]
           prot = ras[0];
           #print 'ras=',ras
           rcos=1
           rsin=0
           #generate farnesyl
           leaflet = leaflet_dict[ras[1]]
           lipid='CYF'
           pos = far_mxy.pop(0)
           atoms    = list(zip(lipidsa[lipid][1].split(),lipidsx[lipidsa[lipid][0]],lipidsy[lipidsa[lipid][0]],lipidsz[lipidsa[lipid][0]]))
           at,ax,ay,az = list(zip(*[i for i in atoms if i[0] != "-"]))
           # The z-coordinates are spaced at 0.3 nm,
           # starting with the first bead at 0.15 nm
           h = phaseHeight.height(pos[0]+lipdx/2,pos[1]+lipdy/2,pbcx,pbcy); #default is 0.15
           az       = [h + leaflet*(0.15 + 0.3*(i-min(az))) for i in az ]
           if (leaflet==1):
             lipidZbd[1][pos[2]][pos[3]] = max(az)
           else:
             lipidZbd[-1][pos[2]][pos[3]] = min(az)

           xx       = list(zip( ax,ay ))
           nx       = [rcos*i-rsin*j+pos[0]+lipdx/2+random.random()*kick for i,j in xx]
           ny       = [rsin*i+rcos*j+pos[1]+lipdy/2+random.random()*kick for i,j in xx]
           #print >>sys.stderr, nx

           #shift
           radius = 8  #10 lipid square radius
           membd=[ lipidZbd[leaflet][(pos[2]+ii+lipids_x)%lipids_x][(pos[3]+jj+lipids_y)%lipids_y] for ii in range(-radius,radius+1) for jj in range(-radius,radius+1) ]
           restz = [i[6] for i in prot.atoms if not m.match(i[1]+i[0]) ]
           zmin = min(restz)
           zmax = max(restz)
           if ras[1] == 'l':
#             zshift = memz_min-zmax-2.
             zshift = min(membd)-zmax-0.3
           elif ras[1] == 'u':
#             zshift = memz_max-zmin+2.
             zshift = max(membd)-zmin+0.3

           far = [i[4:6] for i in prot.atoms if m.match(i[1]+i[0]) ]
           if (len(far) != 4 ):
             print("not clear what to do with a farnesyl with ",len(far)," atmos\n")
           xy = [ sum(i)/len(i) for i in zip(*far)]

           xshift = pos[0] + lipdx/2 - xy[0]
           yshift = pos[1] + lipdy/2 - xy[1]

           tmp = prot.coord
           prot.coord = [(q[0]+xshift,q[1]+yshift,q[2]+zshift) for q in tmp]

           #overwrite back farnesyl
           far = [i for i,item in enumerate(prot.atoms) if m.match(item[1]+item[0]) ]
           for i,x,y,z in zip(far,nx,ny,az):
               prot.coord[i] = ( x , y , z )
               #print >>sys.stderr, "x,y,z=",x,y,z



           atid += len(prot)
           protein.coord.extend(prot.coord)
           protein.atoms.extend(((x[0],x[1],x[2]+resi,x[3],x[4],x[5]) for x in prot.atoms))

           resi = protein.atoms[-1][2]
    #end of if rasl

    # Now move everything to the center of the box before adding solvent
    mz  = pbcz/2
    z   = [ i[2] for i in protein.coord+membrane.coord ]
    mz -= (max(z)+min(z))/2
    #print >>sys.stderr, "mz=",mz
    protein += (0,0,mz)
    membrane += (0,0,mz)

# The box dimensions are now (likely) set.
# If a protein was given, it is positioned in the center of the
# rectangular brick.

#xyz   = protein.coord+membrane.coord
#minxyz = [min(i) for i in zip(*xyz)]
#maxxyz = [max(i) for i in zip(*xyz)]
#pbcx = maxxyz[0] - minxyz[0]
#pbcy = maxxyz[1] - minxyz[1]
#pbcz = maxxyz[2] - minxyz[2]

# Set the lattice vectors
if ("rectangular".startswith(options["-pbc"].value) or
    "square".startswith(options["-pbc"].value) or
    "cubic".startswith(options["-pbc"].value)):
    box    = [[pbcx,0,0],[0,pbcy,0],[0,0,pbcz]]
elif ((not lipL) or (not options["-lf"])):
    # Rhombic dodecahedron with square XY plane
    box    = [[pbcx,0,0],[0,pbcy,0],[0.5*pbcx,0.5*pbcx,pbcz]]
elif "hexagonal".startswith(options["-pbc"].value):
    box    = [[pbcx,0,0],[math.sin(math.pi/6)*pbcx,pbcy,0],[0,0,pbcz]]
else: # optimal packing; rhombic dodecahedron with hexagonal XY plane
    box    = [[pbcx,0,0],[math.sin(math.pi/6)*pbcx,pbcy,0],[pbcx/2,pbcy/3,pbcz]]

# Override lattice vectors if they were set explicitly
box[0] = pbcSetX or box[0]
box[1] = pbcSetY or box[1]
box[2] = pbcSetZ or box[2]

grobox = (box[0][0],box[1][1],box[2][2],
          box[0][1],box[0][2],box[1][0],
          box[1][2],box[2][0],box[2][1])





################
## 3. SOLVENT ##
################

# Charge of the system so far

last = None
mcharge = 0
for j in membrane.atoms:
    if not j[0].strip().startswith('v') and j[1:3] != last:
        mcharge += charges.get(j[1].strip(),0)
    last = j[1:3]

last = None
pcharge = 0
for j in protein.atoms:
    if not j[0].strip().startswith('v') and j[1:3] != last:
        pcharge += charges.get(j[1].strip(),0)
    last = j[1:3]

# We have strange protein so manually add protein carge
#mcharge = sum([charges.get(i[0].strip(),0) for i in set([j[1:3] for j in membrane.atoms])])
pcharge = sum([charges.get(i[0].strip(),0) for i in set([j[1:3] for j in protein.atoms if not j[0].strip().startswith('v')])])

charge  = mcharge + pcharge

if protein:
    print("; Charge of protein: %f" % pcharge, file=sys.stderr)
if membrane:
    print("; Charge of membrane: %f" % mcharge, file=sys.stderr)
print("; Total charge: %f" % charge, file=sys.stderr)

if solv:

    # Set up a grid
    d        = 1/options["-sold"].value

    nx,ny,nz = int(1+d*pbcx),int(1+d*pbcy),int(1+d*pbcz)
    dx,dy,dz = pbcx/nx,pbcy/ny,pbcz/nz
    excl,hz  = int(nz*options["-excl"].value/pbcz), int(0.5*nz)

    #For now, populate the water grid with lipid molecules
    #convert lipid grid to solvent grid
    #use the ratio of the sizes

    #for each water grid
    #calculate coordi
    #floor and ceil to find nearest lipid grid point
    #pick max for the height



    zshift   = 0
    if membrane:
        memz   = [i[2] for i in membrane.coord]
        midz   = (max(memz)+min(memz))/2
        hz     = int(nz*midz/pbcz)  # Grid layer in which the membrane is located
        zshift = (hz+0.5)*nz - midz # Shift of membrane middle to center of grid layer


    # Initialize a grid of solvent, spanning the whole cell
    # Exclude all cells within specified distance from membrane center
#    radius=1
#    solvExU = [[ 0 for j in xrange(ny)] for i in xrange(nx)]
#    solvExL = [[ 0 for j in xrange(ny)] for i in xrange(nx)]

#    print 'lipids_x=%d lipids_y=%d\n' % (lipids_x,lipids_y)
#    for i in xrange(nx):
#      for j in xrange(ny):
#        lii = int(float(i)/(lipdx*d))
#        lij = int(float(j)/(lipdy*d))
#        solvExU[i][j] = (max([ lipidZbd[ 1][(lii+ii+lipids_x)%lipids_x][(lij+jj+lipids_y)%lipids_y] for ii in range(-radius,radius+1) for jj in range(-radius,radius+1) ]) + mz + 0.1 ) / pbcz * nz;
#        solvExL[i][j] = (min([ lipidZbd[-1][(lii+ii+lipids_x)%lipids_x][(lij+jj+lipids_y)%lipids_y] for ii in range(-radius,radius+1) for jj in range(-radius,radius+1) ]) + mz - 0.1 ) / pbcz * nz;
#        print 'solvExU i=%d j=%d lii=%d lij=%d %f\n' % (i,j,lii,lij,solvExU[i][j])
#        print 'solvExL',i,' ',j,' ',solvExL[i][j]
#    grid   = [[[i < solvExL[k][j] or i > solvExU[k][j] for i in xrange(nz)] for j in xrange(ny)] for k in xrange(nx)]

    def hconv(ii,jj):
      #return (phaseHeight.height(float(ii)/nx*pbcx,float(jj)/ny*pbcy,pbcx,pbcy) + mz)/pbcz*nz
      return (phaseHeight.height(float(ii)/nx*pbcx+0.5*lipdx,float(jj)/ny*pbcy+0.5*lipdy,pbcx,pbcy) + mz)/pbcz*nz

    grid   = [[[i < hconv(k,j)-excl or i > hconv(k,j)+excl for i in range(nz)] for j in range(ny)] for k in range(nx)]

    # Flag all cells occupied by protein or membrane
    for x,y,z in protein.coord+membrane.coord:
        if z >= pbcz:
            x -= box[2][0]
            y -= box[2][1]
            z -= box[2][2]
        if z < 0:
            x += box[2][0]
            y += box[2][1]
            z += box[2][2]
        if y >= pbcy:
            x -= box[1][0]
            y -= box[1][1]
        if y < 0:
            x += box[1][0]
            y += box[1][1]
        if x >= pbcx:
            x -= box[0][0]
        if x < 0:
            x += box[0][0]
        grid[int(nx*x/rx)][int(ny*y/ry)][int(nz*z/rz)] = False

    # Set the center for each solvent molecule
    kick = options["-solr"].value
    grid = [ (R(),(i+0.5+R()*kick)*dx,(j+0.5+R()*kick)*dy,(k+0.5+R()*kick)*dz)
             for i in range(nx) for j in range(ny) for k in range(nz) if grid[i][j][k] ]

    # Sort on the random number
    grid.sort()

    # 'grid' contains all positions on which a solvent molecule can be placed.
    # The number of positions is taken as the basis for determining the salt concentration.
    # This is fine for simple salt solutions, but may not be optimal for complex mixtures
    # (like when mixing a 1M solution of this with a 1M solution of that

    # First get names and relative numbers for each solvent
    solnames, solnums = list(zip(*[ parse_mol(i) for i in solv ]))
    solnames, solnums = list(solnames), list(solnums)
    totS       = float(sum(solnums))

    # Set the number of ions to add
    nna, ncl = 0, 0
    if options["-salt"]:

        # If the concentration is set negative, set the charge to zero
        if options["-salt"].value.startswith("-"):
            charge = 0
            options["-salt"].value = -float(options["-salt"].value)
        else:
            options["-salt"].value = float(options["-salt"].value)

        # Determine charge to use, either determined or given on command line
        # HII we keep lipid charge auto calculate but add supplied charge
        if options["-charge"].value != "0":
            #charge = (options["-charge"].value != "auto") and int(options["-charge"].value) or charge
            charge = charge + int(options["-charge"].value)
        #else:
        #    charge = charge

        # Determine number of sodium and chloride to add
        concentration = options["-salt"].value
        nsol = ("SPC" in solnames and 1 or 4)*len(grid)

        # HII this is a dirty hack needs to be fixed
        # add salt independent on charge then addjust +Na for -salt > 0 or +Cl for -salt < 0
        #ncl  = max(max(0,charge),int(.5+.5*(concentration*nsol/(27.7+concentration)+charge)))
        #nna  = ncl - charge
        ncl  = max(0,int(.5+.5*(concentration*nsol/(27.7+concentration))))
        nna  = ncl
        #if options["-charge"].value != "0":
        #    tempCharge = int(options["-charge"].value)
        #    if tempCharge > 0:
        #        nna = nna + tempCharge
        #    else:
        #        ncl = ncl - tempCharge
        if charge < 0:
            nna = nna - charge
        else:
            ncl = ncl + charge
    # Correct number of grid cells for placement of solvent
    ngrid   = len(grid) - nna - ncl
    num_sol = [int(ngrid*i/totS) for i in solnums]


    # Add salt to solnames and num_sol
    if nna:
        solnames.append("NA")
        num_sol.append(nna)
        solv.append("NA")
    if ncl:
        solnames.append("CL")
        num_sol.append(ncl)
        solv.append("CL")


    # Names and grid positions for solvent molecules
    solvent    = list(zip([s for i,s in zip(num_sol,solnames) for j in range(i)],grid))


    # Extend the list of molecules (for the topology)
    molecules.extend(list(zip(solnames,num_sol)))


    # Build the solvent
    sol = []
    for resn,(rndm,x,y,z) in solvent:
        resi += 1
        solmol = solventParticles.get(resn)
        if solmol and len(solmol) > 1:
            # Random rotation (quaternion)
            u,  v,  w       = random.random(), 2*math.pi*random.random(), 2*math.pi*random.random()
            s,  t           = math.sqrt(1-u), math.sqrt(u)
            qw, qx, qy, qz  = s*math.sin(v), s*math.cos(v), t*math.sin(w), t*math.cos(w)
            qq              = qw*qw-qx*qx-qy*qy-qz*qz
            for atnm,(px,py,pz) in solmol:
                qp = 2*(qx*px + qy*py + qz*pz)
                rx = x + qp*qx + qq*px + qw*(qy*pz-qz*py)
                ry = y + qp*qy + qq*py + qw*(qz*px-qx*pz)
                rz = z + qp*qz + qq*pz + qw*(qx*py-qy*px)
                sol.append(("%5d%-5s%5s%5d"%(resi%1e5,resn,atnm,atid%1e5),(rx,ry,rz)))
                atid += 1
        else:
            sol.append(("%5d%-5s%5s%5d"%(resi%1e5,resn,solmol and solmol[0][0] or resn,atid%1e5),(x,y,z)))
            atid += 1
else:
    sol = []


## Write the output ##

# Open the output stream
oStream = options["-o"] and open(options["-o"].value,"w") or sys.stdout

# Print the title
if membrane.atoms:
    if options["-uf"]:
       #patch_x,patch_y, lipU_file = parse_patch()
       #title  = "INSANE! Membrane UpperLeaflet"+":"+str(patch_x)+"x"+str(patch_y)+"\r\n"+"\n".join(["    "+str(i) for i in lipU_file])
       title  = "INSANE! Membrane UpperLeaflet from file:"+options["-uf"].value
    else:
       title  = "INSANE! Membrane UpperLeaflet>"+":".join(lipU)+"="+":".join([str(i) for i in numU])

    #title +="\r\n"

    if options["-lf"]:
       #patch_x,patch_y, lipL_file = parse_patch(options["-lf"].value)
       #title += "        Membrane LowerLeaflet"+":"+str(patch_x)+"x"+str(patch_y)+"\r\n"+"\n".join(["    "+str(i) for i in lipL_file])
       title += "        Membrane LowerLeaflet from file:"+options["-lf"].value
    else:
       title += " LowerLeaflet>"+":".join(lipL)+"="+":".join([str(i) for i in numL])

    if protein:
        title = "Protein in " + title
else:
    title = "Insanely solvated protein."

print(title, file=oStream)

# Print the number of atoms
print("%5d"%(len(protein)+len(membrane)+len(sol)), file=oStream)

# Print the atoms
id = 1
if membrane:
#    print "membrane is prepared"
    for i in range(len(membrane)):
        at,rn,ri = membrane.atoms[i][:3]
        x,y,z    = membrane.coord[i]
        # HII 2013.11.13 to to fix .gro format error
        #oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri,rn,at,id%1e5,x,y,z))
        oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri%1e5,rn,at,id%1e5,x,y,z))
        id += 1
if protein:
    for i in range(len(protein)):
        at,rn,ri = protein.atoms[i][:3]
        x,y,z    = protein.coord[i]
        oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri,rn,at,id,x,y,z))
        id += 1
if sol:
    # Print the solvent
    print("\n".join([i[0]+"%8.3f%8.3f%8.3f"%i[1] for i in sol]), file=oStream)

# Print the box
print("%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n"%grobox, file=oStream)

if options["-p"]:
    # Write a rudimentary topology file
    top = open(options["-p"].value,"w")
    print('#include "martini.itp"\n', file=top)
    print('[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%title, file=top)
    print("\n".join("%-10s %7d"%i for i in molecules), file=top)
#    if protein:
#        print >>top, "%-10s %5d"%("Protein",1)
    top.close()
else:
    print("\n".join("%-10s %7d"%i for i in molecules), file=sys.stderr)
