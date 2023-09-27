#!/bin/bash


 sed 's/ZN+ ION/ZN   ZN/' system.pdb | sed 's/MG+ ION/MG   MG/' | sed 's/NA+ ION/NA   NA/' | sed 's/CL- ION/CL   CL/' > raf-mem.pdb
