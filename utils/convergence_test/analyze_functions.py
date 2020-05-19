#!/usr/bin/python
# -*- coding: utf-8 -*-

# each analyze-function takes the output-lines of a strukti-run 

# extract the L1 error of the last timestep
def get_last_L1_error(lines) :
   tmp = lines[-6].split(":")[1]
   return [float(x) for x in tmp.split()]

# extract the L2 error of the last timestep
def get_last_L2_error(lines) :
   tmp = lines[-5].split(":")[1]
   return [float(x) for x in tmp.split()]

# extract the L_inf error of the last timestep
def get_last_Linf_error(lines) :
   tmp = lines[-4].split(":")[1]
   return [float(x) for x in tmp.split()]

# extract the computation time
def get_comp_time(lines) :
   tmp = lines[-2].split(":")[1]
   x=tmp.split()
   return float(x[0]) 
