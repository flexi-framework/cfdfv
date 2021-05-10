#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import subprocess

def execute_cfdfv(flexi_path, prm_path, analyze_fcts = None, log = True, \
      mpi_procs = 1) :
   if mpi_procs == 1 :
      cmd = []
   else :
      cmd = ["mpirun", "-np", "%d" % mpi_procs]
   cmd.append(flexi_path)
   cmd.append(prm_path)
   p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
   p.wait()

   lines = p.stdout.readlines()
   # for line in lines :
      # print line
   if log :
      log_path = os.path.splitext(os.path.basename(prm_path))[0] + ".log"
      f = open(log_path, 'w')
      for line in lines :
         f.write(line)
      f.close()

   if analyze_fcts :
      results = []
      if type(analyze_fcts) != list : 
         return analyze_fcts(lines)
      for analyze_fct in analyze_fcts :
         results.append(analyze_fct(lines))
      return results
