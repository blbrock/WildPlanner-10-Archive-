# ------------------------------------------------------------------------------
#
# Copyright 2011, 2012, 2013 Brent L. Brock and the Craighead Institute
#
# This file is part of Wild Planner.
#
# Wild Planner is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Wild Planner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Wild Planner in the file named LICENSE.TXT.  If not, see <http://www.gnu.org/licenses/>.
#
# ------------------------------------------------------------------------------
#
#compare_resistance.py Parses all *.out files in a directory into
#single table.
import sys, os, glob, arcpy, numpy as np 

p = sys.argv[1]
searchSub = sys.argv[2]
removePath = sys.argv[3]
outFile = sys.argv[4]

# prepare variables for processing
p = p.replace("\\","\\\\")
outFile = outFile.replace("\\", "\\\\")

line = "File name,Mean Ohms\n"

if searchSub == 'true':
    arcpy.AddMessage("Searching subfolders of " + p + "...")
    outfiles = [os.path.join(root, name) 
                 for root, dirs, files in os.walk(p) 
                 for name in files 
                 if name.endswith((".out"))]
else:
    outfiles = glob.glob('*.out')

for files in outfiles: 
    if files.endswith(".out"): 
        arcpy.AddMessage("Processing " + files + "...")
        name = files.split('_resistances_3columns')[0]
        if removePath == 'true':
            name = os.path.basename(name)
        numlines = sum(1 for line in files)
        with open(files) as f: 
            for i, l in enumerate(f): 
                pass
        numline = i + 1

        if numline == 1:
            means = open(files, 'r').readline()
            means = means.split(" ")
#            means = means[len(means) - 1]
        else:
            data = np.loadtxt(files)
            means = data.mean(0)
        means = means[len(means) - 1]
        line = line + name + "," + str(round(float(means),4)) + "\n"

if os.path.exists(outFile):
    os.remove(outFile)
with open(outFile, 'a') as f: f.write(line)
