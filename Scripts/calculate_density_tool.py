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
# calculate_density.py (version 1.0 beta)
# Created on: Mon Nov 24 2008 04:08:09 PM
#
# Written by Brent L. Brock, Landscape Ecologist, Craighead Environmental Research Institute
#
# Calculates approximate housing density that will satisfy specified wildlife conservation targets assuming a random distribution of
# houses on the landscape
#
#This module depends on functionlib.py which contains functions required for processing this script
# 
# Usage: calculate_density <outWorkspace> <outTable> <Number_of_Iterations_to_Run> <aExtent> <minCoreSize> <distDistance> <numHouses> 
# ---------------------------------------------------------------------------

# Import system modules
import sys, os, arcgisscripting
import functionlib as fl

# Create the Geoprocessor object
gp = arcgisscripting.create()
##gp.Workspace = None


# Set the necessary product code
gp.SetProduct("ArcInfo")

# Load required toolboxes...
gp.overwriteoutput = 1

# Script arguments...
table = gp.GetParameterAsText(0)
table = str(table).rsplit("\\",1)
outWorkspace = table[0]
outTable = table[1]

i = gp.GetParameterAsText(1) # Number of iterations to run
aExtent = gp.GetParameterAsText(2)
## FIX:  check for polygon topology
constraintLayer = gp.GetParameterAsText(3)
minCoreSize = gp.GetParameterAsText(4) 
distDistance = gp.GetParameterAsText(5)
target = gp.GetParameterAsText(6)
outputUnits = gp.GetParameterAsText(11)

##i = "3" # Number of iterations to run
##aExtent = r"w:\gis\inputs\boundary.shp"
#### FIX:  check for polygon topology
##constraintLayer = r"w:\gis\inputs\constraint.shp"
##minCoreSize = "3000000" #units in square meters
##distDistance = "500" #units in meters
##target = "50" # initial number of houses to run

# Create temporary workspaces
tWorkspace = fl.CreateTempWorkspace(outWorkspace)
sWorkspace = fl.CreateTempWorkspace(outWorkspace)
gp.ScratchWorkspace = tWorkspace


# Local variables...

area = minCoreSize.split(" ")
minCoreSize = str(fl.ConvertAreaToMeters(area[0], area[1]))
if outputUnits == "" or outputUnits == "#":
    outputUnits = area[1]
gp.AddMessage("Minimum Patch Size: " + str(minCoreSize) + " square meters...")
del area

distance = distDistance.split(" ")
distDistance = str(fl.ConvertDistanceToMeters(distance[0], distance[1]))
gp.AddMessage("Disturbance Distance: " + str(distDistance) + " meters...")
del distance

cExtent = tWorkspace+"\\aExtent_copy.shp" #Temporary copy of analysis extent for manipulation
nTable = outWorkspace+"\\"+outTable #Final output table
i = int(i)
endSim = False
target = int(target)
low = target - (target * 0.05)
high = target + (target * 0.05)
currentDirectory = sys.path[0]


# If constraint layer is specified, make copy and clip to aExtent
if constraintLayer == "#" or constraintLayer == "":
    constraintLayer = aExtent
else:
    gp.MakeFeatureLayer_management(constraintLayer, "tempLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
    gp.Clip_analysis("tempLayer", aExtent, sWorkspace + "\\clip_constraint.shp", "")
    constraintLayer = sWorkspace + "\\clip_constraint.shp"
    
##Calculate Total Area and store in variable 'tArea'##
gp.AddMessage("Caculating area of analysis extent...")
# Make copy of analysis extent layer to modify
try:
    gp.Copy_management(aExtent, cExtent, "ShapeFile")
except:
    a = gp.MakeFeatureLayer_management(aExtent, "aLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
    gp.CopyFeatures(a, cExtent)
    
# Add 'Area' Field and calculate Area
fl.CalcArea(cExtent)
fl.CalcArea(constraintLayer)
# Calculate total area of analysis extent
tArea = fl.CalcTotalArea(gp.ScratchWorkspace, cExtent)
gp.AddMessage("Total Analysis Area = " + str(tArea))
cArea = fl.CalcTotalArea(gp.ScratchWorkspace, constraintLayer)
gp.AddMessage("Total Buildable Area = " + str(cArea))
prcntBuildable = int(cArea/tArea)
  
    
## End Calculate Total Area ##

MAX = int(tArea/float(minCoreSize)) #maximum number of houses for simulation
MIN = 1 #minimum number of houses for simulation


# If number of houses in simulation will be very large, give user a chance to abort
if MAX > 50000:
    from Tkinter import *
    import tkMessageBox
    root = Tk()
    root.withdraw()
    if tkMessageBox.askokcancel(
    "Housing Density",
    "The parameters entered result in a maximum number of houses of " + str(MAX) + " and may take a long time to process.  Click 'Cancel' to quit or 'OK' to continue."):
        pass
    else:
        gp.AddMessage("Simulation aborted by user.  Cleaning up temporary files...")
        fl.CleanFiles(sWorkspace)
        fl.CleanFiles(tWorkspace)
        gp.Delete_management(tWorkspace,"")
        gp.Delete_management(sWorkspace,"")
        sys.exit()
        

## Create New Table to store final results if it doesn't exist... ##

#If the nTable does not exist, create it.
gp.AddMessage("Creating output table " + nTable + "...")
if not gp.Exists(nTable):
    fl.MakeTable(outWorkspace, outTable)
##    # Process: Create Table...
##    gp.CreateTable_management(outWorkspace, outTable, "", "")
##    # Process: Add Fields...
##    gp.AddField_management(nTable, "SUM_SUM_AR", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
##    gp.AddField_management(nTable, "MAX_SUM_AR", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
##    gp.AddField_management(nTable, "MIN_SUM_AR", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
##    gp.AddField_management(nTable, "HOUSES", "LONG", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
##    gp.AddField_management(nTable, "PRCNT_AREA", "DOUBLE", "18", "4", "", "", "NON_NULLABLE", "NON_REQUIRED", "")

# Calculate initial number of houses for random pattern simulations
gp.AddMessage("Calculating number of houses for first simulation run...")
numHouses = MAX / 2 # Set initial seed for number of houses
if numHouses == 0:
    numHouses = 1
numHouses = str(numHouses)

### Check if target is attainable with given parameters
##if prcntBuildable < target:
##    gp.AddWarning("Buildable area is only " + str(prcntBuildable) + "% and is less than the target of " + str(target))
##    gp.AddWarning("Adjust target minimum possible...")
##    target = prcntBuildable
##    numHouses = (int(tArea/int(minCoreSize)))


# Get result of monte carlo simulation
gp.AddMessage("Begin simulation 1...")
result = fl.RunSimulation(i,gp.ScratchWorkspace, aExtent, constraintLayer, numHouses, minCoreSize, distDistance, nTable, tArea, "", "", True)


# Test whether result does not equal target +/- 5%
# If result is outside target range, use "principle of halves" to adjust number of houses
# until a solution is reached

count = 1
lastSim = 0

while result < low or result >= high and not endSim:
    
    if result >= high:
##        if count > 20:
##            gp.AddWarning("Could not reach a solution after 10 simulation runs.  Setting number of houses to MAX for final simulation...")
##            numHouses = MAX
###            result = fl.RunSimulation(i,gp.ScratchWorkspace, aExtent, constraintLayer, numHouses, minCoreSize, distDistance, nTable, tArea, "", "", True)
##            endSim = True
##            break
##        
##        else:
##            MIN = int(numHouses)
##            numHouses = ((MAX - MIN)/2)+ MIN
##            gp.AddMessage("Begin simulation " + str(count + 1) + "...")
##        result = fl.RunSimulation(i,gp.ScratchWorkspace, aExtent, constraintLayer, numHouses, minCoreSize, distDistance, nTable, tArea, "", "", True)
##        count = count + 1

        MIN = int(numHouses)
        numHouses = ((MAX - MIN)/2)+ MIN
        gp.AddMessage("Begin simulation " + str(count + 1) + "...")
        result = fl.RunSimulation(i,gp.ScratchWorkspace, aExtent, constraintLayer, numHouses, minCoreSize, distDistance, nTable, tArea, "", "", True)
        count = count + 1
        #Store descriptors for output if target cannot be reached
        minmax = "maximum"
        obtain = "exceeded"
        # Double the maximum number of houses if numHouses is within 5% of MAX and result is > than high
        if abs(MAX - numHouses) < MAX * 0.05:
            MAX = MAX * 2
        
    if result < low:
        MAX = int(numHouses)
        numHouses = MIN + ((MAX - MIN)/2)
        gp.AddMessage("Begin simulation " + str(count + 1) + "...")
        result = fl.RunSimulation(i,gp.ScratchWorkspace, aExtent, constraintLayer, numHouses, minCoreSize, distDistance, nTable, tArea, "", "", True)
        count = count + 1
        minmax = "minimum"
        obtain = "not obtainable"
        if numHouses < 2:
            endSim = True
            break
            #Store descriptors for output if target cannot be reached
            
            

    # if changing numHouses does not change result more than 1%, end simulation and report results.
    if abs(result - lastSim) < 1:
        endSim = True
        break
    lastSim = result

# Calculate density from the solution
gp.AddMessage("Calculating density...")
density = float(tArea)/int(numHouses)
density = int(round(fl.ConvertSquareMetersToOther(density, outputUnits)))
cDensity = float(cArea)/int(numHouses)
cDensity = int(round(fl.ConvertSquareMetersToOther(cDensity, outputUnits)))
gp.SetParameterAsText(7, density)
gp.SetParameterAsText(8, cDensity)

# Get cores and houses of last iteration result
gp.AddMessage("Copying shapefiles from last iteration to: " + outWorkspace)
rootName = str(outTable).split(".")
rootName = rootName[0]

try:
    outShape1 = gp.Copy_management(tWorkspace + "\\random_temp.shp", outWorkspace + "\\" + rootName + "_example_points.shp", "ShapeFile")
    outShape2 = gp.Copy_Management(tWorkspace + "\\xxcoresLayer.shp", outWorkspace + "\\" + rootName + "_example_cores.shp", "Shapefile")
    gp.SetParameterAsText(9, outShape1)
    gp.SetParameterAsText(10, outShape2)
    params = gp.GetParameterInfo()

    
    ## FIX: symbology is not displaying in ArcMap
    # Set the symbology of the output. 
    #   output      = the output value
    #   params[2] = the output parameter
    #
    gp.AddMessage("The current python script directory is: " + currentDirectory)
    params[8].symbology = currentDirectory + "\\houses.lyr"
    params[9].symbology = currentDirectory + "\\cores.lyr"
    
except:
    gp.AddMessage("Could not copy final outputs to " + outWorkspace)
    
# Clean up temporary workspace
fl.CleanFiles(sWorkspace)
fl.CleanFiles(tWorkspace)
gp.Delete_management(tWorkspace,"")
gp.Delete_management(sWorkspace,"")

if not endSim:
# Display normal result in messagebox
    from Tkinter import *
    import tkMessageBox
    root = Tk()
    root.withdraw()
    tkMessageBox.showinfo(
        "Housing Density",
        "Target: " + str(target) + "% of area conserved as core habitat\n Overall density: " + str(density) + " " + outputUnits + " per house\n Density within buildable area: " + str(cDensity) + " " + outputUnits + " per house")
else:
    # Display truncated result in messagebox
    from Tkinter import *
    import tkMessageBox
    root = Tk()
    root.withdraw()
    tkMessageBox.showinfo(
        "Housing Density",
        "The target of " + str(target) + "% " + obtain + " even at " + minmax + " buildout.\n Which resulted in the following results:\n Area conserved: " + str(int(result)) + "% of avialable habitat\n Overall density: " + str(density) + " " + outputUnits + " per house\n Density within buidable area: " + str(cDensity) + " " + outputUnits + " per house")


    
