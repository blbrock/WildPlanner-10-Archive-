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
# evaluate_habitat_cores.py (version 1.0 beta)
# Created on: Mon Nov 24 2008 04:08:09 PM
#
# Written by Brent L. Brock, Landscape Ecologist, Craighead Environmental Research Institute
#
# Calculates percent of core habitat remaining within evaluation area and generates a shapefile
# of the results
#
#This module depends on functionlib.py which contains functions required for processing this script
# 
# ---------------------------------------------------------------------------

# Import system modules
import sys, os, random, arcpy
import functionlib as fl

# Create the Geoprocessor object
#gp = arcgisscripting.create()

# Set the necessary product code
arcpy.SetProduct("ArcInfo")
arcpy.overwriteoutput = 1

# Check out spatial analyst license
arcpy.AddMessage("\tChecking availability of spatial analyst license...")  
if arcpy.CheckExtension("spatial") == "Available":
    arcpy.CheckOutExtension("spatial")         
else:
    arcpy.AddWarning("\tSpatial Analyst is needed to run this tool.  Processing aborted...")
    raise "LicenseError"

# Script arguments...

aExtent = sys.argv[1]
pointLayer = sys.argv[2]
roadLayer = sys.argv[3]
minCoreSize = sys.argv[4] 
stinfDistance = sys.argv[5]
rdinfDistance = sys.argv[6]
outputUnits = sys.argv[7]
outShapefile = sys.argv[8]
s = str(outShapefile).rsplit("\\",1)
outWorkspace = s[0]
outTable = s[1]


# Create temporary workspaces
scratchWS = arcpy.env.scratchWorkspace
if scratchWS:
    #tWorkspace = fl.CreateTempWorkspace(scratchWS)
    #sWorkspace = fl.CreateTempWorkspace(scratchWS)
    sWorkspace = "in_memory"
    tWorkspace = "in_memory"
else:
    # Set output workspace
    arcpy.Workspace = outWorkspace
    tWorkspace = fl.CreateTempWorkspace(outWorkspace)
    sWorkspace = fl.CreateTempWorkspace(outWorkspace)
arcpy.ScratchWorkspace = tWorkspace

# Local variables...
area = minCoreSize.split(" ")
minCoreSize = str(fl.ConvertAreaToMeters(area[0], area[1]))
if outputUnits == "" or outputUnits == "#":
    outputUnits = area[1]
arcpy.AddMessage("Minimum Patch Size: " + str(minCoreSize) + " square meters...")
del area

stdistance = stinfDistance.split(" ")
stinfDistance = str(fl.ConvertDistanceToMeters(stdistance[0], stdistance[1]))
arcpy.AddMessage("Structure Distance: " + str(stinfDistance) + " meters...")
del stdistance

rddistance = rdinfDistance.split(" ")
rdinfDistance = str(fl.ConvertDistanceToMeters(rddistance[0],rddistance[1]))
arcpy.AddMessage("Road Distance: " + str(rdinfDistance) + " meters...")
del rddistance

cExtent = tWorkspace+"\\aExtent_copy.shp" #Temporary copy of analysis extent for manipulation
## nTable = outWorkspace+"\\"+outTable #Final output table
nTable = "in_memory\\xxnTable" #Temporary output table
i = 1 # Number of iterations to run

constraintLayer = aExtent

### If constraint layer is specified, make copy and clip to aExtent
##if constraintLayer == "#":
##    constraintLayer = aExtent

currentDirectory = sys.path[0]

numHouses = fl.GetNumHouses(pointLayer, aExtent)
if numHouses == 0:
    arcpy.AddError("No structures found within analysis area.  Verify that all inputs are in the same projection...")
else:
    arcpy.AddMessage(str(numHouses) + " structures found within analysis area")


# numHouses = fl.GetNumHouses(pointLayer, aExtent)

##Calculate Total Area and store in variable 'tArea'##
arcpy.AddMessage("Caculating area of analysis extent...")
# Make copy of analysis extent layer to modify
##try:
##    arcpy.Copy_management(aExtent, cExtent, "ShapeFile")
##except:
##    a = arcpy.MakeFeatureLayer_management(aExtent, "aLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
##    arcpy.CopyFeatures_management(a, cExtent)

cExtent = arcpy.MakeFeatureLayer_management(aExtent, "cLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
    
# Add 'Area' Field and calculate Area
fl.CalcArea(cExtent)
tArea = fl.CalcTotalArea(arcpy.env.scratchWorkspace, "cLayer")
arcpy.AddMessage("Total Analysis Area = " + str(tArea))
## End Calculate Total Area ##

## Create New Table to store final results if it doesn't exist... ##

#If the nTable does not exist, create it.
arcpy.AddMessage("Creating output table " + nTable + "...")
if not arcpy.Exists(nTable):
    # Process: Create Table...
    fl.MakeTable(nTable.split("\\")[0], nTable.split("\\")[1])

# Get result 
arcpy.AddMessage("Begin simulation ...")
sim = fl.RunSimulation(i,tWorkspace, aExtent, constraintLayer, numHouses, minCoreSize, stinfDistance, rdinfDistance, nTable, tArea, pointLayer, roadLayer, False, "", "", "")
result = sim[0]

# Calculate density in acres from the solution
arcpy.AddMessage("Calculating density...")
density = float(tArea)/int(numHouses)
density = int(round(fl.ConvertSquareMetersToOther(density, outputUnits)))
result = int(round(result, 1))

# Get cores and houses of last iteration result
arcpy.AddMessage("Copying shapefiles from last iteration...")

# Join area calculations to core output table
arcpy.AddMessage("\tAppending area calculations to core shapefile table...")
arcpy.env.qualifiedFieldNames = "UNQUALIFIED"
arcpy.AddField_management(nTable, "CORE", "SHORT", "", "", "", "", "", "", "")
arcpy.CalculateField_management (nTable, "CORE", "1", "VB", "")
arcpy.Dissolve_management (tWorkspace + "\\xxcoresLayer", tWorkspace + "\\xxcoresLayer2", "CORE") 
arcpy.MakeFeatureLayer_management (tWorkspace + "\\xxcoresLayer2", "jLayer")
arcpy.AddJoin_management ("jLayer", "CORE", nTable, "CORE")
outShape2 = arcpy.CopyFeatures_management("jLayer", outShapefile)
arcpy.AddField_management(outShape2, "AREA", "DOUBLE", "", "", "", "", "", "", "")
arcpy.AddField_management(outShape2, "TOT_AREA", "DOUBLE", "", "", "", "", "", "", "")
arcpy.CalculateField_management (outShape2, "AREA", "!shape.area@" + outputUnits + "!", "PYTHON", "")
arcpy.AddField_management(outShape2, "PRCNT_TOT", "DOUBLE", "", "", "", "", "", "", "")
tArea = fl.ConvertSquareMetersToOther(tArea, outputUnits)
arcpy.CalculateField_management (outShape2, "TOT_AREA", str(tArea), "PYTHON", "")
arcpy.CalculateField_management (outShape2, "PRCNT_TOT", "(!AREA! / " + str(tArea) + ") * 100", "PYTHON", "")
arcpy.AddField_management(outShape2, "AREA_UNITS", "TEXT", "", "", "50")
arcpy.CalculateField_management(outShape2, "AREA_UNITS", '"' + outputUnits + '"', "PYTHON")

#Delete unwanted fields
arcpy.AddMessage('\tCleaning up output table...')
l = [u'Field1', u'OID_', u'SUM_SUM_AR', u'MAX_SUM_AR', u'MIN_SUM_AR', u'"RCNT_AREA', u'CORE_1']
for f in l:
    if arcpy.ListFields(outShape2, f):        
        arcpy.DeleteField_management(outShape2, f)

#Send output to ArcMap
arcpy.SetParameterAsText(7, outShapefile)
params = arcpy.GetParameterInfo()
    
# Clean up temporary workspace
try:
    fl.CleanFiles(sWorkspace)
    #fl.CleanFiles(tWorkspace)
    #arcpy.Delete_management(tWorkspace,"")
    arcpy.Delete_management(sWorkspace,"")
except:
    pass
try:
    arcpy.Delete_management("in_memory")
except:
    pass

### Display normal result in messagebox
##from Tkinter import *
##import tkMessageBox
##root = Tk()
##root.withdraw()
##tkMessageBox.showinfo(
##    title="Housing Density", \
##    message="Target: " + str(result) + "% of area conserved as core habitat\n Overall density: " + str(density) + " " + outputUnits + " per house")
##

    
