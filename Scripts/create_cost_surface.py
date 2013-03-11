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
# corridor_mapper.py
# Created on: Thu Apr 14 2011 03:58:21 PM
# Author: Brent L. Brock, Craighead Institute
# Usage: corridor_mapper <vshdLayer> <lcLayer> 
# ---------------------------------------------------------------------------


# Import system modules
import sys, string, os, arcpy
from functionlib import CreateTempWorkspace, CleanFiles, ConvertMetersToOther, Rescale

# Create the Geoprocessor object

# Set the necessary product code
arcpy.SetProduct("ArcInfo")
arcpy.overwriteoutput = 1

arcpy.AddMessage("\tChecking availability of spatial analyst license...")  
if arcpy.CheckExtension("spatial") == "Available":
    arcpy.CheckOutExtension("spatial")         
else:
    arcpy.AddWarning("\tSpatial Analyst is needed to run this tool.  Processing aborted...")
    raise "LicenseError"    

# Script arguments...
pntLayer = sys.argv[1]
lcLayer = sys.argv[2]
lcSource = sys.argv[3]
spPref = sys.argv[4]
linkLayer = sys.argv[5]
outRaster = sys.argv[6]

# Set up temporary workspace
scratchWS = arcpy.env.scratchWorkspace
if scratchWS:
    tWorkspace = CreateTempWorkspace(scratchWS)
else:
    # Set output workspace
    r = outRaster.rsplit(os.sep,1)
    outWorkspace = r[0]
    arcpy.Workspace = outWorkspace
    tWorkspace = CreateTempWorkspace(outWorkspace)

arcpy.env.extent = arcpy.Describe(linkLayer).catalogPath
arcpy.env.cellSize = arcpy.Describe(linkLayer).catalogPath
CellSize = str(arcpy.env.cellSize)

# Set path for density reclass table
dFile = sys.path[0] + os.sep + "ReclassTables" + os.sep + "density.txt"

# Locate path for correct reclass tables
lcSource = lcSource.replace("'", "")
if lcSource == "Montana Landcover;National Landcover Dataset (NLCD)":
    lcSource = lcSource.split(";")[0]
    
if lcSource == "Montana Landcover":
    fPath = sys.path[0] + os.sep + "ReclassTables" + os.sep + "MT_landcov"
    
elif lcSource == "National Landcover Dataset (NLCD)":
    fPath = sys.path[0] + os.sep + "ReclassTables" + os.sep + "NLCD"

# Find the correct reclass table
if spPref == "FOREST-SHRUB":
    rFile = fPath + os.sep + "ForestClasses.txt"
elif spPref == "GRASSLAND-OPEN":
    rFile = fPath + os.sep + "ForestClasses.txt"
elif spPref == "WETLAND-MARSH":
    rFile = fPath + os.sep + "WetlandClasses.txt"

# Local variables...
reclass = tWorkspace + os.sep + "xxrcls"
##east_cst = tWorkspace + os.sep + "xxeastcst"
##west_cst = tWorkspace + os.sep + "xxwestcst"
##north_cst = tWorkspace + os.sep + "xxnorthcst"
##south_cst = tWorkspace + os.sep + "xxsouthcst"
##corTempNS = tWorkspace + os.sep + "xxcorns"
##corTempEW = tWorkspace + os.sep + "xxcorew"
##resclNS = tWorkspace + os.sep + "xxcornsr"
##resclEW = tWorkspace + os.sep + "xxcorewr"
##corTemp = tWorkspace + os.sep + "xxcor"
##corTemp2 = tWorkspace + os.sep + "xxcor2"
##lcProj = tWorkspace + os.sep + "xxlcproj"

# Make sure inputs are in same projection

lcProjName = arcpy.Describe(lcLayer).spatialreference.name
linkProjName = arcpy.Describe(linkLayer).spatialreference.name

if not lcProjName == linkProjName:
    arcpy.AddWarning("Landcover and Linkage Layer projections do not match. Attempting to reproject " + lcLayer + "...")
    arcpy.ProjectRaster_management(lcLayer, lcProj, linkLayer)
    lcLayer = lcProj

# Create cost surface
arcpy.AddMessage("Creating cost surface raster...")
r = arcpy.Describe(pntLayer).SpatialReference.LinearUnitName
if not r == "Meter":
    d = str(ConvertMetersToOther("908", r))
else:
    d = "908"
    
arcpy.AddMessage("\tCalculating structure density...")
pntDensity = arcpy.sa.PointDensity(pntLayer, "NONE", "30", "CIRCLE " + d + " MAP", "SQUARE_MILES")
arcpy.AddMessage("\tAssigning house density costs...")
pntDensity = arcpy.sa.ReclassByASCIIFile(pntDensity, dFile)
#
pntDensity.save(tWorkspace + os.sep + "xxdensity")
# Create cost surface
arcpy.AddMessage("\tAssigning habitat costs...")
reclass = arcpy.sa.ReclassByASCIIFile(lcLayer, rFile)
#
reclass.save(tWorkspace + os.sep + "xxrcls")


arcpy.AddMessage("\tCombining habitat and house density costs...")
xxplus = arcpy.sa.Plus(pntDensity, reclass)
#
reclass.save(tWorkspace + os.sep + "xxplus")

xxCombine = arcpy.sa.Con(xxplus, "50", xxplus, "Value > 50")

##maxval = arcpy.GetRasterProperties_management (xxplus, "MAXIMUM")
##maxval = float(str(maxval))    
##maxcost = maxval * 10
costRaster = arcpy.sa.Con(linkLayer, "50", xxplus, "Value = 0")
costRaster.save(outRaster) 


#Add outputs to display
arcpy.SetParameterAsText(5, outRaster)

# Clean up temporary workspace
CleanFiles(tWorkspace)
arcpy.Delete_management(tWorkspace,"")










    
    






