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
direction = sys.argv[6]
outRaster = sys.argv[7]

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

# Set the extent to link Layer
# arcpy.Extent = arcpy.Describe(linkLayer).extent

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
east_cst = tWorkspace + os.sep + "xxeastcst"
west_cst = tWorkspace + os.sep + "xxwestcst"
north_cst = tWorkspace + os.sep + "xxnorthcst"
south_cst = tWorkspace + os.sep + "xxsouthcst"
corTempNS = tWorkspace + os.sep + "xxcorns"
corTempEW = tWorkspace + os.sep + "xxcorew"
resclNS = tWorkspace + os.sep + "xxcornsr"
resclEW = tWorkspace + os.sep + "xxcorewr"
corTemp = tWorkspace + os.sep + "xxcor"
corTemp2 = tWorkspace + os.sep + "xxcor2"
lcProj = tWorkspace + os.sep + "xxlcproj"

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
costRaster.save(tWorkspace + os.sep + "costraster1") 

# Create source and destination polygons (NSEW)
arcpy.AddMessage("Creating corridor destination polygons...")
maxval = arcpy.GetRasterProperties_management (linkLayer, "MAXIMUM")
xxblank = arcpy.sa.Reclassify(linkLayer, "VALUE", "0 " + str(maxval) + " 1", "DATA")
arcpy.RasterToPolygon_conversion(xxblank, tWorkspace + os.sep + "xxpoly", "SIMPLIFY", "VALUE")

# Added feature envelope to polygon to clean unwanted vertices from output
arcpy.FeatureEnvelopeToPolygon_management(tWorkspace + os.sep + "xxpoly.shp", tWorkspace + os.sep + "xxpoly2.shp", "SINGLEPART")
arcpy.FeatureToLine_management(tWorkspace + os.sep + "xxpoly2.shp", tWorkspace + os.sep + "xxline.shp", "", "ATTRIBUTES")
arcpy.SplitLine_management(tWorkspace + os.sep + "xxline.shp", tWorkspace + os.sep + "xxsplit.shp")

fc = arcpy.Buffer_analysis(tWorkspace + os.sep + "xxsplit.shp", tWorkspace + os.sep + "buf.shp", CellSize + " Meters", "RIGHT", "FLAT", "NONE", "")
fc = arcpy.MakeFeatureLayer_management(fc, "pLayer")

rows = arcpy.SearchCursor(fc) 

for row in rows:
    p = row.getValue("FID")
    arcpy.MakeFeatureLayer_management(fc, "fLayer")
    arcpy.SelectLayerByAttribute_management(fc, "NEW_SELECTION", "\"FID\" = " + str(p))
    arcpy.MakeFeatureLayer_management(fc, "pLayer_" + str(p))
    

west = "pLayer_0"
north = "pLayer_1"
east = "pLayer_2"
south = "pLayer_3"

del row
del rows
del fc

# Run cost-distance analysis
arcpy.AddMessage("Creating cost-distance rasters...")
# arcpy.ClearEnvironment("extent")
if direction == "East-West" or direction == "East-West;North-South":
    originEW = tWorkspace + os.sep + "xxoriginEW.shp"
    arcpy.Merge_management([east, west], originEW)
    east_cst = arcpy.sa.CostDistance(east, costRaster)
    west_cst = arcpy.sa.CostDistance(west, costRaster)
    arcpy.AddMessage("Modeling east-west linkage zones...")
    corTempEW = arcpy.sa.Corridor(east_cst, west_cst)
    if not direction == "East-West;North-South":
        corTemp = corTempEW
        origin = originEW
        
if direction == "North-South" or direction == "East-West;North-South":
    originNS = tWorkspace + os.sep + "xxoriginNS.shp"
    arcpy.Merge_management([north, south], originNS)
    north_cst = arcpy.sa.CostDistance(north, costRaster)
    south_cst = arcpy.sa.CostDistance(south, costRaster)
    arcpy.AddMessage("Modeling north-south linkage zones...")
    corTempNS = arcpy.sa.Corridor(north_cst, south_cst)
    if not direction == "East-West;North-South":
        corTemp = corTempNS
        origin = originNS

if direction == "East-West;North-South":
    arcpy.AddMessage("Combining linkage zones...")
    corTempEW.save(tWorkspace + os.sep + "xxcorTempEW")
    corTempEW = tWorkspace + os.sep + "xxcorTempEW"
    corTempNS.save(tWorkspace + os.sep + "xxcorTempNS")
    corTempNS = tWorkspace + os.sep + "xxcorTempNS"
    arcpy.AddWarning("Attempting to rescale " + corTempEW)
    #Type = arcpy.GetRasterProperties (corTempEW, "MINIMUM")
    Rescale(corTempEW, resclEW, tWorkspace)
    arcpy.AddWarning("Attempting to rescale " + corTempNS)
    Rescale(corTempNS, resclNS, tWorkspace)
    arcpy.Mosaic_management(resclEW, resclNS, "MINIMUM")
    corTemp = resclNS
    origin = arcpy.Merge_management([originEW, originNS], tWorkspace + os.sep + "xxorigin.shp")

# Create connected area mask
arcpy.AddMessage("Creating connected area mask to remove isolated habitat patches...")
sliceRaster = arcpy.sa.Slice(corTemp, 5, "NATURAL_BREAKS")
#
sliceRaster.save(tWorkspace + os.sep + "sliceRaster")
movementAreas = arcpy.sa.SetNull(linkLayer, linkLayer, "Value = 0")
#
movementAreas.save(tWorkspace + os.sep + "movement")
region = arcpy.sa.RegionGroup(movementAreas, "EIGHT")
#
region.save(tWorkspace + os.sep + "region")
arcpy.RasterToPolygon_conversion(region, tWorkspace + os.sep + "xxregion", "NO_SIMPLIFY", "VALUE")
arcpy.MakeFeatureLayer_management(tWorkspace + os.sep + "xxregion.shp", 'region_lyr')
arcpy.SelectLayerByLocation_management('region_lyr', 'intersect', origin)
arcpy.MakeFeatureLayer_management('region_lyr', 'connected')
arcpy.AddMessage("Removing compromised areas from analysis area...")
corNull = arcpy.sa.SetNull(linkLayer, corTemp, "Value = 0")
#
corNull.save(tWorkspace + os.sep + "corNull")
arcpy.AddMessage("Removing isolated patches from analysis area...")
corTemp2 = arcpy.sa.ExtractByMask(corNull, "connected")
#
corTemp2.save(tWorkspace + os.sep + "corTemp2")
arcpy.AddMessage("Ranking linkage zones...")
sliceRaster = arcpy.sa.Slice(corTemp2, 5, "NATURAL_BREAKS")
conRaster = arcpy.sa.Con(linkLayer, "0", sliceRaster, "Value = 0")
outIsNull = arcpy.sa.IsNull(conRaster)
rankRaster = arcpy.sa.Con(outIsNull, "0", conRaster, "Value = 1")
rankRaster.save(outRaster)


#Add outputs to display
arcpy.SetParameterAsText(6, outRaster)

# Clean up temporary workspace
arcpy.Delete_management("pLayer")
arcpy.Delete_management(south)
arcpy.Delete_management(west)
arcpy.Delete_management(north)
arcpy.Delete_management(east)
CleanFiles(tWorkspace)
arcpy.Delete_management(tWorkspace,"")










    
    






