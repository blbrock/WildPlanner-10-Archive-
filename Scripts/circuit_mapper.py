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
# circuit_mapper.py
# Created on: Oct 21 2012 03:58:21 PM
# Author: Brent L. Brock, Craighead Institute 
# ---------------------------------------------------------------------------


# Import system modules
import sys, string, os, fileinput, arcpy
import local_params as lp
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
source = sys.argv[6]
sourceID = sys.argv[7]
outName = sys.argv[8]
output = sys.argv[9]
scriptpath = sys.path[0]

# Set the extent and cell size to link Layer
arcpy.env.extent = arcpy.Describe(linkLayer).extent
arcpy.env.cellSize = '"%s"' % arcpy.Describe(linkLayer).catalogPath
CellSize = str(arcpy.env.cellSize)

arcpy.env.cellSize = round(float(arcpy.env.cellSize), 12)
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
elif spPref == "NO-PREFERENCE":
    rFile = "none"

# Local variables...
lcProj = "xxlcproj"

# Make sure inputs are in same projection
linkProjName = arcpy.Describe(linkLayer).spatialreference.name
if not spPref == "NO-PREFERENCE":
    lcProjName = arcpy.Describe(lcLayer).spatialreference.name

    if not lcProjName == linkProjName:
        arcpy.AddWarning("Landcover and Linkage Layer projections do not match. Attempting to reproject " + lcLayer + "...")
        arcpy.ProjectRaster_management(lcLayer, lcProj, linkLayer)
        lcLayer = lcProj
        
# Create cost surface
arcpy.AddMessage("Creating cost surface raster...")
desc = arcpy.Describe(pntLayer) #.SpatialReference.LinearUnitName
r = desc.SpatialReference.LinearUnitName
if not r == "Meter":
    d = str(ConvertMetersToOther("908", r))
else:
    d = "908"

pntProjName = desc.spatialreference.name
# linkProjName = arcpy.Describe(linkLayer).spatialreference.name

if not pntProjName == linkProjName:
    arcpy.AddWarning("Structure Layer and Linkage Layer projections do not match. Attempting to reproject " + pntLayer + "...")
    arcpy.ProjectRaster_management(pntLayer, pntProj, linkLayer)
    pntLayer = pntProj
    
arcpy.AddMessage("\tCalculating structure density...")
pntDensity = arcpy.sa.PointDensity(pntLayer, "NONE", CellSize, "CIRCLE " + d + " MAP", "SQUARE_MILES")
arcpy.AddMessage("\tAssigning house density costs...")
pntDensity = arcpy.sa.ReclassByASCIIFile(pntDensity, dFile)

# Create cost surface
if not rFile == "none":
    arcpy.AddMessage("\tAssigning habitat costs...")
    reclass = arcpy.sa.ReclassByASCIIFile(lcLayer, rFile)

    arcpy.AddMessage("\tCombining habitat and house density costs...")
    xxplus = arcpy.sa.Plus(pntDensity, reclass)
    
else:
    arcpy.AddWarning("\tNo habitat preference specified... \n\tCalculating movement resistance from housing density only...")
    xxplus = pntDensity

xxCombine = arcpy.sa.Con(xxplus, "50", xxplus, "Value > 50")
costRaster = arcpy.sa.Con(linkLayer, "50", xxCombine, "Value = 0")
cell = arcpy.GetRasterProperties_management(costRaster, "CELLSIZEX").getOutput(0)

if not str(cell) == CellSize:
    arcpy.AddWarning("\tResampling Cost Surface to " + CellSize + " meter cells...")
    xxResample = arcpy.Resample_management(costRaster, 'in_memory\\xxResample', CellSize, "CUBIC").getOutput(0)
    costRaster = xxResample

# Extract cost raster by linkage layer to make sure cells outside analysis area are removed
########################### Need to FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!############################################
## Stating with arcgis 10.1 costRaster is being produced with negative values and 0 hidden outside the analysis area. These do not
## display and are removed with extract by mask in toolbox. But the line below does not work for some reason.
costRaster = arcpy.sa.ExtractByMask(costRaster, linkLayer)


## This is old rank linkages code
### Create source and destination polygons (NSEW)
##arcpy.AddMessage("Creating corridor destination polygons...")
##maxval = arcpy.GetRasterProperties_management (linkLayer, "MAXIMUM")
##xxblank = arcpy.sa.Reclassify(linkLayer, "VALUE", "0 " + str(maxval) + " 1", "DATA")
##arcpy.RasterToPolygon_conversion(xxblank, tWorkspace + "xxpoly", "SIMPLIFY", "VALUE")
##
### Added feature envelope to polygon to clean unwanted vertices from output
##arcpy.FeatureEnvelopeToPolygon_management(tWorkspace + "xxpoly.shp", tWorkspace + "xxpoly2.shp", "SINGLEPART")
##arcpy.FeatureToLine_management(tWorkspace + "xxpoly2.shp", tWorkspace + "xxline.shp", "", "ATTRIBUTES")
##arcpy.SplitLine_management(tWorkspace + "xxline.shp", tWorkspace + "xxsplit.shp")
##
##fc = arcpy.Buffer_analysis(tWorkspace + "xxsplit.shp", tWorkspace + "buf.shp", CellSize + " Meters", "RIGHT", "FLAT", "NONE", "")
##fc = arcpy.MakeFeatureLayer_management(fc, "pLayer")
##
##rows = arcpy.SearchCursor(fc) 
##
##for row in rows:
##    p = row.getValue("FID")
##    arcpy.MakeFeatureLayer_management(fc, "fLayer")
##    arcpy.SelectLayerByAttribute_management(fc, "NEW_SELECTION", "\"FID\" = " + str(p))
##    arcpy.MakeFeatureLayer_management(fc, "pLayer_" + str(p))
##    
##
##west = "pLayer_0"
##north = "pLayer_1"
##east = "pLayer_2"
##south = "pLayer_3"
##
##del row
##del rows
##del fc

### Run cost-distance analysis
##arcpy.AddMessage("Creating cost-distance rasters...")
### arcpy.ClearEnvironment("extent")
##if direction == "East-West" or direction == "East-West;North-South":
##    originEW = tWorkspace + "xxoriginEW.shp"
##    arcpy.Merge_management([east, west], originEW)
##    east_cst = arcpy.sa.CostDistance(east, costRaster)
##    west_cst = arcpy.sa.CostDistance(west, costRaster)
##    arcpy.AddMessage("Modeling east-west linkage zones...")
##    corTempEW = arcpy.sa.Corridor(east_cst, west_cst)
##    if not direction == "East-West;North-South":
##        corTemp = corTempEW
##        origin = originEW
##        
##if direction == "North-South" or direction == "East-West;North-South":
##    originNS = tWorkspace + "xxoriginNS.shp"
##    arcpy.Merge_management([north, south], originNS)
##    north_cst = arcpy.sa.CostDistance(north, costRaster)
##    south_cst = arcpy.sa.CostDistance(south, costRaster)
##    arcpy.AddMessage("Modeling north-south linkage zones...")
##    corTempNS = arcpy.sa.Corridor(north_cst, south_cst)
##    if not direction == "East-West;North-South":
##        corTemp = corTempNS
##        origin = originNS
##
##if direction == "East-West;North-South":
##    arcpy.AddMessage("Combining linkage zones...")
##    corTempEW.save(tWorkspace + "xxcorTempEW")
##    corTempEW = tWorkspace + "xxcorTempEW"
##    corTempNS.save(tWorkspace + "xxcorTempNS")
##    corTempNS = tWorkspace + "xxcorTempNS"
##    arcpy.AddWarning("Attempting to rescale " + corTempEW)
##    #Type = arcpy.GetRasterProperties (corTempEW, "MINIMUM")
##    Rescale(corTempEW, resclEW, tWorkspace)
##    arcpy.AddWarning("Attempting to rescale " + corTempNS)
##    Rescale(corTempNS, resclNS, tWorkspace)
##    arcpy.Mosaic_management(resclEW, resclNS, "MINIMUM")
##    corTemp = resclNS
##    origin = arcpy.Merge_management([originEW, originNS], tWorkspace + "xxorigin.shp")
##
### Create connected area mask
##arcpy.AddMessage("Creating connected area mask to remove isolated habitat patches...")
##sliceRaster = arcpy.sa.Slice(corTemp, 5, "NATURAL_BREAKS")
###
##sliceRaster.save(tWorkspace + "sliceRaster")
##movementAreas = arcpy.sa.SetNull(linkLayer, linkLayer, "Value = 0")
###
##movementAreas.save(tWorkspace + "movement")
##region = arcpy.sa.RegionGroup(movementAreas, "EIGHT")
###
##region.save(tWorkspace + "region")
##arcpy.RasterToPolygon_conversion(region, tWorkspace + "xxregion", "NO_SIMPLIFY", "VALUE")
##arcpy.MakeFeatureLayer_management(tWorkspace + "xxregion.shp", 'region_lyr')
##arcpy.SelectLayerByLocation_management('region_lyr', 'intersect', origin)
##arcpy.MakeFeatureLayer_management('region_lyr', 'connected')
##arcpy.AddMessage("Removing compromised areas from analysis area...")
##corNull = arcpy.sa.SetNull(linkLayer, corTemp, "Value = 0")
###
##corNull.save(tWorkspace + "corNull")
##arcpy.AddMessage("Removing isolated patches from analysis area...")
##corTemp2 = arcpy.sa.ExtractByMask(corNull, "connected")
###
##corTemp2.save(tWorkspace + "corTemp2")
##arcpy.AddMessage("Ranking linkage zones...")
##sliceRaster = arcpy.sa.Slice(corTemp2, 5, "NATURAL_BREAKS")
##conRaster = arcpy.sa.Con(linkLayer, "0", sliceRaster, "Value = 0")
##outIsNull = arcpy.sa.IsNull(conRaster)
##rankRaster = arcpy.sa.Con(outIsNull, "0", conRaster, "Value = 1")
##rankRaster.save(outName)

#Set up inputs for Circuitscape
arcpy.AddMessage("Performing Circuitscape Modeling...")

#Convert inputs to ASCII
arcpy.AddMessage("\tConverting Cost Surface to ASCII...")
# desc = arcpy.Describe(costRaster)

costName = outName + "_cost.asc"
arcpy.RasterToASCII_conversion(costRaster, costName)

desc = arcpy.Describe(source)
dType = desc.dataType
sourceName = outName + "_source.asc"

# Check if source patches are raster and convert if necessary
if dType == "ShapeFile" or dType == "FeatureLayer":
    arcpy.AddMessage("\tConverting Source Patches to ASCII...")
    arcpy.FeatureToRaster_conversion(source, sourceID, "in_memory\\xxsrcrast", CellSize) # Tried to use .getOutput(0) to set sourceRaster but get inconsistent results and errors.
    sourceRaster = "in_memory\\xxsrcrast"

    if arcpy.Describe(sourceRaster).SpatialReference.name == 'Unknown':
        sr = arcpy.Describe(source).SpatialReference
        arcpy.DefineProjection_management(sourceRaster, sr)
else:
    sourceRaster = source

# Make sure source patches are in correct projection
SrcPrj = arcpy.Describe(sourceRaster).spatialreference.name
if not SrcPrj == linkProjName:
    arcpy.AddMessage("\tReprojecting Source Patch Layer to " + linkProjName + "...")
    arcpy.ProjectRaster_management (sourceRaster, "in_memory\\xxsrcrast2", linkLayer) # Tried to use .getOutput(0) to set sourceRaster but get inconsistent results and errors.
    sourceRaster = "in_memory\\xxsrcrast2"
   
arcpy.RasterToASCII_conversion(sourceRaster, sourceName)


#Copy circuitscape.ini file to temp workspace and modify
arcpy.AddMessage("\tGenerating Circuitscape setup file...")
cs_ini = os.path.dirname(outName) + os.sep + "cs_ini.ini"
fileToSearch  = scriptpath + os.sep + "circuitscape.ini"
 
tempFile = open( cs_ini, 'w' ) 
 
for line in fileinput.input( fileToSearch ):
    tempFile.write(line.replace( "%NODES%", sourceName ) )
tempFile.close() 

tempFile = open( cs_ini, 'r+' )
                                 
for line in fileinput.input( cs_ini): 
    tempFile.write( line.replace( "%OUTPUT%", outName ) )
tempFile.close()

tempFile = open( cs_ini, 'r+' )
for line in fileinput.input( cs_ini):
    tempFile.write( line.replace( "%HABITAT%", costName ) )
tempFile.close()

# Run Circuitscape
arcpy.AddMessage("\tExecuting Circuitscape...")
cs_ini = cs_ini.replace("\\", "\\\\")

# Check location of cs_run.exe and change path if found in a different location
filename = "cs_run.exe"
n = False

if os.path.exists(lp.cs_path):
    arcpy.AddMessage('\tcs_run.exe found...')
else:
    arcpy.AddWarning("\tCircuitscape executable not found at default location. \n\tSearching C: drive for file...")
    print "I'm looking for the file..."
    for root, dirs, names in os.walk("c:\\"):
        if filename in names:
            n = os.path.join(root, filename)
            n = n.replace('\\', '\\\\')
            arcpy.AddMessage("\tCircuitscape executable found at: " + n)
            arcpy.AddMessage("\tUpdating local copy of script with new Circuitscape path...")
            print n

            # paramFile  = sys.argv[0]
            old = os.path.dirname(sys.argv[0]) + os.sep + "local_params.py"
            new = os.path.dirname(sys.argv[0]) + os.sep + "local_params.txt"
            f = open( old, 'r' ) 
            temp = open(new, 'w')
            for line in f:
                if "cs_path = " in line:
                    l = line.replace(line, "cs_path = '" + n) + "'\n"
                else:
                    l = line
                print l
                temp.write(l)
            f.close()

            temp.close()
            os.remove(old)
            os.rename(new, old)
            break
    else:
        arcpy.AddError('Circuitscape executable not found. Make sure Circuitscape is insalled \nor manually edit the cs_path variable in local_params.py')
        
# Execute Circuitscape
os.system('"' + lp.cs_path + '" ' + cs_ini)

inASCII = outName + "_cum_curmap.asc"
outRaster = outName + "_cum_curmap.img"
curTemp = os.path.dirname(outName) + os.sep + "curTemp.img"

arcpy.AddMessage("\tConverting Circuit Map to Raster...")
arcpy.ASCIIToRaster_conversion(inASCII, curTemp, "FLOAT")
if arcpy.Describe(curTemp).SpatialReference.name == 'Unknown':
    sr = arcpy.Describe(linkLayer).SpatialReference
    arcpy.DefineProjection_management(curTemp, sr)

# Process Circuitscape output and bring into map
arcpy.AddMessage("\tSetting Source Areas to NULL...")
outCon = arcpy.sa.Con(arcpy.sa.IsNull(sourceRaster),-9999, sourceRaster)

curRaster = arcpy.sa.SetNull(outCon, curTemp, "Value <> -9999")
curRaster = arcpy.sa.SetNull(curRaster, curRaster, "Value = 0")
# outRaster = curMapName
curRaster.save(outRaster)
desc = arcpy.Describe(linkLayer)
sr = desc.SpatialReference
arcpy.DefineProjection_management(outRaster, sr)
arcpy.AddMessage("\tComputing Final Raster Statistics...")
arcpy.CalculateStatistics_management(outRaster)

# delete resistance output file if it exists
name = outName + "_resistances_3columns.out"
if os.path.exists(name):
    os.remove(name)
n = name.split(".")
os.rename(n[0], name)
#Add outputs to display
arcpy.SetParameterAsText(8, outRaster)

# Clean up temporary workspace
os.remove(cs_ini)
arcpy.Delete_management(curTemp)
arcpy.Delete_management("in_memory")
