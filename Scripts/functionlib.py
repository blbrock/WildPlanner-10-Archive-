# ------------------------------------------------------------------------------
#
# Copyright 2008, 2009, 2010, 2011, 2012, 2013 Brent L. Brock and the Craighead Institute
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
# ---------------------------------------------------------------------------
# functionlib.py
#
# Created December 5 2008 
#
# Written by Brent L. Brock, Landscape Ecologist, Craighead Environmental Research Institute
#
# A library of functions required to run housing density simulations for wildlife conservation
#
#These functions are intended for use with calculate_density
# 
# ---------------------------------------------------------------------------

# Import system modules
import os, random, arcgisscripting, sys, arcpy
from math import pi
gp = arcgisscripting.create()
# Set the necessary product code
gp.SetProduct("ArcInfo")
gp.overwriteoutput = 1

## Create 'Area' field for polygon shapes
## Note: in arcgis 9.2 there appears to be a problem with using virtual layers as input
## to this function.  The listfields method will return false regardless of whether
## an 'Area' field actually exists.  If the field does exist, an error will be thrown
## when the function tries to create a new area field.  Running this same script in arcgis 9.3
## but using an arcgis 9.2 geoprocessor seems to work.
def CalcArea(inputFC):
        arcpy.AddMessage("\tCalculating area of polygon cores...")
        if not arcpy.ListFields(inputFC,"AREA"):
                arcpy.AddField_management(inputFC, "Area", "DOUBLE", "18", "4")

        # Calculate 'Area' field to shape area
        expression = "float(!shape.area!)"
        arcpy.CalculateField_management(inputFC, "Area", expression, "PYTHON")

## Calculate the mean distance between random points relative to a minimum distance threshold
def CalcCorridorDistance(pointLayer, workspace):
        arcpy.AddMessage("\tCalculating distances between structures...")

        if workspace == '':
                workspace = arcpy.env.Workspace

        elif workspace == "in_memory":
                sumTable = workspace + os.sep + "Core_summary"
                iTable = workspace + os.sep + "Iteration_summary"
        else:
                sumTable = workspace + os.sep + "Core_summary.dbf"
                iTable = workspace + os.sep + "Iteration_summary.dbf"
        cLayer = workspace + os.sep + "coresLayer"
        
        whereClause = '"NEAR_DIST" <> -1'
        
        # Calculate distance between points
        arcpy.Near_analysis (pointLayer, pointLayer) 
   
        # Process: Make Feature Layer ...
        arcpy.MakeFeatureLayer_management(pointLayer, cLayer, "", "", "Input_FID Input_FID VISIBLE NONE")
        
        # Process: Select Near Distance <> -1 (distance between a point and itself)...
        arcpy.SelectLayerByAttribute_management(cLayer, "NEW_SELECTION", whereClause)
    
        # Process: Summary Statistics on selected set...
        arcpy.Statistics_analysis(cLayer, iTable, "NEAR_DIST MEAN")

        arcpy.Delete_management(cLayer,"")

        return iTable #send iteration table back for further processing

# Calculates the mean of all records for the input table ('inTable') and specified field ('theField')
def CalcMeans(inTable, workspace, theField, numHouses):
        arcpy.AddMessage("\tCalculating means for " + theField + " in table: " + inTable + "...")
        if workspace == '':
                workspace = "in_memory"

        elif workspace == "in_memory":
                mTable = workspace + os.sep + "mean_value"
                
        else:
                mTable = workspace + os.sep + "mean_value.dbf"
        
        
        # Get value from summary table
        
        if not numHouses == "":
                arcpy.Statistics_analysis(inTable, mTable, theField + " MEAN", "HOUSES")
                theField = "MEAN_" + theField
                theField = arcpy.ValidateFieldName(theField, workspace)
                rows = arcpy.SearchCursor(mTable, "HOUSES = " + str(numHouses))
        else:
                arcpy.Statistics_analysis(inTable, mTable, theField + " MEAN")
                # Get value from summary table
                theField = "MEAN_" + theField
                theField = arcpy.ValidateFieldName(theField, workspace)
                rows = arcpy.SearchCursor(mTable)
                
        arcpy.AddMessage("\tGetting value for " + theField)
        row = rows.next()
        theMean = row.getValue(theField)
        
##      Delete the row and cursor
        del row, rows
        return theMean

## Calculate total area from a table containing an 'Area' field
def CalcTotalArea(workspace, extentLayer):
        arcpy.AddMessage("\tCalculating total area of " + extentLayer + "...")
        if workspace == "":
                workspace = arcpy.env.Workspace
        aTable = workspace + os.sep + "aTable.dbf"
        
        # Process: Summary Statistics...
        arcpy.Statistics_analysis(extentLayer, aTable, "AREA SUM", "")
        # Get value from summary table
        rows = arcpy.SearchCursor(aTable)
        row = rows.next()
        tArea = row.getValue("SUM_AREA")
        # Delete the row and cursor
        del row, rows
        try:
                arcpy.Delete_management(aTable,"")
        except:
                arcpy.AddWarning("\tUnable to delete " + aTable + "...")

        return tArea

## Clean up tables and feature classes in temporary workspace
def CleanFiles(workspace):
        gp.AddMessage("\tDeleting temporary files...")
        gp.Workspace = workspace
        # Get a list of the tables in the workspace and delete.
        tables = gp.ListTables()
##  This sections includes syntax for argis processor 9.3 but does not work with older versions
##        for t in tables:
        t = tables.Next()
        while t:
                gp.Delete_management(t,"")
                t = tables.Next()

        # Get a list of the feature classes in the workspace and delete.
        fcs = gp.ListFeatureClasses()
        fc = fcs.Next()
## The following line is version 9.3 syntax
##        for fc in fcs:
        while fc:
                gp.Delete_management(fc,"")
                fc = fcs.Next()

# Get a list of the rasters in the workspace and delete.
        rss = gp.ListRasters()
        rs = rss.Next()
## The following line is version 9.3 syntax
##        for rs in eas:
        while rs:
                gp.Delete_management(rs,"")
                rs = rss.Next()

## Convert area units into square meters
def ConvertAreaToMeters(areaValue, units):
    areaValue = float(areaValue)
    units = units.replace(" ", "")
    
    if units == 'SquareMeters':
        pass
    elif units == 'Acres':
        areaValue = areaValue * 4047
    elif units == 'Ares':
        areaValue = areaValue * 100
    elif units == 'Hectares':
        areaValue = areaValue * 10000
    elif units == 'SquareCentimeters':
        areaValue = areaValue * 0.0001
    elif units == 'SquareDecimeters':
        areaValue = areaValue * 0.01
    elif units == 'SquareFeet':
        areaValue = areaValue * 0.0929
    elif units == 'SquareInches':
        areaValue = areaValue * 0.0006452
    elif units == 'SquareKilometers':
        areaValue = areaValue * 1000000
    elif units == 'SquareMillimeters':
        areaValue = areaValue * 0.000001
    elif units == 'SquareMiles':
        areaValue = areaValue * 2590000
    elif units == 'SquareYards':
        areaValue = areaValue * 0.8361
    else:
        gp.AddError("\tUnknown units " + units + " could not convert to square meters")
    return areaValue

## Convert square meters to other area units
def ConvertSquareMetersToOther(areaValue, units):
    areaValue = float(areaValue)
    units = units.replace(" ", "")
    
    if units == 'Acres':
        areaValue = areaValue * 0.000247105381467165
    elif units == 'Ares':
        areaValue = areaValue * 0.01
    elif units == 'Hectares':
        areaValue = areaValue * 0.0001
    elif units == 'SquareCentimeters':
        areaValue = areaValue * 10000
    elif units == 'SquareDecimeters':
        areaValue = areaValue * 100
    elif units == 'SquareFeet':
        areaValue = areaValue * 10.76
    elif units == 'SquareInches':
        areaValue = areaValue * 1550
    elif units == 'SquareKilometers':
        areaValue = areaValue * 0.000001
    elif units == 'SquareMeters':
        areaValue = areaValue * 1
    elif units == 'SquareMillimeters':
        areaValue = areaValue * 1000000
    elif units == 'SquareMiles':
        areaValue = areaValue * 0.0000003861
    elif units == 'SquareYards':
        areaValue = areaValue * 1.196
    else:
        gp.AddError("\tUnknown units " + units + " could not convert from square meters")
    return areaValue

## Convert distance units to meters
def ConvertDistanceToMeters(distanceValue, units):
    distanceValue = float(distanceValue)
    units = units.replace(" ", "")
    
    if units == 'Centimeters':
        distanceValue = distanceValue * 0.01
    elif units == 'DecimalDegrees':
        gp.AddError("\tCannot convert from Decimal Degrees.  Choose different units.")
    elif units == 'Feet':
        distanceValue = distanceValue * 0.3048
    elif units == 'Inches':
        distanceValue = distanceValue * 0.0254
    elif units == 'Kilometers':
        distanceValue = distanceValue * 1000
    elif units == 'Meters':
        distanceValue = distanceValue * 1
    elif units == 'Millimeters':
        distanceValue = distanceValue * 0.001
    elif units == 'Miles':
        distanceValue = distanceValue * 1609
    elif units == 'NauticalMiles':
        distanceValue = distanceValue * 1852
    elif units == 'Points':
        distanceValue = distanceValue * 0.0003514
    elif units == 'Yards':
        distanceValue = distanceValue * 0.9144
    else:
        gp.AddError("\tUnknown units " + units + " could not convert to meters")
    return distanceValue

## Convert meters to other linear units
def ConvertMetersToOther(distanceValue, units):
    distanceValue = float(distanceValue)
    units = units.replace(" ", "")
    
    if units == 'Centimeters':
        distanceValue = distanceValue * 100
    elif units == 'DecimalDegrees':
        distanceValue = 6378137.0 * math.pi * distanceValue / 180.0
    elif units == 'Feet':
        distanceValue = distanceValue * 3.281
    elif units == 'Inches':
        distanceValue = distanceValue * 39.37
    elif units == 'Kilometers':
        distanceValue = distanceValue * 0.001
    elif units == 'Meters':
        distanceValue = distanceValue * 1
    elif units == 'Millimeters':
        distanceValue = distanceValue * 1000
    elif units == 'Miles':
        distanceValue = distanceValue * 0.0006214
    elif units == 'NauticalMiles':
        distanceValue = distanceValue * 0.00054
    elif units == 'Points':
        distanceValue = distanceValue * 2845
    elif units == 'Yards':
        distanceValue = distanceValue * 1.094
    else:
        gp.AddError("\tUnknown units " + units + " could not convert to meters")
    return distanceValue

# Create road cost surface
def CreateRoadCost (inDEM, inSlope, inRoads, outWorkspace):
        gp.AddMessage("Creating road cost surface...")
              
        ##    # Set up workspaces
        ##    if outWorkspace == '' or outWorkspace == '#':
        ##        r = outRdCost.rsplit(os.sep,1)
        ##        outWorkspace = r[0]
        ##    else:
        ##        outWorkspace = dName[0]

        # Create temporary workspaces
        tWorkspace = CreateTempWorkspace(outWorkspace)

        # Specify output filenames
        dName = inRoads.rsplit(os.sep,1)
        c = len(dName) - 1
        dName = dName[c]
        dNamePath = outWorkspace + os.sep + dName[:7]
        outRdCost = dNamePath + "_rdcst"
        outBacklink = dNamePath + "_bklnk"

        # Local variables...
        outContour = tWorkspace + os.sep + "outContour.shp"
        Contour_1 = tWorkspace + os.sep + "cRast"
        Reclass_Cont2 = tWorkspace + os.sep + "Reclass_Cont2"

        # Set cell size to input DEM
        CstRaster = tWorkspace + os.sep + "CstRaster"
        cellSize = str(gp.GetRasterProperties (inDEM, "CELLSIZEX"))
        
## UNCOMMENT THIS SECTION TO SPEED PROCESSING ##       
##        if cellSize < 90:
##                cellSize = "90"
##        else:
##                cellSize = str(cellSize)

        # Process: Contour...
        gp.AddMessage("\tGenerating contours from " + inDEM + "...")
        gp.Contour_sa(inDEM, outContour, "50", "0", "1")

        # Process: Add Field...
        gp.AddMessage("\tCalculating Value Field...")
        gp.AddField_management(outContour, "VALUE", "SHORT", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")

        # Process: Calculate Field...
        gp.CalculateField_management(outContour, "VALUE", "1", "VB", "")

        # Process: Polyline to Raster...
        gp.AddMessage("\tConverting contours to raster...")
        gp.PolylineToRaster_conversion(outContour, "VALUE", Contour_1, "MAXIMUM_LENGTH", "NONE", cellSize)

        # Process: Reclassify...
        gp.Reclassify_sa(Contour_1, "VALUE", "1 1;NODATA 0", Reclass_Cont2, "DATA")
        
        # Generate slope raster if it does not exist      

        if not gp.Exists(inSlope):
                gp.AddMessage("\tCalculating slope...")
                if inSlope == "" or inSlope == "#":
                        inSlope = dNamePath + "_slope"
                gp.Slope_sa(inDEM, inSlope, "PERCENT_RISE","")

        # Process: Con...
        gp.Con_sa(Reclass_Cont2, "1", CstRaster, inSlope, "VALUE = 1")

        # Process: Cost Distance...
        gp.AddMessage("\tGenerating cost-distance raster...")
        gp.CostDistance_sa(inRoads, CstRaster, outRdCost, "", outBacklink)

        # Clean up temporary workspace
        CleanFiles(tWorkspace)
        gp.Delete_management(tWorkspace,"")

        return outRdCost, outBacklink

## Create a temporary workspace
def CreateTempWorkspace(outWorkspace):
        # Set up temporary workspace
        if outWorkspace == "":
                outWorkspace = gp.Workspace
                desc = gp.Describe(outWorkspace)
                if not desc.workspaceType == 'FileSystem':
                        gp.AddError("\tCannot create temporary workspace in a geodatabase.\n\tChange output workspace environment settings and try again...")
        if outWorkspace == "" and gp.Workspace == None:
                gp.AddError("\tOutput Workspace is not defined...")
        gp.Workspace = outWorkspace
        if gp.ScratchWorkspace == None:
            gp.ScratchWorkspace = outWorkspace
        desc = gp.Describe(gp.ScratchWorkspace)
        if not desc.workspaceType == 'FileSystem':
                gp.ScratchWorkspace = outWorkspace
        WorkspaceName = os.sep + "xx" + str(random.randrange(0,1000000))
        Workspace = gp.ScratchWorkspace + WorkspaceName
        gp.AddMessage("\tCreating temporary workspace \n\t" + Workspace + "...")
        if not gp.Exists(Workspace):
                gp.CreateArcInfoWorkspace(gp.ScratchWorkspace, WorkspaceName)
        else:
                gp.AddWarning("\t" + Workspace + " already exists...")
        return Workspace

## Simulate undisturbed habitat areas created by a given input point pattern.
## Results are clipped to boundaries of polygons in aExtent.
def GenerateCores (pointLayer, workspace, aExtent, stinfDistance, rdinfDistance, lineLayer):
        arcpy.AddMessage("Generating habitat cores...")
        #Check out necessary extensions...
        arcpy.CheckOutExtension("spatial")         

        #Set workspace
        if workspace == '':
                workspace = arcpy.Workspace
                
##        #Create local variables
##        pBuff = workspace + os.sep + "pointBuffer.shp"
##        lBuff = workspace + os.sep + "lineBuffer.shp"
##        polyLayer = workspace + os.sep + "polyLayer.shp"
##        #lineLayer = workspace + os.sep + "lineLayer.shp"
##        eraseLayer1 = workspace + os.sep + "eraseLayer1.shp"
##        eraseLayer2 = workspace + os.sep + "eraseLayer2.shp"
##        clip = workspace + os.sep + "clipLayer.shp"
##        habitat = workspace + os.sep + "habitatPatches.shp"

        #Create local variables
        pBuff = workspace + os.sep + "pointBuffer"
        lBuff = workspace + os.sep + "lineBuffer"
        polyLayer = workspace + os.sep + "polyLayer"
        #lineLayer = workspace + os.sep + "lineLayer.shp"
        eraseLayer1 = workspace + os.sep + "eraseLayer1"
        eraseLayer2 = workspace + os.sep + "eraseLayer2"
        clip = workspace + os.sep + "clipLayer"
        habitat = workspace + os.sep + "habitatPatches"
        
        empty = False
        
        # Set the analysis extent to aExtent
        desc = arcpy.Describe(aExtent)
        extent = desc.Extent
        arcpy.Extent = extent
        polyLayer = aExtent
        pLayer = arcpy.MakeFeatureLayer_management(pointLayer, "pLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
        pointLayer = arcpy.SelectLayerByLocation_management(pLayer, "COMPLETELY_WITHIN", aExtent, "", "NEW_SELECTION")

## NOTE: Vector processing eliminated to improve performance.
##        try:
##        
##                # Process: Buffer points by disturbance distance...
##                arcpy.AddMessage("\tBuffering structures by disturbance distance...")
##                arcpy.Buffer_analysis(pointLayer, pBuff, buffDistance, "FULL", "ROUND", "NONE", "")
##                # Process: Buffer lines (roads) by disturbance distance ...
##                arcpy.AddMessage("\tBuffering roads by disturbance distance...")
##                arcpy.Buffer_analysis(lineLayer, lBuff, buffDistance, "FULL", "ROUND", "NONE", "")
##        
####                try:
##                # Process: Erase point (houses) disturbance buffers from thiessen polygons ...
##                arcpy.AddMessage("\tRemoving point buffers from analysis extent...")
##                arcpy.Erase_analysis(polyLayer, pBuff, eraseLayer1, "")
##                arcpy.Delete_management(pBuff,"")
##                # Process: Erase line buffers from thiessen polygons (note: polygons have now been erased for both houses and roads) ...
##                arcpy.AddMessage("\tRemoving road buffers from analysis extent...")
##                arcpy.Erase_analysis(eraseLayer1, lBuff, eraseLayer2, "")
##                arcpy.Delete_management(lBuff,"")
##
####                except:
####                        # erase fails, assume extent is empty
####                        empty = True
##        except:
            # if buffer or erase fails, try running processing in raster format.
##        arcpy.AddWarning("\tVector processing failed due to memory limitations... \n\tAttempting to continue processing using raster analysis...")

##        EucDist1 = workspace + os.sep + "EucDist1"           
##        SetNull1 = workspace + os.sep + "SetNull1"             
##        EucDist2 = workspace + os.sep + "EucDist2"              
##        SetNull2 = workspace + os.sep + "SetNull2"
        arcpy.AddMessage("\tCalculating euclidian distance from structures...")
        EucDist1 = arcpy.sa.EucDistance(pointLayer, "", "30", "")
        arcpy.AddMessage("\tSetting areas within " + str(stinfDistance) + " meters of structures to NULL...")
        SetNull1 = arcpy.sa.SetNull(EucDist1, "1", "Value < " + str(stinfDistance))
        arcpy.AddMessage("\tCalculating euclidian distance from roads...")
        EucDist2 = arcpy.sa.EucDistance(lineLayer, "", "30", "")
        arcpy.AddMessage("\tSetting areas within " + str(rdinfDistance) + " meters of roads to NULL...")
        SetNull2 = arcpy.sa.SetNull(EucDist2, SetNull1, "Value < " + str(rdinfDistance))

        del pointLayer
        del pLayer
        try:
            arcpy.AddMessage("\tConverting raster to polygons...")
            arcpy.RasterToPolygon_conversion(SetNull2, eraseLayer2, "NO_SIMPLIFY")
        except:
            empty = True

        try: 
                # Process: Clip cores by analysis extent ...
                arcpy.Clip_analysis(eraseLayer2, aExtent, clip, "")
        except:
                empty = True
                
        if not empty:
                
                # Make sure features are single part for proper area calculations
                arcpy.AddMessage("\tGenerating single part features...")
                arcpy.MultipartToSinglepart_management (clip, habitat)
                
                return habitat
        else:
                #if extent is empty return false
                return False

# Generate Roads
def GenRoads (inPoints, RdCost, backlink, outRoads):

        gp.AddMessage("Generating simulated road network to structures...")

        # Set up workspaces
        if outRoads == "" or outRoads == "#":
                desc = gp.Describe(RdCost)
                outWorkspace = desc.Path
                dName = desc.BaseName
                dNamePath = outWorkspace + os.sep + dName[:7]
                outRoads = dNamePath + "_sim_roads"
        
        r = outRoads.rsplit(os.sep,1)
        outWorkspace = r[0]

        # Create temporary workspaces
        tWorkspace = CreateTempWorkspace(outWorkspace)

        # Local variables...
        Output_raster = tWorkspace + os.sep + "Costrast"

        # Process: Cost Path...
        gp.AddMessage("\tGenerating least-cost paths for roads...")
        gp.CostPath_sa(inPoints, RdCost, backlink, Output_raster, "", "")

        # Set Workspace for permanent output
        gp.Workspace = outWorkspace

        # Process: Raster to Polyline...
        gp.AddMessage("\tConverting raster to polyline featureclass...")
        outRoads = gp.RasterToPolyline_conversion(Output_raster, outRoads, "ZERO", "0", "SIMPLIFY", "VALUE")


        # Clean up temporary workspace
        CleanFiles(tWorkspace)
        gp.Delete_management(tWorkspace,"")

        return outRoads

## Get the number of structure points from a point layer
def GetNumHouses(pointLayer, aExtent):
        gp.AddMessage("\tCalculating the number of structures within analysis extent...")
        gp.MakeFeatureLayer_management(pointLayer, "pLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
        gp.SelectLayerByLocation("pLayer", "COMPLETELY_WITHIN", aExtent, "", "NEW_SELECTION")
        numHouses = gp.GetCount_management("pLayer")
        return numHouses

## Identify core patches based on minimum size criteria
def identifyCores(theTable, workspace, minCoreSize):
        arcpy.AddMessage("\tIdentifying habitat cores...")

        if workspace == '':
                workspace = arcpy.env.Workspace

        elif workspace == "in_memory":
                sumTable = workspace + os.sep + "Core_summary"
                iTable = workspace + os.sep + "Iteration_summary"
        else:
                sumTable = workspace + os.sep + "Core_summary.dbf"
                iTable = workspace + os.sep + "Iteration_summary.dbf"

##        cLayer = workspace + os.sep + "coresLayer"
        cLayer = "coresLayer"
##        sumTable = workspace+os.sep + "Core_summary.dbf"
##        iTable = workspace+os.sep + "Iteration_summary.dbf"
        whereClause = '"AREA" < ' + str(minCoreSize)
             
        
        
        # Process: Add 'Area' and 'Prcnt Area'Fields ...
        if not arcpy.ListFields(theTable,"CORE"):
##        if not arcpy.ListFields(theTable,"CORE").Next():
##                arcpy.AddField_management(theTable, "CORE", "SHORT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
                arcpy.AddField_management(theTable, "CORE", "SHORT")
                
        # Process: Calculate Area...
        CalcArea(theTable)
        
        # Process: Make Feature Layer ...
        arcpy.MakeFeatureLayer_management(theTable, cLayer, "", "", "Input_FID Input_FID VISIBLE NONE")
        
        # Process: Select Area GT Core...
        arcpy.SelectLayerByAttribute_management(cLayer, "NEW_SELECTION", whereClause)

        # Process: Identify Non-Cores...
        arcpy.CalculateField_management(cLayer, "CORE", "0", "PYTHON", "")

        # Process: Switch Selection...
        arcpy.SelectLayerByAttribute_management(cLayer, "SWITCH_SELECTION", "")

##        global numCores
        numCores = arcpy.GetCount_management(cLayer)

        if numCores > 0:

                # Process: Identify Cores...
                arcpy.CalculateField_management(cLayer, "CORE", "1", "VB", "")
                
                # Process: Summary Statistics on selected set...
                arcpy.Statistics_analysis(cLayer, sumTable, "AREA SUM;AREA MIN;AREA MAX;AREA MEAN;AREA RANGE;AREA STD", "CORE")

                # Process: Summary Statistics (5)...
                arcpy.Statistics_analysis(sumTable, iTable, "SUM_AREA SUM;SUM_AREA MAX;SUM_AREA MIN", "")

                # Process: Clear Selection...
                arcpy.SelectLayerByAttribute_management(cLayer, "CLEAR_SELECTION", "")

        else:
                # Make feature class for later use...
                arcpy.AddWarning("\tNo cores contained in simulation...")

                 # Process: Clear Selection...
                arcpy.SelectLayerByAttribute_management(cLayer, "CLEAR_SELECTION", "")
                
                # Process: Summary Statistics on selected set...
                arcpy.Statistics_analysis(cLayer, sumTable, "AREA SUM;AREA MIN;AREA MAX;AREA MEAN;AREA RANGE;AREA STD", "CORE")

                # Process: Summary Statistics (5)...
                arcpy.Statistics_analysis(sumTable, iTable, "SUM_AREA SUM;SUM_AREA MAX;SUM_AREA MIN", "")

        # Make feature class for later use...
        arcpy.AddMessage("\tCreating temporary feature class from layer...")
        arcpy.CopyFeatures_management(cLayer, "in_memory\\xxcoresLayer")
        arcpy.Delete_management(cLayer,"")

        return numCores,iTable #send iteration table back for further processing

## Make output table and add appropriate fields
def MakeTable(outWorkspace, outTable):
        nTable = outWorkspace + os.sep + outTable
        arcpy.CreateTable_management(outWorkspace, outTable)
        # Process: Add Fields...
        arcpy.AddField_management(nTable, "SUM_SUM_AR", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(nTable, "MAX_SUM_AR", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(nTable, "MIN_SUM_AR", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(nTable, "HOUSES", "LONG", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(nTable, "PRCNT_AREA", "DOUBLE", "18", "4", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        return nTable

## Create random pattern within a constraint layer
def RandomPattern(Workspace, aExtent, constraint, numHouses):
        # Process: Create Random Points...
        gp.AddMessage("\tCreating random pattern of " + str(numHouses) + " houses...")
        rPointsName = "random_temp.shp" #Random points for analysis        
        rPoints = Workspace + os.sep + rPointsName

        # This nonsense works around a persistent schema lock that appears on the 3rd iteration
        # and prevents rPonts from being overwritten or deleted.
        if gp.Exists(rPoints):
                try:
                        gp.AddMessage("\tAttempting to delete previous replication output...")
                        gp.Delete_management(rPoints,"")
                        gp.AddMessage(rPoints + " succesfully deleted...")
                except:
                        n = 1
                        gp.AddWarning("\tUnable to delete previous replication output...\n\tCreating unique filename for output...")
                        desc = gp.Describe(rPoints)
                        bName = desc.BaseName
                        path = desc.Path
                        while gp.Exists(rPoints) == 1:
                                rPointsName = bName + str(n) + ".shp"
                                rPoints = path + os.sep + rPointsName
                                n = n + 1
                        
        rPoints = gp.CreateRandomPoints_management(Workspace, rPointsName, constraint, aExtent, numHouses, "0 Meters", "POINT", "0")
        
        return rPoints

#Rescale input rasters to a range of 0-100
def Rescale(inRaster,outRaster,tWorkspace):
    
    tmpRaster = tWorkspace + os.sep + 'tmpRast1'
    tmpRaster2 = tWorkspace + os.sep + 'tmpRast2'
    
    gp.AddMessage("Rescaling " + inRaster + "...")        
    Type = gp.GetRasterProperties (inRaster, "VALUETYPE")
    if Type == '1':
        gp.AddMessage("\tConverting integer raster to floating point...")
        gp.Float_sa(inRaster, tmpRaster)
        inRaster = tmpRaster
    Rmin = gp.GetRasterProperties (inRaster, "MINIMUM")
    
    #if the input has negative values. \
    if float(Rmin < 0):
        gp.AddWarning("\tValues less than zero found. Adjusting by adding " + str(Rmin) + "to input...")
        Rmin = str(abs(Rmin))
        gp.Plus_sa(inRaster, Rmin, tmpRaster)
        inRaster = tmpRaster
    elif float(Rmin > 0):
        Rmin = str(Rmin)
        gp.AddMessage("\tSubtracting " + Rmin + " from input raster...")
        gp.Minus_sa(inRaster, Rmin, tmpRaster)
        inRaster = tmpRaster
        
    Rmax = str(gp.GetRasterProperties (inRaster, "MAXIMUM"))
    gp.AddMessage("\tDividing input raster by it's maximum value...")
    gp.AddMessage("\tMaximum Value = " + Rmax)
    gp.Divide_sa(inRaster, Rmax, tmpRaster2)
    gp.AddMessage("\tMultiplying input raster by 100...")
    gp.Times_sa(tmpRaster2, '100', outRaster)
    gp.AddMessage("\tRescale complete...")
    return(outRaster)


# Run Monte Carlo Simulation for corridor density
def RunCorridorSimulation (i,tWorkspace, aExtent, constraint, numHouses, nTable, theField, pointLayer):

        #create local variables

        w = tWorkspace.rsplit(os.sep + "",1)
        outWorkspace = w[0]
        theField = "MEAN_" + theField
        theField = gp.ValidateFieldName(theField)

        if gp.Exists(tWorkspace + os.sep + "random_temp.shp"):
                gp.Delete_management(tWorkspace + os.sep + "random_temp.shp","")


        #Iterate model for monte carlo simulation of random points
        for a in range(i):
                gp.AddMessage("\tRunning " + str(a + 1) + " of " + str(i) + " iterations...")

                CleanFiles(tWorkspace)

                rPoints = RandomPattern(tWorkspace, aExtent, constraint, numHouses)

                # Identify core areas based on user-defined area requirements
                iTable = CalcCorridorDistance(rPoints, tWorkspace)

                if iTable:
                        # Add Houses Field
                        gp.AddField_management(iTable, "HOUSES", "LONG", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")

                        # Process: Calculate Number of Houses...
                        gp.CalculateField_management(iTable, "HOUSES", str(numHouses), "VB", "")
                        
                        # Process: Append...
                        gp.Append_management(os.path.join(tWorkspace, iTable), nTable, "NO_TEST", "", "")

                else:
                        return False
                
        # Determine mean area conserved under simulation run
        gp.AddMessage("\tCalculating mean corridor width...")
        meanDist = CalcMeans(nTable, tWorkspace, theField, numHouses)
        return meanDist

# Run Monte Carlo Simulation
# NOTE: random variable is no longer used but may be converted to a 'pattern' variable to allow users to select other point patterns
def RunSimulation (i,tWorkspace, aExtent, constraint, numHouses, minCoreSize, stinfDistance, rdinfDistance, nTable, tArea, pointLayer, roadLayer, random, outRdCost, backlink, existRoads):

        #create local variables

##        w = tWorkspace.rsplit(os.sep + "",1)
##        outWorkspace = w[0]

        #Iterate model for monte carlo simulation of random points
        for a in range(i):
                arcpy.AddMessage("\tRunning " + str(a + 1) + " of " + str(i) + " replications...")

                if pointLayer == "" or pointLayer == "#":
                        pointDelete = 1                                        
                        rPoints = RandomPattern(tWorkspace, aExtent, constraint, numHouses)
                        desc = arcpy.Describe(rPoints)
                        t = desc.DataType
                        c = desc.CatalogPath
                else:
                        pointDelete = 0
                        rPoints = pointLayer

                # Generate simulated road network
                if roadLayer == '':
                        lineDelete = 1
                        lineLayer = tWorkspace + os.sep + "lineLayer.shp"
                        simRoads = GenRoads (rPoints, outRdCost, backlink, "")
                        arcpy.AddMessage("Merging simulated roads with existing roads...")
##                        arcpy.AddWarning("CreateObject")
##                        vt = arcpy.createobject("ValueTable")
##                        arcpy.AddWarning("AddRow")
##                        vt.AddRow(existRoads)
##                        arcpy.AddWarning(vt)
##                        vt.AddRow(simRoads)
##                        arcpy.AddWarning(vt)
##                        roadLayer = arcpy.Merge_management (vt, lineLayer, "")
                        roadLayer = arcpy.Merge_management (existRoads + ";" + simRoads, lineLayer, "")
                else:
                        lineDelete = 0

                # Process point pattern to simulate habitat cores
                sCores = GenerateCores(rPoints, tWorkspace, aExtent, stinfDistance, rdinfDistance, roadLayer) #, inDEM, existRoads)
                # If sCores = False, then extent is empty and percent area conserved equals zero.
                if sCores:

                        # Identify core areas based on user-defined area requirements
                        iTable = identifyCores(sCores, tWorkspace, minCoreSize)
                        numCores = iTable[0]
                        iTable = iTable[1]
                        
                        # Add Houses Field
                        arcpy.AddField_management(iTable, "HOUSES", "LONG", "", "", "", "", "", "NON_REQUIRED", "")

                        # Process: Add Prcnt Area Field...
                        arcpy.AddField_management(iTable, "PRCNT_AREA", "DOUBLE", "18", "4", "", "", "", "NON_REQUIRED", "")

                else:
                        iTable = MakeTable(tWorkspace, "Iteration_summary")
                        numCores = 0
                        
                # Process: Calculate Number of Houses...
                arcpy.CalculateField_management(iTable, "HOUSES", str(numHouses), "PYTHON", "")

                # Process: Calculate Percent Area...
                if numCores > 0:
                        arcpy.CalculateField_management(iTable, "PRCNT_AREA", "(!MAX_SUM_AREA! /" + str(tArea) + ")*100", "PYTHON", "")
                else:
                        arcpy.CalculateField_management(iTable, "PRCNT_AREA", "0", "PYTHON", "")

                # Process: Append...
                arcpy.Append_management(iTable, nTable, "NO_TEST", "", "")

                #Reset temporary variables
                if pointDelete == 1:
                        del pointLayer
                        pointLayer = ''

                if lineDelete == 1:
                        del roadLayer
                        roadLayer = ''
                
        # Determine mean area conserved under simulation run
        arcpy.AddMessage("\tCalculating percent area conserved...")
        meanPrcntArea = CalcMeans(nTable, tWorkspace, "PRCNT_AREA", numHouses)
        arcpy.AddMessage(str(int(meanPrcntArea)) + "% of area conserved with " + str(numHouses) + " houses...\n")
        
        return meanPrcntArea, rPoints
