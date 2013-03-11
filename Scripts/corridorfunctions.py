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
# functionlib.py (version 1.0 beta)
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
import os, random, arcgisscripting, sys
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
        gp.AddMessage("Calculating area of polygon cores...")
        if not gp.ListFields(inputFC,"AREA").Next():
                gp.AddField_management(inputFC, "Area", "DOUBLE", "18", "4", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        # Calculate 'Area' field to shape area
        expression = "float(!shape.area!)"
        gp.CalculateField_management(inputFC, "Area", expression, "PYTHON")





## Calculate total area from a table containing an 'Area' field
def CalcTotalArea(workspace, extentLayer):
        gp.AddMessage("Calculating total area of " + extentLayer +"...")
        if workspace == "":
                workspace = gp.Workspace
        aTable = workspace + "\\aTable.dbf"
        
        # Process: Summary Statistics...
        gp.Statistics_analysis(extentLayer, aTable, "AREA SUM", "")
        # Get value from summary table
        rows = gp.SearchCursor(aTable)
        row = rows.Next()
        tArea = row.GetValue("SUM_AREA")
        # Delete the row and cursor
        del row, rows
        try:
                gp.Delete_management(aTable,"")
        except:
                gpAddWarning("Unable to delete " + aTable + "...")
        
        return tArea





# Calculates the mean of all records for the input table ('inTable') and specified field ('theField')
def CalcMeans(inTable, workspace, theField, numHouses):
        gp.AddMessage("Calculating means...")
        if workspace == '':
                workspace = gp.Workspace
        mTable = workspace + "\\mean_value.dbf"
        
##        try:
        gp.Statistics_analysis(inTable, mTable, theField + " MEAN", "HOUSES")
        # Get value from summary table
        theField = "MEAN_" + theField
        theField = theField[0:10]
        rows = gp.SearchCursor(mTable, "HOUSES = " + str(numHouses))
        row = rows.Next()
        theMean = row.GetValue(theField) 
##        
        
##        Delete the row and cursor
        del row, rows
##        except:
##                gp.AddWarning("Unable to calculate mean  area.")
##                gp.AddWarning("   This is mostly likely due to zero cores created in simulation.")
##                gp.AddWarning("   Assuming mean = 0 and continuing simulations...")
##                theMean = 0
        return theMean





## Clean up tables and feature classes in temporary workspace
def CleanFiles(workspace):
        gp.AddMessage("Deleting temporary files...")
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






#### Convert to acres.  Converts square meters to acreage
##def ConvertAcres(meters):
##        gp.AddMessage("Converting area units to acres...")
##        acres = meters * 0.000247105381467165 
##        return float(acres)




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
        gp.AddError("Unknown units " + units + " could not convert to square meters")
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
        gp.AddError("Unknown units " + units + " could not convert from square meters")
    return areaValue





## Convert distance units to meters
def ConvertDistanceToMeters(distanceValue, units):
    distanceValue = float(distanceValue)
    units = units.replace(" ", "")
    
    if units == 'Centimeters':
        distanceValue = distanceValue * 0.01
    elif units == 'DecimalDegrees':
        gp.AddError("Cannot convert from Decimal Degrees.  Choose different units.")
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
        gp.AddError("Unknown units " + units + " could not convert to meters")
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
        gp.AddError("Unknown units " + units + " could not convert to meters")
    return distanceValue



## Create a temporary workspace
def CreateTempWorkspace(outWorkspace):
        # Set up temporary workspace
        if outWorkspace == "":
                outWorkspace = gp.Workspace
        if outWorkspace == "" and gp.Workspace == None:
                gp.AddError("Output Workspace is not defined...")
        gp.Workspace = outWorkspace
        if gp.ScratchWorkspace == None:
            gp.ScratchWorkspace = outWorkspace
        WorkspaceName = "\\xx" + str(random.randrange(0,1000000))
        Workspace = gp.ScratchWorkspace + WorkspaceName
        gp.AddMessage("Creating temporary workspace " + Workspace + "...")
        if not gp.Exists(Workspace):
                gp.CreateArcInfoWorkspace(gp.ScratchWorkspace, WorkspaceName)
        else:
                gp.AddWarning(Workspace + " already exists...")
        return Workspace



## Simulate undisturbed habitat areas created by a given input point pattern.
## Results are clipped to boundaries of polygons in aExtent.
def GenerateCores (pointLayer, workspace, aExtent, buffDistance, roads):
        gp.AddMessage("Generating habitat cores...")

        if workspace == '':
                workspace = gp.Workspace

        pBuff = workspace + "\\pointBuffer.shp"
        lBuff = workspace + "\\lineBuffer.shp"
        polyLayer = workspace + "\\polyLayer.shp"
        lineLayer = workspace + "\\lineLayer.shp"
        eraseLayer1 = workspace + "\\eraseLayer1.shp"
        eraseLayer2 = workspace + "\\eraseLayer2.shp"
        clip = workspace + "\\clipLayer.shp"
        habitat = workspace + "\\habitatPatches.shp"
        
        # Set the analysis extent to aExtent
        desc = gp.describe(aExtent)
        extent = desc.Extent
        gp.Extent = extent

        
        
        
        
        if roads == '':
                # Process: Buffer points by disturbance distance...
                gp.Buffer_analysis(pointLayer, pBuff, buffDistance, "FULL", "ROUND", "NONE", "")
                # Process: Create Thiessen Polygons from random points...
                gp.CreateThiessenPolygons_analysis(pointLayer, polyLayer, "ONLY_FID")
                # Process: Erase point (houses) disturbance buffers from thiessen polygons ...
                gp.AddMessage("Removing point buffers from analysis extent...")
                gp.Erase_analysis(polyLayer, pBuff, eraseLayer1, "")
                # Process: Convert polygons to lines as surrogate for roads ...
                gp.PolygonToLine_management(polyLayer, lineLayer)
        else:
                polyLayer = aExtent
                lineLayer = roads
                # Process: Buffer points by disturbance distance...
                gp.AddMessage("Buffering structures by disturbance distance...")
                gp.MakeFeatureLayer_management(pointLayer, "pLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
                gp.SelectLayerByLocation("pLayer", "COMPLETELY_WITHIN", aExtent, "", "NEW_SELECTION")
                gp.Buffer_analysis("pLayer", pBuff, buffDistance, "FULL", "ROUND", "NONE", "")
                # Process: Erase point (houses) disturbance buffers from thiessen polygons ...
                gp.AddMessage("Removing point buffers from analysis extent...")
                gp.Erase_analysis(polyLayer, pBuff, eraseLayer1, "")             

        # Process: Buffer lines (roads) by disturbance distance ...
        gp.AddMessage("Buffering roads by disturbance distance...")
        gp.Buffer_analysis(lineLayer, lBuff, buffDistance, "FULL", "ROUND", "NONE", "")

        # Process: Erase line buffers from thiessen polygons (note: polygons have now been erased for both houses and roads) ...
        gp.AddMessage("Removing road buffers from analysis extent...")
        gp.Erase_analysis(eraseLayer1, lBuff, eraseLayer2, "")

        # Process: Clip cores by analysis extent ...
        gp.Clip_analysis(eraseLayer2, aExtent, clip, "")

        # Make sure features are single part for proper area calculations
        gp.AddMessage("Generating single part features...")
        gp.MultipartToSinglepart_management (clip, habitat)
        
        return habitat




## Get the number of structure points from a point layer
def GetNumHouses(pointLayer, aExtent):
        gp.AddMessage("Calculating the number of structures within analysis extent...")
        gp.MakeFeatureLayer_management(pointLayer, "pLayer", "", "", "Input_FID Input_FID VISIBLE NONE")
        gp.SelectLayerByLocation("pLayer", "COMPLETELY_WITHIN", aExtent, "", "NEW_SELECTION")
        numHouses = gp.GetCount_management("pLayer")
        return numHouses






## Identify core patches based on minimum size criteria
def identifyCores(theTable, workspace, minCoreSize):
        gp.AddMessage("Identifying habitat cores...")

        if workspace == '':
                workspace = gp.Workspace

        cLayer = workspace + "\\coresLayer"
        sumTable = workspace+"\\Core_summary.dbf"
        iTable = workspace+"\\Iteration_summary.dbf"
        whereClause = '"AREA" < ' + str(minCoreSize)
             
        
        
        # Process: Add 'Area' and 'Prcnt Area'Fields ...   
        if not gp.ListFields(theTable,"CORE").Next():
                gp.AddField_management(theTable, "CORE", "SHORT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
                
        # Process: Calculate Area...
        CalcArea(theTable)
        
        # Process: Make Feature Layer ...
        gp.MakeFeatureLayer_management(theTable, cLayer, "", "", "Input_FID Input_FID VISIBLE NONE")
        
        # Process: Select Area GT Core...
        gp.SelectLayerByAttribute_management(cLayer, "NEW_SELECTION", whereClause)

        # Process: Identify Non-Cores...
        gp.CalculateField_management(cLayer, "CORE", "0", "VB", "")

        # Process: Switch Selection...
        gp.SelectLayerByAttribute_management(cLayer, "SWITCH_SELECTION", "")

        numCores = gp.GetCount_management(cLayer)

        if numCores > 0:

                # Process: Identify Cores...
                gp.CalculateField_management(cLayer, "CORE", "1", "VB", "")
                
                # Process: Summary Statistics on selected set...
                gp.Statistics_analysis(cLayer, sumTable, "AREA SUM;AREA MIN;AREA MAX;AREA MEAN;AREA RANGE;AREA STD", "CORE")

                # Process: Summary Statistics (5)...
                gp.Statistics_analysis(sumTable, iTable, "SUM_AREA SUM;SUM_AREA MAX;SUM_AREA MIN", "")

                # Process: Clear Selection...
                gp.SelectLayerByAttribute_management(cLayer, "CLEAR_SELECTION", "")
                
                # Make feature class for later use...
                gp.AddMessage("Creating temporary feature class from layer...")
                gp.CopyFeatures(cLayer, workspace + "\\xxcoresLayer.shp")

                gp.Delete_management(cLayer,"")

                return iTable #send iteration table back for further processing
        else:
                # Make feature class for later use...
                gp.AddWarning("Maximum buildout has been reached...")
                gp.AddWarning("Ending simulation and summarinzing results...")
                gp.AddMessage("Creating temporary feature class from layer...")
                gp.CopyFeatures(cLayer, workspace + "\\xxcoresLayer.shp")
                return False





# Run Monte Carlo Simulation
def RunSimulation (i,tWorkspace, aExtent, constraint, numHouses, minCoreSize, distDistance, nTable, tArea, pointLayer, roadLayer, random):

        #create local variables

        w = tWorkspace.rsplit("\\",1)
        outWorkspace = w[0]


        #Iterate model for monte carlo simulation of random points
        for a in range(i):
                gp.AddMessage("Running " + str(a + 1) + " of " + str(i) + " iterations...")

                CleanFiles(tWorkspace)

                if random == True:
                        # Process: Create Random Points...
                        gp.AddMessage("Creating random pattern of " + str(numHouses) + " houses...")
                        rPointsName = "random_temp.shp" #Random points for analysis
                        rPoints = tWorkspace + "\\" + rPointsName
                        gp.CreateRandomPoints_management(tWorkspace, rPointsName, constraint, aExtent, numHouses, "0 Meters", "POINT", "0")
                else:
                        rPoints = pointLayer

                # Process point pattern to simulate habitat cores
                sCores = GenerateCores(rPoints, tWorkspace, aExtent, distDistance, roadLayer)

                # Identify core areas based on user-defined area requirements
                iTable = identifyCores(sCores, tWorkspace, minCoreSize)

                if iTable:

                        # Add Houses Field
                        gp.AddField_management(iTable, "HOUSES", "LONG", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")

                        # Process: Calculate Number of Houses...
                        gp.CalculateField_management(iTable, "HOUSES", str(numHouses), "VB", "")

                        # Process: Add Prcnt Area Field...
                        gp.AddField_management(iTable, "PRCNT_AREA", "DOUBLE", "18", "4", "", "", "NON_NULLABLE", "NON_REQUIRED", "")

                        # Process: Calculate Percent Area...
                        gp.CalculateField_management(iTable, "PRCNT_AREA", "([MAX_SUM_AR] /" + str(tArea) + ")*100", "VB", "")

                        # Process: Append...
                        gp.Append_management(iTable, nTable, "NO_TEST", "", "")

                else:
                        return False

                
        # Determine mean area conserved under simulation run
        gp.AddMessage("Calculating percent area conserved...")
        meanPrcntArea = CalcMeans(nTable, tWorkspace, "PRCNT_AREA", numHouses)
        gp.AddMessage(str(int(meanPrcntArea)) + "% of area conserved with " + str(numHouses) + " houses...")
        return meanPrcntArea

