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
# build_cost_surface.py (version 1.0 beta)
# Created on: Mon Nov 24 2008 04:08:09 PM
#
# Written by Brent L. Brock, Landscape Ecologist
# Craighead Institute
# bbrock@craigheadresearch.org
#
# Rescales input raster so output raster values are 0-100.  If the invert values option is selected, the
# output raster values will be inverted so the maximum value of the input will be zero in output and the input
# minimum value will be 100 in the output.
#
#This module depends on functionlib.py which contains functions required for processing this script
# 
# ---------------------------------------------------------------------------

# Import system modules
import sys, os, arcpy
# from functionlib import CreateTempWorkspace, CleanFiles
# from link_functions import HabToCost

# Create the Geoprocessor object
#gp = arcgisscripting.create(9.3)
##arcpy.Workspace = None

# Set the necessary product code
arcpy.SetProduct("ArcInfo")

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

arcpy.AddMessage("Setting analysis environment...")

# Load required toolboxes...
arcpy.overwriteoutput = 1

# Script arguments...

cost1 = sys.argv[1]
rescale = sys.argv[2]
Min = sys.argv[3]
Max = sys.argv[4]
invertCost1 = sys.argv[5]
outRaster = sys.argv[6]

# Convert max and min to multplier value
Max = str(float(Max)-float(Min))
# Set up workspaces
r = outRaster.rsplit(os.sep,1)
outWorkspace = r[0]

###-------------------------------CLASSES------------------------------------###

class RoundError(Exception): pass

# --------------------- Functions -------------------------#
# Invert habitat quality layer to create cost surface
def HabToCost (habLayer, costSurf):
        # Convert Habitat model to cost surface
    arcpy.AddMessage("Inverting values in " + habLayer + "...")
    Habmax = arcpy.GetRasterProperties_management (habLayer, "MAXIMUM").getOutput(0)
##    invHab =  "in_memory\\inv_hab"
##    plusHab = "in_memory\\plus_hab"
    invHab = arcpy.sa.Times(habLayer, -1)
    plusHab = arcpy.sa.Plus(invHab, float(Habmax))
    
    #ArcGIS may introduce negative values through rounding errors so this is
    #a crude way to correct the problem.
    
    Costmin = arcpy.GetRasterProperties_management (plusHab, "MINIMUM").getOutput(0)
    if float(Costmin) < 0:
        if float(Costmin) < 0 and float(Costmin) > -.1:
            arcpy.AddWarning("The cost surface contains negative values. Rounding error assumed. Minimum Value will be truncated to zero...")
            arcpy.sa.Con(plusHab, "0", costSurf, plusHab, "Value < 0")
        else:
            string = 'Input contains negative values that cannot be attributed to rounding errors. Check inputs'
            roundError = 'Error'
            raise RoundError(string)
        
    else:
        arcpy.CopyRaster_management(plusHab,costSurf)

    return costSurf


#Rescale input rasters to a range of 1-100
def Rescale(inRaster,outRaster, Min, Max):
    
##    tmpRaster = tWorkspace + os.sep + 'tmpRast1'
##    tmpRaster2 = tWorkspace + os.sep + 'tmpRast2'
##    tmpRaster3 = tWorkspace + os.sep + 'tmpRast3'
##    tmpRaster4 = tWorkspace + os.sep + 'tmpRast4'
    
    arcpy.AddMessage("Rescaling " + inRaster + "...")        
    Type = arcpy.GetRasterProperties_management (inRaster, "VALUETYPE").getOutput(0)
    if Type == '1':
        arcpy.AddMessage("\tConverting integer raster to floating point...")
        inRaster = arcpy.sa.Float(inRaster)
        # inRaster = tmpRaster
    Rmin = arcpy.GetRasterProperties_management (inRaster, "MINIMUM").getOutput(0)
    
    #if the input has negative values. \
    if float(Rmin < 0):
        arcpy.AddWarning("\tValues less than zero found. Adjusting by adding " + Rmin + "to input...")
        Rmin = str(abs(Rmin))
        inRaster = arcpy.sa.Plus(inRaster, float(Rmin))
        #arcpy.Delete_management("in_memory\\tmpRaster")
        #inRaster = tmpRaster
        #arcpy.Delete_management(tmpRaster)
    elif float(Rmin > 0):
        arcpy.AddMessage("\tSubtracting " + Rmin + " from input raster...")
        inRaster = arcpy.sa.Minus(inRaster, float(Rmin))
    #arcpy.Delete_management("in_memory\\tmpRaster")
        #inRaster = tmpRaster4
        
    Rmax = arcpy.GetRasterProperties_management (inRaster, "MAXIMUM").getOutput(0)
    arcpy.AddMessage("\tDividing input raster by it's maximum value...")
    arcpy.AddMessage("\tMaximum Value = " + Rmax)
    inRaster = arcpy.sa.Divide(inRaster, float(Rmax))
    arcpy.AddMessage("\tMultiplying input raster by " + Max + "...")
    #arcpy.Delete_management("in_memory\\tmpRaster2")
    inRaster = arcpy.sa.Times(inRaster, float(Max))
    arcpy.AddMessage("\tAdding " + Min + " to input raster...")
    plusRaster = arcpy.sa.Plus(inRaster, float(Min))
    plusRaster.save(outRaster)
    arcpy.AddMessage("\tRescale complete...")
    return(outRaster)


###------------------------------------ END FUNCTIONS -------------------------


# Set up temporary workspace
##scratchWS = arcpy.scratchWorkspace
##if scratchWS:
##    desc = arcpy.Describe(scratchWS)
##    if desc.workspaceType == u'FileSystem':
##        tWorkspace = CreateTempWorkspace(scratchWS)
##    else:
##        tWorkspace = CreateTempWorkspace(outWorkspace)
##else:
##    tWorkspace = CreateTempWorkspace(outWorkspace)
##arcpy.Workspace = outWorkspace
##tmpRaster = tWorkspace + os.sep + 'tmpRast'


# Rescale inputs if requested.  Rescaled inputs will be stored in outRescale.
if rescale == 'true':
    cost2 = Rescale(cost1,outRaster,Min,Max)
    cost1 = cost2


# Invert any inputs with this option selected
if invertCost1 == 'true':
    cost2 = HabToCost(cost1,outRaster)
    
arcpy.AddMessage("Outraster is: " + outRaster)
# Add Cost Raster to Display

cost2 = outRaster

try:
    arcpy.SetParameterAsText(5, outRaster)
    params = arcpy.GetParameterInfo()
except:
    arcpy.getMessages(2)

# Clean up temporary workspace
# CleanFiles(tWorkspace)
arcpy.Delete_management("in_memory")
del arcpy

