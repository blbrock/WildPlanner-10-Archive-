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
# gen_roads.py
# Created on: Thu Mar 03 2011 03:23:33 PM
#   (generated by ArcGIS/ModelBuilder)
# Usage: gen_roads <inDEM> <inSlope> <backlink> <inRoads> <inPoints> <roadCost> <Output_polyline_features> 
#
#
#
# Written by Brent L. Brock, Landscape Ecologist
# Craighead Institute
# bbrock@craigheadresearch.org
#
# Creates a hypothetical road connection between planned or hypothetical structure locations and existing roads.
# Roads are placed along the least-cost path that minimizes gradient and distance.
#
#This module depends on functionlib.py which contains functions required for processing this script
# 
# ---------------------------------------------------------------------------


# Import system modules
import sys, string, os, arcgisscripting
from functionlib import CreateTempWorkspace, CleanFiles, CreateRoadCost, GenRoads

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

# Script arguments...
inDEM = sys.argv[1]
inSlope = sys.argv[2]
inRoads = sys.argv[3]
inPoints = sys.argv[4]
outRoads = sys.argv[5]
appendRoads = sys.argv[6]
aExtent = sys.argv[7]
# Set output workspace
r = outRoads.rsplit(os.sep,1)
outWorkspace = r[0]
gp.Workspace = outWorkspace
gp.Extent = aExtent

rdcst = CreateRoadCost(inDEM, inSlope, inRoads, "")
outRdCost = rdcst[0]
backlink = rdcst[1]

if appendRoads == "true":
    gp.AddMessage("Merging new roads with existing roads...")
    if gp.ScratchWorkspace:
        newRoads = gp.ScratchWorkspace + os.sep + "xxnewroads.shp"
    else:
        newRoads = gp.Workspace + os.sep + "xxnewroads.shp"
    gp.MakeFeatureLayer_management(inRoads,"existRoads")
    gp.AddField_management("existRoads", "SIM_RD", "SHORT")
    gp.CalculateField_management ("existRoads", "SIM_RD", "0")
    GenRoads (inPoints, outRdCost, backlink, newRoads)
    gp.MakeFeatureLayer_management(newRoads,"newRoads")
    gp.AddField_management("newRoads", "SIM_RD", "SHORT")
    gp.CalculateField_management ("newRoads", "SIM_RD", "1")
    gp.Merge_management ("existRoads;newRoads", outRoads)
    gp.Delete_management(newRoads)
#    gp.Merge_management (inRoads + ";" + newRoads, outRoads)
else:
    GenRoads (inPoints, outRdCost, backlink, outRoads)

try:
    gp.SetParameterAsText(4, outRoads)
    params = gp.GetParameterInfo()
except:
    gp.getMessages(2)
