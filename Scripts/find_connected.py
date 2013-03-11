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
import sys, os, arcpy
from itertools import combinations, chain

link = sys.argv[1]
fullBuild = sys.argv[2]
source = sys.argv[3]
pField = sys.argv[4]
curSurface = sys.argv[5]
baseName = sys.argv[6]

outConnect = baseName + '_connect.shp'
outLCP = baseName + '_lcp.shp'
outPAreas = baseName + '_pAreas.shp'

###-------------------------------CLASSES------------------------------------###

class RoundError(Exception): pass

# --------------------- Functions -------------------------#
# Invert habitat quality layer to create cost surface
def invertRaster (habLayer, costSurf):
        # Convert Habitat model to cost surface
    arcpy.AddMessage("Inverting values in " + habLayer + "...")
    Habmax = arcpy.GetRasterProperties_management (habLayer, "MAXIMUM").getOutput(0)
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

# Creates a list of pair-wise combinations of unique values in input field
def getCombinations(inTable, pField):
    # Get a list of unique values in pFields
    valueList = []
    rows = arcpy.SearchCursor(inTable)
    for row in rows:
        valueList.append(row.getValue(pField))

    uniqueSet = set(valueList)
    uniqueList = list(uniqueSet)
    uniqueList.sort()

    del rows
    del row

    # Process each pair-wise combination of patches to find connected areas

    combs = list(combinations(uniqueList, 2))
    return uniqueList, combs

# Function findConnected generatess a polygon layer of uncompromised regions
# connecting patch pairs.
def findConnected(link,source,pField,outConnect):
    global connectList, noConnectList
    connectList = list()

    # Create a list of patch pairs that are not connected.
    # noConnectList is currently created but not used.
    noConnectList = list()

    # If source layer is raster, convert to polygon
    if 'Raster' in arcpy.Describe(source).dataType:
        arcpy.AddMessage('Converting source patch raster to polygons...')
        source = arcpy.RasterToPolygon_conversion(source, "in_memory\\source", "NO_SIMPLIFY", pField )
        pField = 'grid_code'

    # Convert linkage layer to shapefile
    arcpy.AddMessage('Converting linkage layer to polygons...')
    link = attExtract = arcpy.sa.ExtractByAttributes(link, "VALUE = 1") 
    linkPoly = arcpy.RasterToPolygon_conversion(link, "in_memory\\linkPoly", "NO_SIMPLIFY", "VALUE")
    linkPoly = arcpy.Erase_analysis(linkPoly, source, "in_memory\\linkErase", "")
    linkPoly = arcpy.MultipartToSinglepart_management(linkPoly, "in_memory\\linkSingle")
    arcpy.MakeFeatureLayer_management(linkPoly, "linkLayer")
    arcpy.MakeFeatureLayer_management(source, "sourceLayer")

    arcpy.AddMessage('Processing patch pairs to find connected linkages...')

    # Get a list of pair-wise combinations of valued in pField
    combs = getCombinations(source, pField)[1]
    
    arcpy.AddMessage('\tThe following patch pairs will be processed: ' + str(combs) + '...')

    for c in combs:
        name = 'connect_' + str(int(c[0])) + '_' + str(int(c[1]))    
            
        # Select linkage areas touching first patch
        arcpy.AddMessage('\t\t\tSelecting linkages connected to first patch...')
        arcpy.SelectLayerByAttribute_management("sourceLayer", "NEW_SELECTION", pField + " = " + str(c[0]) )
        arcpy.SelectLayerByLocation_management("linkLayer", "BOUNDARY_TOUCHES", "sourceLayer", "", "NEW_SELECTION" )

        #Select linkage areas touching first AND second patch
        arcpy.AddMessage('\t\t\tSelecting linkages connected to second patch...')
        arcpy.SelectLayerByAttribute_management("sourceLayer", "NEW_SELECTION", pField + " = " + str(c[1]) )
        arcpy.SelectLayerByLocation_management("linkLayer", "BOUNDARY_TOUCHES", "sourceLayer", "", "SUBSET_SELECTION" )

        #Copy selected features to new feature class
        
        selected = int(arcpy.GetCount_management("linkLayer").getOutput(0))
        if int(arcpy.GetCount_management("linkLayer").getOutput(0)) > 0:
            arcpy.AddMessage('\t\tConnected linkage between patch ' +  str(int(c[0])) + ' and ' + str(int(c[1])) + ' found...')
            arcpy.AddMessage("\t\t\tProcessing linkage '" + name + "'...")
            arcpy.MakeFeatureLayer_management("linkLayer", name )
            arcpy.AddField_management("linkLayer", "CONNECTION", "TEXT", "", "", "50")
            arcpy.CalculateField_management("linkLayer", "CONNECTION", '"' + name + '"' )
            connectList.append(name)
        else:
            arcpy.AddWarning('\t\tNo linkage between patch ' +  str(int(c[0])) + ' and ' + str(int(c[1])) + ' skipping patch pair...')
            noConnectList.append(name)
            
        arcpy.SelectLayerByAttribute_management("linkLayer", "CLEAR_SELECTION", "" )
        arcpy.SelectLayerByAttribute_management("sourceLayer", "CLEAR_SELECTION", "" )
    if len(connectList) > 0:
        if len(connectList) > 1:
            arcpy.AddMessage('Merging connected linkages : ' + str(connectList) + '...')
            cMerge = arcpy.Merge_management(connectList, "in_memory\\cMerge")
            arcpy.AddMessage('\tDissolving features...')
            outConnect = arcpy.Dissolve_management("in_memory\\cMerge", outConnect, '', "", "SINGLE_PART", "")

        else:
           arcpy.AddMessage('\tCopying connected linkage features to shapefile...')
           outConnect = arcpy.CopyFeatures_management(connectList[0], outConnect)
    else:
           arcpy.AddWarning('No connected linkages between any patch pairs were found...')
           arcpy.AddWarning('Landscape connectivity may be compromised for this species...')
           arcpy.AddWarning('Priorities are based on restoration needed to restore connectivity to the landscape...')
           arcpy.SelectLayerByAttribute_management("linkLayer", "CLEAR_SELECTION", "" )
           outConnect = arcpy.CopyFeatures_management("linkLayer", outConnect)
           
    arcpy.Delete_management("in_memory")
    return (outConnect)

def extractPriorities(link, curSurface, LCPlayer, fullBuild):
    mosList = []

    ########
    #workspace = 'C:\\workspace\\test'
    arcpy.AddMessage('Analyzing priority linkage areas...')
    arcpy.AddMessage('\tExtracting fuzzy member set of linkage bottlenecks...')
    fuzzyCurrent = arcpy.sa.FuzzyMembership ( curSurface, arcpy.sa.FuzzyMSLarge(0.9, 1))
    # curExtract = arcpy.sa.SetNull (curSurface, curSurface, 'Value = 0')
    fuzzySTD = arcpy.sa.FocalStatistics(curSurface,arcpy.sa.NbrCircle(1000, "MAP"),"STD", "DATA")
    #fuzzySTD = arcpy.sa.FocalStatistics(curSurface,arcpy.sa.NbrRectangle(3, 3, "CELL"),"STD", "DATA")
    fuzzySTD = arcpy.sa.FuzzyMembership(fuzzySTD, arcpy.sa.FuzzyLarge("", 1))
    #fuzzyCurrent.save(os.path.join('c:\\workspace\\fuzzy', 'fuzzcur.img'))
    arcpy.AddMessage('\tExtracting fuzzy member set of linkage areas vulnerable to loss...')
    fuzzyLink = arcpy.sa.GreaterThan (link, fullBuild)
    fuzzyLink = arcpy.sa.Float(fuzzyLink)
    #fuzzyLink.save(os.path.join('c:\\workspace\\fuzzy', 'fuzzlink.img'))
    
##    fuzzyLink = arcpy.sa.FuzzyMembership ( lost_link, arcpy.sa.FuzzyLarge(0.5, 5))
    values = getCombinations(LCPlayer, 'PATH_RNK')[0]
    LCPlayer = arcpy.MakeFeatureLayer_management(LCPlayer, "lcpLayer")
    
    arcpy.AddMessage('\tExtracting fuzzy member set of linkage path regions...')
    for v in values:
        if v == 1:
            string = 'Primary Linkage Path...'
        elif v==2:
            string = 'Secondary Linkage Paths...'
        elif v==3:
            string = 'Potential Restoration Paths...'
        arcpy.AddMessage('\t\tProcessing ' + string)
        pName = 'prior_' + str(int(v))
        #pName = os.path.join(workspace, 'prior_' + str(v))
        arcpy.SelectLayerByAttribute_management("lcpLayer", "NEW_SELECTION", "PATH_RNK = " + str(v) )
        cd = arcpy.sa.CostDistance(LCPlayer, curSurface )
        fuzzyDistance = arcpy.sa.FuzzyMembership (cd, arcpy.sa.FuzzyMSSmall(0.1, 1), 'VERY')
        #fuzzyDistance.save(os.path.join('c:\\workspace\\fuzzy', 'fuzzdist_' + str(int(v)) + '.img'))
        overlay = arcpy.sa.FuzzyOverlay ([fuzzyCurrent, fuzzyDistance, fuzzySTD], 'AND')
        overlay = arcpy.sa.Slice (overlay, 2, 'EQUAL_INTERVAL')

        # overlay.save(os.path.join(workspace, 'fuzzyoverlay_' + str(int(v)) + '.img'))
        # Do vulnerable areas need to be enumerated separate from bottlenecks??????
        overlay2 = arcpy.sa.FuzzyOverlay ([fuzzyLink, fuzzyDistance, fuzzyCurrent, fuzzySTD], 'AND')
        #overlay2.save(pName + '_2.img')
        overlay2 = arcpy.sa.Slice (overlay2, 3, 'EQUAL_INTERVAL')
        overlay = arcpy.sa.Con(overlay2, 3, overlay, 'Value >= 2')
        #overlay = arcpy.sa.FuzzyOverlay ([overlay, overlay2], 'OR')
        overlay = arcpy.sa.ExtractByAttributes (overlay, 'Value >= 2')
        #overlay = arcpy.sa.Reclassify (overlay, 'Value', arcpy.sa.RemapValue([[3, int(v)]] ))
        overlay = arcpy.sa.Reclassify (overlay, 'Value', arcpy.sa.RemapValue([[2, int(v)],[3, int(v)]] ))
        

        #overlay = arcpy.sa.Con(overlay, int(v))
        arcpy.MakeRasterLayer_management(overlay, pName)
        #overlay.save(pName)
        mosList.append(pName)

    #Add Cleanup routines.  Memory could get full!!!!!!
    
    arcpy.SelectLayerByAttribute_management(LCPlayer, "CLEAR_SELECTION", "" )
    arcpy.MosaicToNewRaster_management(mosList, workspace, "bneck_mos", "", "", "", "1", "MINIMUM")
    arcpy.CalculateStatistics_management(os.path.join(workspace, "bneck_mos"))

    ### add filters here...........................
    arcpy.AddMessage('\tCleaning priority patch boundaries...')
    b_neck = arcpy.sa.BoundaryClean (os.path.join(workspace, "bneck_mos"), "ASCEND", "TWO_WAY")
    # b_neck.save('c:\\workspace\\test\\b_neck_clean.img')
    arcpy.AddMessage('\tExtracting restoration areas...')
    r_patch = arcpy.sa.BooleanAnd (b_neck, link)
    #r_patch.save('c:\\workspace\\test\\r_patch_bool.img')
    r_patch = arcpy.sa.ExtractByAttributes (r_patch, 'Value = 0')
##    r_patch.save('c:\\workspace\\test\\r_patch_extract.img')
    arcpy.AddMessage('\tExtracting areas vulnerable to loss...')
    v_patch = arcpy.sa.GreaterThan (link, fullBuild)
##    v_patch.save('c:\\workspace\\test\\v_patch_GT.img')
    #v_patch = arcpy.sa.Con(v_patch, b_neck, "", 'Value = 1')
    #v_patch.save('c:\\workspace\\test\\v_patch_Con.img')
    v_patch = arcpy.sa.ExtractByAttributes (v_patch, 'Value = 1')
##    v_patch.save('c:\\workspace\\test\\v_patch_extract.img')
    arcpy.AddMessage('\tExtracting secure areas...')
    b_patch = arcpy.sa.Reclassify (b_neck, 'Value', arcpy.sa.RemapRange([[1, 3, 2]] ))

    ###
##    b_patch = arcpy.sa.Times(b_patch, 1)
##    b_patch.save('c:\\workspace\\test\\b_patch_reclass.img')
    ###
##    workspace2 = 'c:\\workspace\\test'
##    b_neck.save(os.path.join(workspace2, 'bneck_mos.img'))

    arcpy.MosaicToNewRaster_management([r_patch, v_patch, b_patch], workspace, "prior_mos", "", "", "", "1", "MINIMUM")
    
    
    p_patch = arcpy.RasterToPolygon_conversion (os.path.join(workspace, "prior_mos"), 'in_memory\\p_patch', 0)
    b_neck = arcpy.RasterToPolygon_conversion (b_neck, 'in_memory\\b_neck', 0)
    
    try:
        arcpy.Delete_management(r_patch)
        arcpy.Delete_management(v_patch)
        arcpy.Delete_management(b_patch)
    except:
        pass
    
    outP = arcpy.Intersect_analysis ([b_neck, p_patch], 'in_memory\\outP')


    arcpy.AddMessage('\tFormatting priority attributes...')

##    fields= arcpy.ListFields(outP)
##    for field in fields:
##        arcpy.AddWarning('field.name: ' + field.name)

    

    table = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'MiscFiles', 'priorityjointable.dbf'))
    arcpy.JoinField_management(outP, 'grid_code', table, "PATH_RANK", "PATH_TYPE")
    arcpy.JoinField_management(outP, 'grid_code_1', table, "PATCH_RANK", "PATCH_TYPE")
    arcpy.AddField_management(outP, "rnk", "SHORT")
    arcpy.CalculateField_management (outP, "rnk", "int(str(!grid_code!) + str( !grid_code_1!))", "PYTHON", "")
    if len(connectList) > 0:
        arcpy.JoinField_management(outP, 'rnk', table, "RNK", "PRIORITY_R")
    else:
        arcpy.JoinField_management(outP, 'rnk', table, "RNK", "REST_RNK")

##    arcpy.AddField_management ('pLayer', 'PATH_TYPE', 'TEXT', "", "", 15)
##    arcpy.AddField_management ('pLayer', 'PATCH_TYPE', 'TEXT', "", "", 15)

    # Get the fields from the input
    fields= arcpy.ListFields(outP)

    # Create a fieldinfo object
    fieldinfo = arcpy.FieldInfo()

    # Iterate through the fields and set them to fieldinfo
    for field in fields:
        if field.name == "grid_code":
            fieldinfo.addField(field.name, "PATH_RANK", "VISIBLE", "")
        elif field.name == "PATH_TYPE":
            fieldinfo.addField(field.name, "PATH_TYPE", "VISIBLE", "")
        elif field.name == "grid_code_1":
            fieldinfo.addField(field.name, "PATCH_RANK", "VISIBLE", "")
        elif field.name == "PATCH_TYPE":
            fieldinfo.addField(field.name, "PATCH_TYPE", "VISIBLE", "")
        elif field.name == "PRIORITY_R" or field.name == "REST_RNK":
            fieldinfo.addField(field.name, "PRIOR_RNK", "VISIBLE", "")
        else:
            fieldInfo = fieldinfo.addField(field.name, field.name, "HIDDEN", "")
        
    arcpy.MakeFeatureLayer_management (outP, 'pLayer', "", "", fieldinfo)

    arcpy.AddField_management("pLayer", "ACRES", "DOUBLE", "", "", "", "", "", "", "")

    arcpy.AddField_management("pLayer", "ACRES", "DOUBLE", "", "", "", "", "", "", "")
    arcpy.CalculateField_management ("pLayer", "ACRES", "!shape.area@acres!", "PYTHON", "")
    arcpy.SelectLayerByAttribute_management("pLayer", "NEW_SELECTION", "ACRES < 5" )
    arcpy.Eliminate_management ("pLayer", "in_memory\\pLayer", "AREA")
    arcpy.Delete_management("pLayer")
    arcpy.MakeFeatureLayer_management ("in_memory\\pLayer", 'outpLayer')
    arcpy.SelectLayerByAttribute_management("outpLayer", "NEW_SELECTION", "ACRES >= 5" )
    
    return('outpLayer')
    

# Rank Paths ranks pathways connecting source patches using cost path analysis
def rankPaths(source, pField, curSurface, outConnect):
    arcpy.AddMessage('Generating ranked cost paths for ' + outConnect + '...')

    cList = []
    zList = []
    rList = []

    # Set intersect tolerance to 3X link layer cell size to prevent Intersect from creating multiple line segments where slivers occur
    interTol = str(3 * int(arcpy.Describe(link).meanCellWidth))

    # If connected region is not empty, extract cost surface by connected region to limit analysis to connected region
    if len(connectList) > 0:
        curSurface2 = arcpy.CopyRaster_management(curSurface, "curSurface2")
        arcpy.AddMessage('Extracting cost surface by connected area...')
        curSurface = arcpy.gp.ExtractByMask_sa(curSurface, outConnect, "cstSurf")
        curSurface = arcpy.Describe(curSurface).name
        curSurface2 = arcpy.Describe(curSurface2).name

    # Create line segment where source patches touch connected regions to use as sources for cost paths

    # Make sure inputs are in same projection

    sourceProjName = arcpy.Describe(source).spatialreference.name
    curProjName = arcpy.Describe(curSurface).spatialreference.name

    if not sourceProjName == curProjName:
        arcpy.AddMessage("\tReprojecting source layer...")
        pSource = arcpy.Project_management (source, os.path.join(arcpy.env.scratchWorkspace, "reproj.shp"), curSurface)
    else:
        pSource = source

    
    arcpy.AddMessage('\tIntersecting source patches with connected area to create source regions...')
    pSource = arcpy.EliminatePolygonPart_management(pSource, "in_memory\\eliminate", "PERCENT", "", 10, "CONTAINED_ONLY")
    try:
        arcpy.Delete_management(os.path.join(arcpy.env.scratchWorkspace, "reproj.shp"))
    except:
        pass
    pSource = arcpy.Intersect_analysis([[pSource, 1], [outConnect, 1]], "in_memory\\intersect", "ALL", interTol, "LINE")
    pSource = arcpy.MultipartToSinglepart_management(pSource, "in_memory\\multipart")
    pSource = arcpy.UnsplitLine_management(pSource, "in_memory\\unsplit", pField)
    pSource = arcpy.MakeFeatureLayer_management(pSource, "pSource")
    # arcpy.CopyFeatures_management(pSource, 'c:\\workspace\\test\\pSource.shp')

    # Calculate least-cost path for each pair-wise combination of source patches
    l = getCombinations(source,pField)
    values = l[0]
    combs = l[1]
    
    # break combination and not connected lists into unique elements and create list of regions with no connections
    if len(connectList) > 0:
        theList = connectList
    else:
        theList = noConnectList
        
    c = list(set(chain.from_iterable(theList)))
    
    
    # Create patch regions and cost distance rasters for each unique value in source patches
    arcpy.AddMessage('\tCreating patch regions and cost distance rasters for each unique value in source patches...')
    for v in values:
        if v  in c:
            v = str(int(v))
            arcpy.AddMessage('\t\tProcessing patch region ' + v + '...')
            arcpy.SelectLayerByAttribute_management(pSource, "NEW_SELECTION", pField + " = " + v )
            arcpy.MakeFeatureLayer_management(pSource, "p_" + v)
            cd = arcpy.sa.CostDistance("p_" + v, curSurface, "", os.path.join(workspace, "bklnk_" + v ))
            arcpy.MakeRasterLayer_management(cd, "CostDist_" + v)
            
            if len(connectList) > 0:
                rd = arcpy.sa.CostDistance("p_" + v, curSurface2, "", os.path.join(workspace, "r_bklnk_" + v )) 
                arcpy.MakeRasterLayer_management(rd, "r_CostDist_" + v)

    # Create least-cost paths for each region pair in both directions
    arcpy.AddMessage('\tGenerating least-cost path for each patch pair combination...')
    
    for c in combs:
        c1 = str(int(c[0]))
        c2 = str(int(c[1]))
        
        if  c in theList:
            arcpy.AddMessage('\t\tCalculating least-cost path from region ' + c1 + ' to region ' + c2 + '...')
            cp = arcpy.sa.CostPath("p_" + c1, "CostDist_" + c2, "bklnk_" + c2, "BEST_SINGLE", "FID")
            cp1 = arcpy.MakeRasterLayer_management(cp, "CP_" + c1 + "_" + c2)
            arcpy.AddMessage('\t\tCalculating least-cost path from region ' + c2 + ' to region ' + c1 + '...')
            cp = arcpy.sa.CostPath("p_" + c2, "CostDist_" + c1, "bklnk_" + c1, "BEST_SINGLE", "FID")
            cp2 = arcpy.MakeRasterLayer_management(cp, "CP_" + c2 + "_" + c1)

            cList.append(str(cp1))
            cList.append(str(cp2))

        else:
            arcpy.AddWarning('\t\tRegions ' + c1 + ' and ' + c2 + ' are not connected.  Skipping cost path for this region pair...')

    # Create combined least-cost path polyline layer
    arcpy.AddMessage('\t\tMosaicing least-cost paths for region pairs...')
    arcpy.MosaicToNewRaster_management(cList, workspace, "lcp_mos", "", "", "", "1", "MAXIMUM")

    for c in cList:
        try:
            arcpy.Delete_management(c)
        except:
            pass
    
    arcpy.CalculateStatistics_management(os.path.join(workspace, "lcp_mos"))
    LCP = arcpy.sa.Con(os.path.join(workspace, "lcp_mos"), "1", "", "VALUE > 0")

    arcpy.Delete_management(os.path.join(workspace, "lcp_mos"))

    # Create least-cost paths by zone
    arcpy.AddMessage('\tGenerating least-cost paths  by zones for each patch pair combination...')
    # Create least-cost paths for each region pair in both directions   
    for c in combs:
        c1 = str(int(c[0]))
        c2 = str(int(c[1]))
        if  c in theList:
            arcpy.AddMessage('\t\tCalculating least-cost path from region ' + c1 + ' to region ' + c2 + '...')
            zp = arcpy.sa.CostPath("p_" + c1, "CostDist_" + c2, "bklnk_" + c2, "EACH_ZONE", "FID")
            zp1 = arcpy.MakeRasterLayer_management(zp, "ZP_" + c1 + "_" + c2)
            arcpy.AddMessage('\t\tCalculating least-cost path from region ' + c2 + ' to region ' + c1 + '...')
            zp = arcpy.sa.CostPath("p_" + c2, "CostDist_" + c1, "bklnk_" + c1, "EACH_ZONE", "FID")
            zp2 = arcpy.MakeRasterLayer_management(zp, "ZP_" + c2 + "_" + c1)

            zList.append(str(zp1))
            zList.append(str(zp2))

    # Create combined least-cost path polyline layer
    arcpy.AddMessage('\t\tMosaicing least-cost paths for region zones...')
    if arcpy.Exists(os.path.join(workspace, "zcp_mos")):
        arcpy.Delete_management(os.path.join(workspace, "zcp_mos"))
    arcpy.MosaicToNewRaster_management(zList, workspace, "zcp_mos", "", "", "", "1", "MAXIMUM")

    for z in zList:
        try:
            arcpy.Delete_management(z)
        except:
            pass
    
    arcpy.CalculateStatistics_management(os.path.join(workspace, "zcp_mos" ))
    ZCP = arcpy.sa.Con(os.path.join(workspace, "zcp_mos" ), "2", "", "VALUE > 0")

    # Create least-cost paths through compromised areas

    if len(connectList) > 0: 
        # Create patch regions and cost distance rasters for each unique value in source patches
        arcpy.AddMessage('\tCalculating costs through restoration zones...')
        
        arcpy.AddMessage('\tGenerating potential restoration paths for each patch pair combination...')
        # Create least-cost paths for each region pair in both directions   
        for c in combs:
            c1 = str(int(c[0]))
            c2 = str(int(c[1]))
            if  c in theList:
                arcpy.AddMessage('\t\tCalculating least-cost path from region ' + c1 + ' to region ' + c2 + '...')
                rp = arcpy.sa.CostPath("p_" + c1, "r_CostDist_" + c2, "r_bklnk_" + c2, "EACH_ZONE", "FID")
                rp1 = arcpy.MakeRasterLayer_management(rp, "RP_" + c1 + "_" + c2)
                arcpy.AddMessage('\t\tCalculating least-cost path from region ' + c2 + ' to region ' + c1 + '...')
                rp = arcpy.sa.CostPath("p_" + c2, "r_CostDist_" + c1, "r_bklnk_" + c1, "EACH_ZONE", "FID")
                rp2 = arcpy.MakeRasterLayer_management(rp, "RP_" + c2 + "_" + c1)

                rList.append(str(rp1))
                rList.append(str(rp2))


        # Create combined least-cost path polyline layer
        arcpy.AddMessage('\t\tMosaicing least-cost paths for region zones...')
        if arcpy.Exists(os.path.join(workspace, "rcp_mos")):
            arcpy.Delete_management(os.path.join(workspace, "rcp_mos"))
        arcpy.MosaicToNewRaster_management(rList, workspace, "rcp_mos", "", "", "", "1", "MAXIMUM")

        for r in rList:
            try:
                arcpy.Delete_management(r)
            except:
                pass
        
        arcpy.CalculateStatistics_management(os.path.join(workspace, "rcp_mos" ))
        RCP = arcpy.sa.Con(os.path.join(workspace, "rcp_mos" ), "3", "", "VALUE > 0")
        mList = [LCP, ZCP, RCP]

    else:
        mList = [LCP, ZCP]

    arcpy.AddMessage('\tCombining least-cost paths by region and least-cost paths by region zones...')
    arcpy.MosaicToNewRaster_management(mList, workspace, "lcp_mos", "", "", "", "1", "MINIMUM")
    LCP = arcpy.RasterToPolyline_conversion(os.path.join(workspace, "lcp_mos"), "LCP", "", "", "NO_SIMPLIFY")

    # Create a fieldinfo object to rename grid_code field
    fieldinfo = arcpy.FieldInfo()
    fieldinfo.addField("GRID_CODE", "PATH_RNK", "VISIBLE", "")
    outLCP = arcpy.MakeFeatureLayer_management(str(LCP), "outLCP", "", "", fieldinfo)
    # arcpy.CopyFeatures_management(outLCP, os.path.join(workspace, outLCP.shp))

    try:
        arcpy.Delete_management(os.path.join(workspace, "lcp_mos"))
        arcpy.Delete_management(os.path.join(workspace, "zcp_mos"))
        arcpy.Delete_management(os.path.join(workspace, "rcp_mos"))
    except:
        pass
    # arcpy.Delete_management("in_memory")
    return (outLCP)

########################### End Functions #####################################

#### NOTE: len(connectList) can be used to test if landscape is connected

outConnect = str(findConnected(link,source,pField,outConnect))

# Convert noConnectList to list of tuples of patch pairs
connectList = [(float(c.split('_')[1]), float(c.split('_')[2])) for c in connectList]
noConnectList = [(float(c.split('_')[1]), float(c.split('_')[2])) for c in noConnectList]

# Temporarily set workspace to scratchWorkspace
oldWorkspace = arcpy.env.workspace

# If ArcGIS version 10.1 or greater, use in_memory workspace for scratch outputs (10.1 supports raster in_memory)
version = arcpy.GetInstallInfo('DESKTOP')
if float(version['Version']) >= 10.1:
    arcpy.env.workspace = 'in_memory'
    workspace = 'in_memory'
else:
    arcpy.env.workspace = arcpy.env.scratchWorkspace
    workspace = arcpy.env.scratchWorkspace

# Create LCPs
invSurface = invertRaster(curSurface, "curSurface")
paths = rankPaths(source, pField, invSurface, outConnect)
outLCP = arcpy.CopyFeatures_management(paths, outLCP )

# Extract priorities

pAreas = extractPriorities(link, curSurface, outLCP, fullBuild)
outPAreas = arcpy.CopyFeatures_management (pAreas, outPAreas)

#Add outputs to display
arcpy.SetParameterAsText(6, outConnect)

params = arcpy.GetParameterInfo()

if len(connectList) > 0:
    params[7].symbology = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'MiscFiles', 'path_rnk.lyr'))
else:
    params[7].symbology = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'MiscFiles', 'rest_path_rnk.lyr'))

arcpy.SetParameterAsText(7, outLCP)

arcpy.SetParameterAsText(8, outPAreas)

arcpy.AddMessage('Restoring current workspace to original location...')
arcpy.env.workspace = oldWorkspace

arcpy.AddMessage('Deleting temporary files...')   
arcpy.Delete_management("in_memory")

