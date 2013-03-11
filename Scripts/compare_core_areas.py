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
#compare_core_area.py Parses all feature class tables created by 'Evaluate Habitat Patches' in a directory into
#single table. Can optionally process all subdirectories of input folder.
# ------------------------------------------------------------------------------

import sys, os, arcpy

## Finds a suitable template file to create core summary table from
def FindTemplate (l, fields):
    template = ''
    for fc in l:
        f = [x.name for x in arcpy.ListFields(fc)]
        if len(set(f).intersection(fields)) == 7:
            arcpy.AddMessage('\tUsing table: ' + fc + ' as template...')
            template = fc
            break
    return template

def NoPoly(template):
    if template == '':
        arcpy.AddError('A suitable table template was not found. Make sure input folder contains at least one output from Evaluate Habitat Patches')
        raise

p = sys.argv[1]
searchSub = sys.argv[2]
outTable = sys.argv[3]

oldWorkspace = arcpy.env.workspace
arcpy.env.workspace = p

fields = [u'HOUSES', u'PRCNT_AREA', u'AREA', u'TOT_AREA', u'PRCNT_TOT', u'AREA_UNITS', u'CORE']

if searchSub== 'true':
    subFolders = [x[0] for x in os.walk(p)]
else:
    subFolders = p.split(',')

# Find a table to use as a template for output
arcpy.AddMessage('Finding table to use as template for output...')
for d in subFolders:
    arcpy.env.workspace = d
    l = arcpy.ListFeatureClasses('', 'Polygon')
    if len(l) > 0:
        template = FindTemplate(l, fields)
    if not str(template) == '':
        break

NoPoly(template)

##try:
##    if template == '':
##        raise NoTemplate()
##except NoTemplate:        
##        arcpy.AddError('A suitable table template was not found. Make sure input folder contains at least one output from Evaluate Habitat Patches')
##        sys.exit

# Make Table to receive data
arcpy.AddMessage('Creating new output table...')
template = arcpy.MakeTableView_management(template, 'template')
sumTable = os.path.join(arcpy.env.scratchWorkspace, 'xxsumTable.dbf')

if arcpy.Exists(sumTable):
    arcpy.Delete_management(sumTable)


arcpy.CreateTable_management(os.path.dirname(sumTable), os.path.basename(sumTable), template)
arcpy.AddField_management(sumTable, "FileName", "TEXT", "", "", "50")
arcpy.AddField_management(sumTable, "CORE_AREA", "DOUBLE", "", "", "")
arcpy.AddField_management(sumTable, "PATCH_AREA", "DOUBLE", "", "", "")
oTable = arcpy.MakeTableView_management(sumTable, 'sumTable')

arcpy.AddMessage('Appending data records to output table...')
for folder in subFolders:
    if searchSub == 'true':
        arcpy.AddMessage('Search Subfolders option selected. Processing directory: ' + folder)
        arcpy.env.workspace = folder
        l = arcpy.ListFeatureClasses()
    for fc in l:
        f = [x.name for x in arcpy.ListFields(fc)]
        if len(set(f).intersection(fields)) == 7:
            arcpy.AddMessage('\tProcessing table: ' + fc + '...')
            tempTable = arcpy.TableToTable_conversion(fc, arcpy.env.scratchWorkspace, 'xxtempTable.dbf')
                    
            table = arcpy.MakeTableView_management(tempTable, "table")
            arcpy.AddField_management(sumTable, "FileName", "TEXT", "", "", "50")
            

            arcpy.AddField_management(table, "CORE_AREA", "DOUBLE", "", "", "")
            arcpy.AddField_management(table, "PATCH_AREA", "DOUBLE", "", "", "")
            arcpy.SelectLayerByAttribute_management("table", "NEW_SELECTION", "CORE = 1" )
            arcpy.CalculateField_management("table", "CORE_AREA", "!AREA!", "PYTHON")
            arcpy.SelectLayerByAttribute_management("table", "NEW_SELECTION", "CORE = 0" )
            arcpy.CalculateField_management("table", "PATCH_AREA", "!AREA!", "PYTHON")
            arcpy.SelectLayerByAttribute_management("table", "CLEAR_SELECTION" )

            # arcpy.SelectLayerByAttribute_management("table", "NEW_SELECTION", ' "CORE" = 1 ')
            arcpy.Append_management(table, oTable, "NO_TEST","","")
            arcpy.SelectLayerByAttribute_management("sumTable", "NEW_SELECTION", '"' + 'FileName' + '"' + ' = ' + "''" )
            fileName = '"' + os.path.splitext(fc)[0] + '"'
            arcpy.CalculateField_management("sumTable", "FileName", fileName)
            arcpy.SelectLayerByAttribute_management("sumTable", "CLEAR_SELECTION" )
            
            del table
            #del sumTable
            arcpy.Delete_management(tempTable)
            # arcpy.Delete_management(oTable)

#Process output table for final output

arcpy.AddMessage('Generating summary statistics...')            
statTable = os.path.join(arcpy.env.scratchWorkspace, 'xxstatTable.dbf')
arcpy.Statistics_analysis(sumTable, statTable, [["HOUSES", "SUM"], ["CORE_AREA", "MAX"], ["PATCH_AREA", "MAX"], ["TOT_AREA", "MAX"]], ["FileName"])
sTable = arcpy.MakeTableView_management(statTable, 'statTable')

for n in ["TOT_HAB", "PCNT_CORE", "PCNT_PATCH", "PCNT_HAB"]:
    arcpy.AddMessage('\tCalculating ' + n + '...')
    arcpy.AddField_management("statTable", n, "DOUBLE", "", "", "")

arcpy.CalculateField_management("statTable", "TOT_HAB", "!MAX_CORE_A! + !MAX_PATCH_!", "PYTHON")
arcpy.CalculateField_management("statTable", "PCNT_CORE", "(!MAX_CORE_A! / !MAX_TOT_AR!) * 100", "PYTHON")
arcpy.CalculateField_management("statTable", "PCNT_PATCH", "(!MAX_PATCH_! / !MAX_TOT_AR!) * 100", "PYTHON")
arcpy.CalculateField_management("statTable", "PCNT_HAB", "(!TOT_HAB! / !MAX_TOT_AR!) * 100", "PYTHON")

arcpy.AddMessage('Formatting summary table...')

# Get the fields from the input
fields= arcpy.ListFields(statTable)

# Create a fieldinfo object
fieldinfo = arcpy.FieldInfo()

# Iterate through the fields and set them to fieldinfo
for field in fields:
    if field.name == "SUM_HOUSES":
        fieldinfo.addField(field.name, "HOUSES", "VISIBLE", "")
    elif field.name == "MAX_CORE_A":
        fieldinfo.addField(field.name, "CORE_AREA", "VISIBLE", "")
    elif field.name == "MAX_PATCH_":
        fieldinfo.addField(field.name, "PATCH_AREA", "VISIBLE", "")
    elif field.name == "MAX_TOT_AR":
        fieldinfo.addField(field.name, "TOT_AREA", "VISIBLE", "")
        
# Rename fields and copy rows to output table
arcpy.MakeTableView_management('statTable', 'outTable', "", "", fieldinfo)
arcpy.CopyRows_management('outTable', outTable)

arcpy.AddMessage('Restoring current workspace to original location...')
arcpy.env.workspace = oldWorkspace
arcpy.Delete_management(sumTable)
arcpy.Delete_management(statTable)


