#!/usr/bin/env python
# version 1.0

__author__ = "Bruce Wollison"
__date__ = "Nov, 2017"

'''
    Description: This script contains functions that write panda
    tables to an excel file. We want to abstract out all excel calls
    to a simple format to minimize how much each script needs to do 
    to write to the excel sheet. As such, we'll have some simple rules
    for how reports are written:
    
    1) all tables have a bolded header, unbolded contents
    2) all written contents have a border of one blank row between them when written to the
            same sheet
    3) all contents will be written from the top to the bottom of the sheet,
            and never from right to left
    These rules should allow us to simplify the API by only adding new sheets, 
    then adding new contents (which would be either pictures or tables) in a 
    sequential organized fashion
    
    Possible issues with these requirements: It prohibits an obvious way to implement notes
    which may be desirable for explaining metrics. A work around for this might be to allow 
    inserting tables to be accompanies by a list or dict object with the comments for the
    header row.
    
    API functions should be include
    a) XLwrite - write to a specific file
    b) newSheet - create a new sheet (with sheet name)
    c) writeTable - insert text box (into current 'cursor' position)
    d) writeImage - insert picture (into current 'cursor' position)
    e) save - save the data to the final excel file
'''

# library imports
import math
import openpyxl as xl
from openpyxl.utils.dataframe import dataframe_to_rows

### wrapping functions from the XL API here so we can keep them abstracted
### as opposed to having an implementation that is specific to an API
### this way we can reimplement with a diff API or in a diff output format with
### hopefully minimal, if any, changing of the other scripts
### this is a lesson learned from using XL APIs in R, it was sometimes necessary
### to use different ones because none are official and all seemed sub-optimal.
### in the future we may want to look at others here:
### https://github.com/kz26/PyExcelerate

# translates a column number that we might get from pandas to
# the column letter format Excel uses. EX: for the first 26 we use a single
# letter then add another, etc.
# we are only going to have this function good for columns 1 through ~700
# if that isn't enough, tough luck, you're probably transposing your table
# wrong anyway
def colNumToExcelColLetter(colNum):
    alphaReps = colNum/26
    if alphaReps > 1:
        # make 2 chars based on Excel formatting, which is A-Z then
        # AA-AZ then BA-BZ etc
        firstChar = chr(math.floor(alphaReps)+64)
        secondChar = chr(colNum-(math.floor(alphaReps)*26) + 65)
        return firstChar + secondChar
    else:
        # note, we are ord by 65 instead of 64 because
        # the colNum is 0 indexed, so we need 0 to be 'A'
        return chr(colNum+65)


class XLwrite:

    # begin a new workbook. Requires 2 strings:
    # XLFileName - string for the name of the excel file to write to
    # XLSheetName - string for the name of the first sheet to begin writing to
    def __init__(self, XLFileName, XLSheetName="Sheet1"):
        self.workbook = xl.Workbook()
        self.filename = XLFileName
        self.currSheet = self.workbook.active
        self.currSheet.title = XLSheetName
        # try and guess the data types of the data we write
        # this prevents numbers from being strings, which is really annoying
        self.workbook.guess_types = True
        # firstSheetElement lets us know if the active sheet has has something
        # written to it yet. If it is False we want to add a blank row
        self.firstSheetElement = True


    # creates a new sheet to begin writing to
    # NewSheetName - string for the name of the new sheet
    def newSheet(self, NewSheetName="New_Sheet"):
        self.currSheet = self.workbook.create_sheet(title=NewSheetName)
        self.firstSheetElement = True

    # sets the style for a column, takes in a column letter
    # and a range of numbers for the rows and sets them to the
    # specified style for the current worksheet
    def setColStyle(self, column, startRow, endRow, style):
        for currRow in range(startRow, endRow):
            cellName = column + str(currRow)
            self.currSheet[cellName].number_format = style


    # writes a table to the current sheet. Skips a row if a table was previously written.
    # pandaTable - a panda data frame to write to the excel file
    # styles shoudl be a dict of {column name : style} keypairs
    def writeTable(self, pandaTable, styles={}):
        firstTableRow = True
        # openpyxl default behavior is to write in the immediate row after the previous one
        # so we need to add blank rows if an items already been written to the sheet
        if not self.firstSheetElement:
            self.currSheet.append([])
        self.firstSheetElement = False
        # loop over the data frame and write it. Openpyxl chokes hard on a lot of pandas
        # data types, so we want to also convert all of them to strings first so that
        # it can deal with them properly <.applymap(str)>.
        for r in dataframe_to_rows(pandaTable.applymap(str), index=False, header=True):
            if not firstTableRow:
                self.currSheet.append(r)
            else:
                self.currSheet.append(r)
                # here we need to apply the header style to the row that we write
                # first find the row we are on then apply the style
                currRow = self.currSheet._current_row
                rowToAlter = self.currSheet.row_dimensions[int(currRow)]
                rowToAlter.font = xl.styles.Font(bold=True, underline='single')
                firstTableRow = False
        # apply styles to columns where necessary, we don't do it if no styles
        # are supplied or the table is empty
        if not styles == {} and not pandaTable.empty:
            # retrieve necessary row values so we know where to start and stop setting
            # the column styles for this table.
            bottomRowInd = self.currSheet._current_row + 1
            topRowInd = bottomRowInd - len(pandaTable.index)
            # check for matching col names between the pandaTable and
            # in supplies styles dict, then set the column style
            for colMatch, styleToSetColTo in styles.iteritems():
                if colMatch in list(pandaTable.columns.values):
                    colLetter = colNumToExcelColLetter(pandaTable.columns.get_loc(colMatch))
                    self.setColStyle(colLetter, topRowInd, bottomRowInd, styleToSetColTo)


    # writes a table to the current sheet. Skips a row if a table was previously written.
    # pandaTable - a panda data frame to write to the excel file
    def writeImage(self, imageFilePath):
        # openpyxl default behavior is to write in the immediate row after the previous one
        # so we need to add blank rows if an items already been written to the sheet
        if not self.firstSheetElement:
            self.currSheet.append([])
        # write the image to the sheet
        img = xl.drawing.image.Image(imageFilePath)
        self.add_image(img)
        self.firstSheetElement = False


    # saves the contents of the workbook
    def save(self):
        self.workbook.save(self.filename)
        print("wrote data to file " + self.filename)

# EOF