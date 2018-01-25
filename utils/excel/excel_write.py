#!/usr/bin/env python
# version 1.0
# tested with python 2.7.13

__author__ = "Bruce Wollison"
__date__ = "Nov, 2017"

'''
    Description: This script contains functions that write panda
    tables to an excel file. We want to abstract out all excel calls
    to a simple format to minimize how much each script needs to do 
    to write to the excel sheet. As such, we'll have some simple rules
    for how reports are written:
    
    1) all tables have a bolded header, unbolded contents
    2) all written contents have one blank row between them when written to the
            same sheet
    3) all contents will be written from the top to the bottom of the sheet,
            and never from right to leftct sc
    These rules should allow us to simplify the API by only adding new sheets, 
    then adding new contents (which would be either pictures or tables) in a 
    sequential organized fashion
    
    Possible issues with these requirements: It prohibits an obvious way to implement excel cell comments
    which may be desirable for explaining metrics. A work around for this might be to allow 
    inserting tables to be accompanies by a list or dict object with the comments for the
    header row. Right now we are not implementing this feature, for the variant report at least we
    are currently writing the column descriptions into their own tab.
    
    Currently table styles are implemented with dict objects, with matching key:value 
    being used to represent column:style pairs. 
    
    API functions should be include
    a) XLwrite - write to a specific file
    b) newSheet - create a new sheet (with sheet name)
    c) writeTable - insert text box (into current 'cursor' position)
    d) writeImage - insert picture (into current 'cursor' position)
    e) save - save the data to the final excel file
'''

# library imports
import math, logging
import openpyxl as xl
from openpyxl.utils.dataframe import dataframe_to_rows

"""
bmw35
wrapping functions from the XL API here so we can keep them abstracted
as opposed to having an implementation that is specific to an API
this way we can reimplement with a diff API or in a diff output format with
hopefully minimal, if any, changing of the other scripts
this is a lesson learned from using XL APIs in R, it was sometimes necessary
to use different ones because none are official and all seemed sub-optimal.
in the future we may want to look at others here (some outlined here):
https://github.com/kz26/PyExcelerate
"""

def colNumToExcelColLetter(colNum):
    """
    bmw35 12-28-2017
    translates a column number that we might get from pandas to
    the column letter format Excel uses. EX: for the first 26 we use a single
    letter then add another, etc.
    we are only going to have this function good for columns 1 through ~700
    if that isn't enough, tough luck, you're probably transposing your table
    wrong anyway
    :param colNum: int value representing the column number to convert to an Excel
      column name
    :return: string, returns 1-2 uppercase letters
    """
    ALPHABET_LENGTH = 26
    alphaReps = colNum/ALPHABET_LENGTH
    if alphaReps > 1:
        # make 2 chars based on Excel formatting, which is A-Z then
        # AA-AZ then BA-BZ etc. No [A-Z][A-Z][A-Z] support.
        firstChar = chr(math.floor(alphaReps)+64)
        secondChar = chr(colNum-(math.floor(alphaReps)*ALPHABET_LENGTH) + 65)
        return firstChar + secondChar
    else:
        # note, we are ord by 65 instead of 64 because
        # the colNum is 0 indexed, so we need 0 to be 'A'
        return chr(colNum+65)


class XLwrite:
    """
    bmw35 01-16-2018
    class XLwrite abstracts the excel writing package used
    to a few basic functions for writing the outputs for reporting
    scripts. Wrapping the package into a module class here prevents
    rewriting every script if we need to change out the python excel
    package at any point
    """

    # globals for this class
    # starting with accepted number styles:
    NUMBER_STYLES = ["0.00%"]

    def __init__(self, XLFileName, XLSheetName="Sheet1"):
        """
        bmw35 12-28-2017
        Constructor. Creates a new workbook object.
        :param XLFileName: string for the name of the excel file to eventually write the workbook to
        :param XLSheetName: string for the name of the first sheet to begin writing to
        :return: XLwrite object new XLwrite object that contains a new workbook with active sheet
        """
        # set the logger for this XL writer object
        # setting it to the same logger as the script logger
        # automagically has it output to the main scripts logging
        # paramters
        self.logger = logging.getLogger('Logger')
        # create a new workbook and what not
        self.workbook = xl.Workbook()
        self.filename = XLFileName
        self.currSheet = self.workbook.active
        self.currSheet.title = XLSheetName
        # try and guess the data types of the data we write
        # this prevents numbers from being strings, which is really annoying
        self.workbook.guess_types = True


    def newSheet(self, NewSheetName="New_Sheet"):
        """
        bmw35 12-28-2017
        creates a new sheet to begin writing to
        :param NewSheetName: string for the name of the new sheet
        :return: adds a new worksheet to the existing workbook object
        """
        self.currSheet = self.workbook.create_sheet(title=NewSheetName)


    def setColStyle(self, column, startRow, endRow, style):
        """
        bmw35 12-28-2017
        sets the style for a column, takes in a column letter
        and a range of numbers for the rows and sets them to the
        specified style for the current worksheet
        :param column: string (1-2 chr) that corresponds to an excel column to apply the style to
        :param startRow: int indicating which row to start applying the style at
        :param endRow: int indicating which row to end applying the style at
        :param style: string that coresponds to a valid cell style
        :return: changes the cell styles for the column provided between the given rows
        """
        for currRow in range(startRow, endRow):
            cellName = column + str(currRow)
            # check that we know the style type
            if style in XLwrite.NUMBER_STYLES:
                self.currSheet[cellName].number_format = style
            else:
                self.logger.error("Error! Style " + style + " unknown!")


    def writeTable(self, pandaTable, styles={}):
        """
        bmw35 12-28-2017
        writes a table to the current sheet. Skips a row if a table was previously written.
        :param pandaTable: a panda data frame to write to the excel file
        :param styles: styles should be a dict of {column name : style} keypairs indicating which columns
           to change to what style
        :return: adds the pandas dataframe to the current worksheet in the appropriate position
        """
        firstTableRow = True
        # openpyxl default behavior is to write in the immediate row after the previous one
        # so we need to add blank rows if an items already been written to the sheet
        if self.currSheet._current_row != 0:
            self.currSheet.append([])
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


    def writeImage(self, imageFilePath):
        """
        bmw35 12-28-2017
        writes an image to the current sheet. Skips a row if anything was previously written in the worksheet
        :param imageFilePath: string of a filepath of an image to write to the excel file
        :return: adds the image at the specific filepath to the current worksheet
        """
        # openpyxl default behavior is to write in the immediate row after the previous one
        # so we need to add blank rows if an items already been written to the sheet
        if self.currSheet._current_row != 0:
            self.currSheet.append([])
        # write the image to the sheet
        img = xl.drawing.image.Image(imageFilePath)
        self.add_image(img)


    def save(self):
        """
        bmw35 12-28-2017
        saves the contents of the workbook
        :return: saves the workbook to a file
        """
        self.workbook.save(self.filename)
        self.logger.info("wrote data to file " + self.filename)

# EOF