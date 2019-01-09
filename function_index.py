# function function_index
# _
# List all Functions from the MACS Toolbox (Python)
# FORMAT function_index
#
# Note: Running this script requires Python 3 and the xlsxwriter module.
# You can install the latter by typing this at the OS command line:
#     pip install Xlsxwriter
# 
# Author: Joram Soch, BCCN Berlin
# E-Mail: joram.soch@bccn-berlin.de
# 
# First edit: 01/02/2018, 22:45 (V1.2/V18)
#  Last edit: 08/01/2019, 12:40 (V1.4/V20)


# Get toolbox file index
#-------------------------------------------------------------------------#
import os
files = [f for f in os.listdir('.') if os.path.isfile(f)]

# Open a new Excel file
#-------------------------------------------------------------------------#
import xlsxwriter
wb = xlsxwriter.Workbook('function_index.xlsx')

# Prepare new worksheet
#-------------------------------------------------------------------------#
ws = wb.add_worksheet('Functions')
bf = wb.add_format({'bold': True, 'bottom': 1, 'bottom_color': '#000000'})
ws.write(0, 0, 'Function Name', bf)          # header line
ws.write(0, 1, 'Function Type', bf)
ws.write(0, 2, 'Function Purpose', bf)
ws.write(0, 3, 'First Edit', bf)
ws.write(0, 4, 'Time', bf)
ws.write(0, 5, 'Ver', bf)
ws.write(0, 6, 'Num', bf)
ws.write(0, 7, 'Last Edit', bf)
ws.write(0, 8, 'Time', bf)
ws.write(0, 9, 'Ver', bf)
ws.write(0,10, 'Num', bf)
ws.write(0,11, 'Code Lines', bf)
ws.set_column('A:A', 25)                     # column widths
ws.set_column('B:B', 20)
ws.set_column('C:C', 70)
ws.set_column('D:D', 10)
ws.set_column('E:G',  5)
ws.set_column('H:H', 10)
ws.set_column('I:K',  5)
ws.set_column('L:L', 10)
ws.freeze_panes(1, 0)                        # fix header line

# Prepare function index
#-------------------------------------------------------------------------#
n = 0
l = 0

# Browse toolbox functions
#-------------------------------------------------------------------------#
for file in files:
    if ('.m' in file and '.md' not in file) or '.py' in file:

        # Get function type
        #-----------------------------------------------------------------#
        n = n + 1
        if file.find('batch_') == 0:
            typ = 'batch function'
        if file.find('MA_') == 0:
            typ = 'model assessment'
        if file.find('MC_') == 0:
            typ = 'model comparison'
        if file.find('MS_') == 0:
            typ = 'model selection'
        if file.find('ME_') == 0:
            typ = 'model estimation'
        if file.find('MD_') == 0:
            typ = 'many distributions'
        if file.find('MF_') == 0:
            typ = 'more functions'
        if file.find('fun') == 0:
            typ = 'meta function'

        # Read function code
        #-----------------------------------------------------------------#
        file_obj = open(file, 'r')
        file_txt = file_obj.readlines()
        purpose  = file_txt[2][2:-1]         # function purpose
        nocl     = len(file_txt)             # number of code lines
        for line in file_txt:
            if line.find('% First edit:') == 0 or line.find('# First edit:') == 0:
                first_edit = line[14:-1]     # first edit
                first_date = first_edit[0:10]
                first_time = first_edit[12:17]
                first_ver  = first_edit[19:first_edit.find('/V')]
                first_num  = first_edit[first_edit.find('/V')+1:-1]
            if line.find('%  Last edit:') == 0 or line.find('#  Last edit:') == 0:
                last_edit = line[14:-1]      # last edit
                last_date = last_edit[0:10]
                last_time = last_edit[12:17]
                last_ver  = last_edit[19:last_edit.find('/V')]
                last_num  = last_edit[last_edit.find('/V')+1:-1]


        # Store function info
        #-----------------------------------------------------------------#
        if '.m' in file:
            ws.write(n, 0, file[:-2])
        else:
            ws.write(n, 0, file)
        ws.write(n, 1, typ)
        ws.write(n, 2, purpose)
        ws.write(n, 3, first_date)
        ws.write(n, 4, first_time)
        ws.write(n, 5, first_ver)
        ws.write(n, 6, first_num)
        ws.write(n, 7, last_date)
        ws.write(n, 8, last_time)
        ws.write(n, 9, last_ver)
        ws.write(n,10, last_num)
        ws.write(n,11, nocl)
        l = l + nocl
       
# Close function index
#-------------------------------------------------------------------------#
ws.write(n+1, 2, 'Total Number of Code Lines')
ws.write(n+1,11, l)
wb.close()
