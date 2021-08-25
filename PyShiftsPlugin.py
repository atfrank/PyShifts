# Copyright Notice
# ================
#
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
#
# ----------------------------------------------------------------------
#               This PyMOL Plugin is Copyright (C) 2016 by
#            Aaron T. Frank <afrankz at umich dot edu>
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
#----------------------------------------------------------------------

# python lib
import os
import sys
import platform
import pandas as pd
if sys.version_info >= (2,4):
    import subprocess # subprocess is introduced in python 2.4
import math
import random
import numpy as np
import tempfile
from scipy import stats
from math import log
from Meter import Meter
from math import sqrt

import tkinter
from tkinter import *

#import tkSimpleDialog
#import tkMessageBox
#import tkFileDialog
#import tkColorChooser

import tkinter.simpledialog as tkSimpleDialog
import tkinter.filedialog as tkFileDialog
import tkinter.messagebox as tkMessageBox
import tkinter.colorchooser as tkColorChooser
from tkinter.filedialog import asksaveasfilename

# pymol lib
try:
    from pymol import cmd
    from pymol.cgo import *
    from pymol import stored
    from pymol import cmd
except ImportError:
    print('ERROR: pymol library cmd not found.')
    sys.exit(1)

# external libraries
try:
    import Pmw
except ImportError:
    print('ERROR: failed to import Pmw. Exit ...')
    sys.exit(1)

try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.neighbors import NearestNeighbors
    from sklearn.cluster import KMeans
    SKLEARN = True
except:
    SKLEARN = False
    print("Warning: Couldn't find sklearn")

# try loading optionally libraries
try:
    import psico.fullinit
    PSICO = True
except:
    PSICO = False
    print("Warning: Couldn't find psico")

try:
    import bme_reweight as bme
    BME = True
except:
    BME = False
    print("Warning: Couldn't find BME")

VERBOSE = False

#################
## here we go
#################
def __init__(self):
    """ PyShifts plugin for PyMol """
    self.menuBar.addmenuitem('Plugin', 'command','PyShifts', label = 'PyShifts', command = lambda s=self : PyShiftsPlugin(s))


#################
## GUI related
#################
class PyShiftsPlugin:
    def __init__(self, app):
        """ Set up user interface and initialize parameters """
        # parameters used by LARMORD
        self.parent = app.root
        self.pymol_sel     = tkinter.StringVar(self.parent)
        self.larmord_bin    = tkinter.StringVar(self.parent)
        self.larmorca_bin    = tkinter.StringVar(self.parent)
        self.larmord_para    = tkinter.StringVar(self.parent)
        self.larmord_ref    = tkinter.StringVar(self.parent)
        self.larmord_acc    = tkinter.StringVar(self.parent)
        self.cs    = tkinter.StringVar(self.parent)
        self.cs2    = tkinter.StringVar(self.parent)
        self.larmord_rlt_dict = {}
        self.larmord_error_all = {}
        self.larmord_error_carbon = {}
        self.larmord_error_proton = {}
        self.larmord_error_nitrogen = {}
        self.larmord_error_color = tkinter.StringVar(self.parent)
        self.larmord_error_scale = tkinter.DoubleVar(self.parent)
        self.larmord_error_height = tkinter.DoubleVar(self.parent)
        self.larmord_error_lsize = tkinter.IntVar(self.parent)
        self.larmord_ndisplayed = tkinter.IntVar(self.parent)
        self.min_size = tkinter.IntVar(self.parent)
        self.larmord_ndisplayed = tkinter.IntVar(self.parent)
        self.mae = {}
        self.measuredCS = {}
        self.predictedCS = []
        self.total_error = {}
        self.proton_error = {}
        self.nitrogen_error = {}
        self.carbon_error = {}
        self.larmord_erros = {}
        self.best_model_indices = []
        self.worst_model_indices = []
        self.larmord_proton_offset = tkinter.DoubleVar(self.parent)
        self.larmord_carbon_offset = tkinter.DoubleVar(self.parent)
        self.larmord_nitrogen_offset = tkinter.DoubleVar(self.parent)
        self.larmord_outlier_threshold = tkinter.DoubleVar(self.parent)
        self.get_shifts_from_larmord = True
        self.get_shifts_from_larmorca = False
        self.get_shifts_from_file = False
        self.weighted_errors = False
        self.larmord_error_sel = 'all'
        self.BME = BME
        self.PSICO = PSICO
        self.SKLEARN = SKLEARN

        # Initialize some variables
        # Define different types of nucleus of interest
        self.proton_list = ["H1","H3","H1'", "H2'", "H3'", "H4'", "H5'",  "H5''", "H2", "H5", "H6", "H8", "H" "HN", "HA"]
        self.carbon_list = ["C1'", "C2'", "C3'", "C4'", "C5'", "C2", "C5", "C6", "C8", "CA", "C", "CB"]
        self.nitrogen_list = ["N1", "N3", "N"]

        self.total_list = self.proton_list + self.carbon_list + self.nitrogen_list

        self.ndisplayed = 10
        if self.ndisplayed > cmd.count_states(self.pymol_sel.get()):
            self.ndisplayed = cmd.count_states(self.pymol_sel.get())

        self.min_size.set(2)
        if self.min_size.get() > cmd.count_states(self.pymol_sel.get()):
            self.min_size.set(cmd.count_states(self.pymol_sel.get()))

        self.sel_obj_list = []
        # there may be more than one seletion or object defined by self.pymol_sel
        # treat each selection and object separately
        # System environment check
        if 'LARMORD_BIN' not in os.environ and 'PYMOL_GIT_MOD' in os.environ:
            if sys.platform.startswith('linux') and platform.machine() == 'x86_32':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'],"Larmord","i86Linux2","larmord")
                os.environ['LARMORD_BIN'] = initialdir_stride
            elif sys.platform.startswith('linux') and platform.machine() == 'x86_64':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'],"Larmord","ia64Linux2","larmord")
                os.environ['LARMORD_BIN'] = initialdir_stride
            else:
                pass
        # Set default file path
        if 'LARMORD_BIN' in os.environ:
            if VERBOSE: print('Found LARMORD_BIN in environmental variables', os.environ['LARMORD_BIN'])
            self.larmord_bin.set(os.environ['LARMORD_BIN'])
            #self.cs.set("/Users/afrankz/Documents/GitHub/PyShifts/test/protein_cs.dat")#for test use only
            #self.cs2.set("/Users/afrankz/Documents/GitHub/RNATransientStates/FluorideRibo/data/shifts.txt") #for test use only
            #self.pymol_sel.set("my_protein")#for test use only
            self.larmord_para.set(os.environ['LARMORD_BIN']+"/../data/larmorD_alphas_betas_rna.dat")
            self.larmord_ref.set(os.environ['LARMORD_BIN']+"/../data/larmorD_reference_shifts_rna.dat")

        else:
            if VERBOSE: print('LARMORD_BIN not found in environmental variables.')
            self.larmord_bin.set('')

        if 'PYSHIFTS_PATH' in os.environ:
            self.larmord_acc.set(os.environ['PYSHIFTS_PATH']+"/test/testAccuracy.dat")
        else:
            print("ERROR: Pyshifts Not Set Up! Please run $. setup.sh first.")
            sys.exit(1)


        if 'LARMORCA_BIN' not in os.environ and 'PYMOL_GIT_MOD' in os.environ:
            if sys.platform.startswith('linux') and platform.machine() == 'x86_32':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'],"Larmorca","i86Linux2","larmorca")
                os.environ['LARMORCA_BIN'] = initialdir_stride
            elif sys.platform.startswith('linux') and platform.machine() == 'x86_64':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'],"Larmorca","ia64Linux2","larmorca")
                os.environ['LARMORCA_BIN'] = initialdir_stride
            else:
                pass
        # Set default file path
        if 'LARMORCA_BIN' in os.environ:
            if VERBOSE: print('Found LARMORCA_BIN in environmental variables', os.environ['LARMORCA_BIN'])
            self.larmorca_bin.set(os.environ['LARMORCA_BIN'])
        else:
            if VERBOSE: print('LARMORCA_BIN not found in environmental variables.')
            self.larmorca_bin.set('')

        # tooltips
        self.balloon = Pmw.Balloon(self.parent)

        self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Exit',),
                                 title = 'PyShifts Plugin for PyMOL',
                                 command = self.execute)
        self.dialog.component('buttonbox').button(0).pack(fill='both',expand = 1, padx=10)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = tkinter.Label(self.dialog.interior(),text = 'PyShifts Plugin for PyMOL\nby  Jingru Xie, Kexin Zhang, and Aaron T. Frank, 2016\n',background = 'black', foreground = 'yellow')
        w.pack(expand = 1, fill = 'both', padx = 8, pady = 5)

        # add progress meter
        self.m = Meter(self.dialog.interior(), relief='ridge', bd=5)
        self.m.pack(expand = 1, padx = 10, pady = 5, fill='x')

        # make a few tabs within the dialog
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        pady=10
        self.notebook.pack(fill = 'both', expand=1, padx=10, pady=pady)

        ######################
        ## Tab : Options Tab
        ######################
        page = self.notebook.add('Compute or Load Shifts')
        self.notebook.tab('Compute or Load Shifts').focus_set()
        group_struc = tkinter.LabelFrame(page, text = 'Setup')
        group_struc.pack(fill='both', expand=True, padx=25, pady=25)
        self.analyzeButton = Pmw.ButtonBox(page,
                            orient = 'horizontal',
                            labelpos = 'w',
                            frame_borderwidth = 2,
                            frame_relief = 'groove',
                            )
        self.analyzeButton.add('Run', command = self.runAnalysis)
        #self.analyzeButton.pack(fill='both', expand=True, padx=25, pady=25)
        self.analyzeButton.button(0).grid(sticky=N, row=0)
        group_struc.grid(row = 0, column = 0)
        self.analyzeButton.grid(sticky=N, row=1)

        self.pymol_sel_ent = Pmw.EntryField(group_struc,
                                       label_text='PyMOL selection/object:',
                                       labelpos='wn',
                                       entry_textvariable=self.pymol_sel
                                       )
        self.balloon.bind(self.pymol_sel_ent, 'Enter the Pymol object for the structure(s) that should be used in this analysis')

        # Larmord_cs2 entry serves as the path to predicted chemical shifts file when in 'Analyze shifts' mode
        # Disabled when in 'Compute and Analyze shifts' mode
        self.cs2_ent = Pmw.EntryField(group_struc,
                                        label_text='Predicted Chemical Shift File:', labelpos='wn',
                                        entry_textvariable=self.cs2,
                                        entry_width=10)
        self.balloon.bind(self.cs2_ent, "Path to predicted chemical shifts data (only applicable in 'Other' mode) \n[Format:model_index, residue_number, residue_name, nucleus_name, CS_value_1, CS_values_2]                                          \nmodel_index - should match the state in the Pymol object to ensure accurate mapping of difference is subsequent analysis                                          \nresidue_number - same as above                                        \nresidue_name - same as above                                         \nnucleus_name - same as above                                    \nCS_values_1 - comparison chemical shifts that, for a given nucleus, may change with model_index                                         \nCS_values_2 - can be any value since it is ignored in this part ")
        self.cs2_but = tkinter.Button(group_struc, text = 'Browse...', command = self.getLarmordCS2)

        self.cs2_but.configure(state = "disabled")
        self.cs2_ent.component('entry').configure(state='disabled')

        # Mode selection: Larmord/ larmorca/ External
        self.mode_radio = Pmw.RadioSelect(group_struc,
                        buttontype = 'radiobutton',
                        orient = 'horizontal',
                        labelpos = 'we',
                        command = self.modeCallBack,
                        label_text = 'Mode',
                        hull_borderwidth = 2,
                        hull_relief = 'ridge',
                )

        self.balloon.bind(self.mode_radio, 'LARMORD mode: chemical shifts will be computed using LARMORD.\nLARMORCA mode: chemical shifts will be computed using larmorca.\nOther mode: chemical shifts to be compared will be read from the user supplied chemical shift file.')
        # Add some buttons to the radiobuttons RadioSelect.
        for text in ('LARMORD', 'LARMORCA', 'Other'):
            self.mode_radio.add(text)
        self.mode_radio.setvalue('LARMORD')

        # Arrange widgets using grid
        pady = 10
        padx = 5
        self.pymol_sel_ent.grid(sticky='we', row=0, column=0, columnspan=5,  pady=pady, padx=padx)
        self.mode_radio.grid(sticky='we', row=1, column=2, columnspan=3, pady=pady, padx=padx)
        self.cs2_ent.grid(sticky='we', row=2, column=0, columnspan=4, pady=pady, padx=padx)
        self.cs2_but.grid(sticky='we', row=2, column=4,  pady=pady, padx=padx)

        page.columnconfigure(0, weight=1)
        page.rowconfigure(0, weight = 1)
        page.rowconfigure(1, weight = 0)

        ######################
        ## Tab : Table Tab
        ######################
        page = self.notebook.add('Error Analysis')
        """
        Create an error Table
        """
        group_table = tkinter.LabelFrame(page)
        group_table.grid(sticky=W+E, row=0)
        self.tableButton = Pmw.ButtonBox(page,
                            orient = 'horizontal',
                            labelpos = 'w',
                            )
        self.tableButton.add('Compare Chemcial Shifts', command = self.runCompare)
        self.tableButton.add('Sort', command = self.runSort)
        self.tableButton.grid(sticky='we', row=1, columnspan=2)
        # disable all other buttons in the current tab before computing predicted chemical shifts
        self.tableButton.button(1).config(state = 'disabled')
        self.tableButton.button(0).config(state = 'disabled')

        self.topheader = tkinter.Label(group_table, text = '')

        fixedFont = Pmw.logicalfont('Fixed')

        # Begin Error Table definition
        self.error_table = Listbox(group_table,
            font = fixedFont,
            selectmode = BROWSE,
            width = 60,
            height = 10,
            bd = 6,
            relief='ridge',
        )
        self.error_table.pack(fill = BOTH, expand = 1)
        self.error_table.bind("<Double-Button-1>", self.sele_from_table)
        # Create the table header
        self.error_table.insert(1, 'Error Table'.center(50))
        self.table_header = 'state MAE R RMSE BME clusters'
        self.table_header = str.split(self.table_header)

        # Create row headers
        headerLine = ' '
        for column in range(len(self.table_header)):
            headerLine = headerLine + ('%-8s ' % (self.table_header[column],))
        self.error_table.insert('end', headerLine)

        """
        Create a chemical shift table with a scroll bar
        """
        group_CStable=Frame(group_table)
        cs_scrollbar = Scrollbar(group_CStable,
                   orient=VERTICAL,
                   bd = 1,
                   width = 16,
                   )
        self.CS_table = Listbox(group_CStable,
            font = fixedFont,
            selectmode = BROWSE,
            width = 60,
            height = 10,
            bd = 6,
            relief='ridge',
        )
        #cs_scrollbar.config(command=self.CS_table.yview)
        #cs_scrollbar.pack(side=RIGHT, fill=Y, )
        self.CS_table.pack(side=LEFT, fill = BOTH, expand = 1)
        self.CS_table.bind("<Double-Button-1>", self.seleRes)
        # Create the table header
        self.CS_table.insert(1, 'Chemical Shift Table'.center(55))
        self.CStable_header = 'resname resid nuclei expCS predCS weighted_error'
        self.CStable_header = str.split(self.CStable_header)

        # Create row headers
        headerLine = ' '
        for column in range(len(self.CStable_header)):
            headerLine = headerLine + ('%-6s ' % (self.CStable_header[column],))
        self.CS_table.insert('end', headerLine)

        # Chemical shift table sort - button box
        self.sort_CStable = Pmw.ButtonBox(group_table,
                        labelpos = 'w',
                        label_text = 'Sort CS table by:',
                        frame_borderwidth = 2,
                        frame_relief = 'groove',
                )
        self.balloon.bind(self.sort_CStable, 'Sort by residue id or scaled absolute error (in ascending order)')
        self.sort_CStable.add('resid', command = self.showCStable)
        self.sort_CStable.add('error', command = self.sort_CStable_error)
        # initialize buttons as disabled
        for button in range(self.sort_CStable.numbuttons()):
            self.sort_CStable.button(button).config(state = 'disabled')

        # Chemical shift table - save button
        self.save_CStable = Pmw.ButtonBox(group_table,
                       # orient = 'vertical',
                        labelpos = 'w',
                        label_text = 'Save to file:',
                        frame_borderwidth = 2,
                        frame_relief = 'groove',
                )
        self.balloon.bind(self.save_CStable, 'Save predicted chemical shift file')
        self.save_CStable.add('Error table', command = self.saveErrortable)
        self.save_CStable.add('CS table', command = self.saveCStable)
        self.save_CStable.add('Save single state', command = self.saveCStableOnestate)
       	 # button to save session
        self.save_CStable.add('Save session', command = self.saveSession)
        # initialize buttons as disabled
        for button in range(self.save_CStable.numbuttons()):
            self.save_CStable.button(button).config(state = 'disabled')

        # Reference chemical shift file
        self.cs_ent = Pmw.EntryField(group_table,
                                      label_text='Chemical Shift File:', labelpos='wn',
                                      entry_textvariable=self.cs,
                                      entry_width=10)
        self.balloon.bind(self.cs_ent, 'This file should contain the reference chemical shifts that will be \ncompared to the chemical shifts computed from the structure(s) using LARMORD or larmorca or reference chemical shifts supplied by the user.      \n[Format:residue_name, residue_number, nucleus_name, CS_value_1, CS_values_2]                                         \nresidue_number - residue number for a given chemical shift (should match the number in the load structure file)                                           \nresidue_name - residue name for a given chemical shift (should match the name in the load structure file)                                          \nnucleus_name - nucleus name for a given chemical shift (should match the name in the load structure file)                                          \nCS_values_1 - reference (measured) chemical shifts                                       \nCS_values_2 - can be any value since it is ignored in this part')
        self.cs_but = tkinter.Button(group_table, text = 'Browse...', command = self.getLarmordCS)


        """
        Radio button for nuclei selection --
        sort error table based on different nuclei types: 'carbon', 'proton', 'nitrogen' or 'all'
        """
        self.sortby_nucleus = Pmw.RadioSelect(group_table,
                        buttontype = 'radiobutton',
                        orient = 'horizontal',
                        labelpos = 'w',
                        command = self.sortByNucleusCallBack,
                        label_text = 'Comparison Nuclei:',
                        hull_borderwidth = 2,
                        hull_relief = 'ridge',
                )
        self.balloon.bind(self.sortby_nucleus, 'Sort using either proton, carbon, nitrogen or all chemical shifts')
        for text in ['proton', 'carbon', 'nitrogen', 'all']:
            self.sortby_nucleus.add(text)
        self.sortby_nucleus.setvalue('all')
        self.larmord_error_sel = 'all'

        """
        Radio button for metric selection --
        sort based on different attributes of selected nuclei error, e.g. MAE of carbon, Pearson correlation coefficients of proton, etc.
        """
        self.sortby_metric = Pmw.RadioSelect(group_table,
                        buttontype = 'radiobutton',
                        orient = 'horizontal',
                        labelpos = 'w',
                        command = self.sortByMetricCallBack,
                        label_text = 'Sorting Metric:',
                        hull_borderwidth = 2,
                        hull_relief = 'ridge',
                )
        self.balloon.bind(self.sortby_metric, 'Error metric to use when sorting states (state = reset)')
        for text in ['MAE', 'pearson', 'RMSE', 'state', 'weight', 'clusters']:
            self.sortby_metric.add(text)
        self.sortby_metric.setvalue('MAE')
        self.larmord_error_metric = 'MAE'

        # Arrange widgets using grid
        padx=1
        pady=1
        #self.topheader.grid(row=0)
        #self.error_table.grid(sticky=W, row=1, column=0, rowspan = 15, columnspan = 30, padx=padx, pady=pady)
        #group_CStable.grid(sticky=W, row=1, column=40, rowspan = 15, columnspan = 50, padx=padx, pady=pady)
        #self.save_CStable.grid(sticky=W, row=16, column=0, rowspan = 1, columnspan = 60, padx=padx, pady=pady)
        #self.sort_CStable.grid(sticky=W, row=16, column=60, rowspan = 1, columnspan = 30, padx=padx, pady=pady)
        #self.cs_ent.grid(sticky='we', row=17, column=0, columnspan=4, pady=pady, padx=padx)
        #self.cs_but.grid(sticky='we', row=17, column=4,  pady=pady, padx=padx)
        #self.sortby_nucleus.grid(sticky='we', row=18, column=0, columnspan = 50, pady=pady)
        #self.sortby_metric.grid(sticky='we', row=18, column=50, columnspan = 30, pady=pady)


        self.error_table.grid(sticky=W, row=1, column=0, rowspan = 5, columnspan = 30, padx=padx, pady=pady)
        group_CStable.grid(sticky=W, row=6, column=0, rowspan = 5, columnspan = 50, padx=padx, pady=pady)
        self.save_CStable.grid(sticky=W, row=11, column=0, rowspan = 1, columnspan = 60, padx=padx, pady=pady)
        self.sort_CStable.grid(sticky=W, row=12, column=0, rowspan = 1, columnspan = 30, padx=padx, pady=pady)
        self.cs_ent.grid(sticky='we', row=13, column=0, columnspan=4, pady=pady, padx=padx)
        self.cs_but.grid(sticky='we', row=13, column=4,  pady=pady, padx=padx)
        self.sortby_nucleus.grid(sticky='we', row=14, column=0, columnspan = 50, pady=pady)
        self.sortby_metric.grid(sticky='we', row=15, column=0, columnspan = 30, pady=pady)


        group_table.columnconfigure(0, weight=1)
        group_table.columnconfigure(50, weight=1)
        group_table.columnconfigure(90, weight=1)

        page.columnconfigure(0, weight=1)

        ###########################
        ## Tab : Advanced Options Tab
        ###########################
        page = self.notebook.add('Advanced Options')
        """
        Larmord file path
        """
        Larmord_file = False
        if Larmord_file:
            group_larmord = tkinter.LabelFrame(page, text = 'LARMORD')
            group_larmord.pack(fill='both', expand=True, padx=5, pady=5)

            enwidth=20
            self.larmord_bin_ent = Pmw.EntryField(group_larmord,
                                          label_text='LARMORD binary:', labelpos='wn',
                                          entry_textvariable=self.larmord_bin,
                                          entry_width=enwidth)
            self.balloon.bind(self.larmord_bin_ent, 'Path to LARMORD binary (only needed if Compute and Analyze mode )')

            self.larmord_bin_but = tkinter.Button(group_larmord, text = 'Browse...', command = self.getLarmordBin)

            self.larmord_para_ent = Pmw.EntryField(group_larmord,
                                          label_text='LARMORD Parameter:', labelpos='wn',
                                          entry_textvariable=self.larmord_para,
                                          entry_width=enwidth)
            self.balloon.bind(self.larmord_para_ent, 'Path to LARMORD parameter file (only needed if Compute and Analyze mode )')

            self.larmord_para_but = tkinter.Button(group_larmord, text = 'Browse...', command = self.getLarmordPara)

            self.larmord_ref_ent = Pmw.EntryField(group_larmord,
                                          label_text='LARMORD Reference Shifts:', labelpos='wn',
                                          entry_textvariable=self.larmord_ref,
                                          entry_width=enwidth)
            self.balloon.bind(self.larmord_ref_ent, 'Path to LARMORD reference chemical shifts file (only needed if Compute and Analyze mode )')

            self.larmord_ref_but = tkinter.Button(group_larmord, text = 'Browse...', command = self.getLarmordRef)

            # arrange widgets using grid
            pady=5
            self.larmord_bin_ent.grid(sticky='we', row=0, column=0, columnspan=3, pady=pady)
            self.larmord_bin_but.grid(sticky='we', row=0, column=3,  pady=pady)
            self.larmord_para_ent.grid(sticky='we', row=1, column=0, columnspan=3,  pady=pady)
            self.larmord_para_but.grid(sticky='we', row=1, column=3,  pady=pady)
            self.larmord_ref_ent.grid(sticky='we', row=2, column=0, columnspan=3,  pady=pady)
            self.larmord_ref_but.grid(sticky='we', row=2, column=3,  pady=pady)

        group_advanc = tkinter.LabelFrame(page, text = 'Advanced Options')
        group_advanc.pack(fill='both', expand=True, padx=5, pady=5)

        """
        Scale options: scale or not/ offset for proton, carbon and nitrogen/ customize larmord or larmorca accuracy file
        """
        group_scale = group_advanc

        # Scale selection: scale/ not scale
        self.wterror_radio = Pmw.RadioSelect(group_scale,
                        buttontype = 'radiobutton',
                        orient = 'horizontal',
                        labelpos = 'w',
                        command = self.scaleCallBack,
                        label_text = 'Scale \nDifferences',
                        hull_borderwidth = 2,
                        hull_relief = 'ridge',
                )
        self.balloon.bind(self.wterror_radio, 'Should chemical differences should be scaled by the expected accuracy of your predictor')
        for text in ('Yes', 'No'):
            self.wterror_radio.add(text)
        self.wterror_radio.setvalue('Yes')
        self.weighted_errors = True

        # 1H and 13C and 15N offset
        self.larmord_proton_offset_ent = Pmw.EntryField(group_scale, value = 0.0,
                                       label_text='1H offset (ppm):',
                                       labelpos='wn',
                                       entry_textvariable=self.larmord_proton_offset,
                                       )
        self.balloon.bind(self.larmord_proton_offset_ent, 'offset to add to the reference 1H chemical shifts -- useful if reference shifts are known to be systematically mis-referenced')

        self.larmord_carbon_offset_ent = Pmw.EntryField(group_scale, value = 0.0,
                                       label_text='13C offset (ppm):',
                                       labelpos='wn',
                                       entry_textvariable=self.larmord_carbon_offset,
                                       )
        self.balloon.bind(self.larmord_carbon_offset_ent, 'offset to add to the reference 13C chemical shifts -- useful if reference shifts are known to be systematically mis-referenced')

        self.larmord_nitrogen_offset_ent = Pmw.EntryField(group_scale, value = 0.0,
                                       label_text='15N offset (ppm):',
                                       labelpos='wn',
                                       entry_textvariable=self.larmord_nitrogen_offset,
                                       )
        self.balloon.bind(self.larmord_nitrogen_offset_ent, 'offset to add to the reference 15N chemical shifts -- useful if reference shifts are known to be systematically mis-referenced')

        # outlier threshold
        self.larmord_outlier_threshold_ent = Pmw.EntryField(group_scale, value=9999.0,
                                                        label_text='Outlier Threshold:',
                                                        labelpos='wn',
                                                        entry_textvariable=self.larmord_outlier_threshold,
                                                        )
        self.balloon.bind(self.larmord_outlier_threshold_ent,
                          'outlier_threshold -- predictions that exhibit errors that outlier_threshold x the expected error will be ignored was determining the overall error of a given structural model')

        # Specify accuracy file
        # load MAE from file if given
        # MAE: estimated mean absolute error between measured and larmord predicted chemical shifts for a certain nucleus type
        self.larmord_acc_ent = Pmw.EntryField(group_scale,
                                      label_text='Accuracy File:', labelpos='wn',
                                      entry_textvariable=self.larmord_acc,
                                      entry_width=10)
        self.balloon.bind(self.larmord_acc_ent, 'Path to customized LARMORD MAE data file\nMAE: estimated mean absolute error between measured and larmord predicted chemical shifts for a certain nucleus type')
        self.larmord_acc_but = tkinter.Button(group_scale, text = 'Browse...', command = self.getLarmordAcc, width=4)

        """
        Rendering options
        -- specify colors, sphere size of residuals and No. of models to display in Pymol
        """
        group_error = group_advanc
        self.error_color_ent = Pmw.EntryField(group_error,
                                       label_text="Colors:",
                                       labelpos='wn', value = "blue_white_red",
                                       entry_textvariable=self.larmord_error_color
                                       )
        self.balloon.bind(self.error_color_ent, "see 'help spectrum' for list of color spectra")

        self.error_scale_ent = Pmw.EntryField(group_error,
                                       label_text='Sphere size:',
                                       labelpos='wn', value = 0.4,
                                       entry_textvariable=self.larmord_error_scale
                                       )
        self.balloon.bind(self.error_scale_ent, "scale for spheres used to render the chemical shift differences")

        self.error_ndisplay_ent = Pmw.EntryField(group_error,
                        labelpos = 'w',
                        label_text='No. models to display:',
                        value = '%s'%(self.ndisplayed),
                        entry_width = 3,
                        entry_textvariable=self.larmord_ndisplayed
                        )
        self.balloon.bind(self.error_ndisplay_ent, "Number of models with lowest error (or highest correlation coefficiencts) to display")

        self.min_size_ent = Pmw.EntryField(group_error,
                        labelpos = 'w',
                        label_text='Number of clusters:',
                        value = 2,
                        entry_width = 3,
                        entry_textvariable=self.min_size
                        )
        self.balloon.bind(self.min_size_ent, "Number of clusters K-means should return")


        # Arrange widgets using grid
        padx = 5
        pady = 5
        self.wterror_radio.grid(sticky='we', row=0, column=0, pady=pady, padx=padx)
        self.larmord_proton_offset_ent.grid(sticky='we', row=1, column=0, columnspan=2, pady=pady, padx=padx)
        self.larmord_carbon_offset_ent.grid(sticky='we', row=2, column=0, columnspan=2, pady=pady, padx=padx)
        self.larmord_nitrogen_offset_ent.grid(sticky='we', row=3, column=0, columnspan=2, pady=pady, padx=padx)
        self.larmord_outlier_threshold_ent.grid(sticky='we', row=4, column=0, columnspan=2, pady=pady, padx=padx)
        self.larmord_acc_ent.grid(sticky='we', row=5, column=0, columnspan=1,  pady=pady, padx=padx)
        self.larmord_acc_but.grid(sticky='we', row=5, column=1, columnspan=1, pady=pady, padx=padx)
        self.error_color_ent.grid(sticky='we', row=6, column=0, columnspan=2, pady=pady, padx=padx)
        self.error_scale_ent.grid(sticky='we', row=7, column=0, columnspan=2, pady=pady, padx=padx)
        self.error_ndisplay_ent.grid(sticky='we', row=8, column=0, columnspan=2, pady=pady, padx=padx)
        self.min_size_ent.grid(sticky='we', row=9, column=0, columnspan=2, pady=pady, padx=padx)

        group_advanc.grid(row=0, column=0)
        page.columnconfigure(0, weight=1)
        page.rowconfigure(0, weight = 1)


        ######################
        ## Tab : About Tab
        ######################
        page = self.notebook.add('About')
        group_about = tkinter.LabelFrame(page, text = 'PyShifts Plugin for PyMOL')
        group_about.grid(sticky='we', row=0,column=0,padx=5,pady=3)
        about_plugin = """ Tool For Comparing and Visualizing Chemical Shift Differences.
                    Jingru Xie <jingrux .at. umich.edu>
                    Kexin Zhang <kexin .at. umich.edu>
                    Aaron T. Frank  <afrankz .at. umich.edu>
                    Please cite this plugin if you use it in a publication.
                        """
        label_about = tkinter.Label(group_about,text=about_plugin)
        label_about.grid(sticky='we', row=0, column=0,  pady=pady)
        self.notebook.setnaturalsize()
        return

    ## Functions related to file and parameter initialization
    def disableAll(self):
        """
        Disable all buttons in the dialog(except for 'Exit')
        Will be called when the user click on 'Predict or Load Chemical Shifts' or 'Compare shifts'
        Purpose: avoid unnecassary repeted running caused by double click
        """
        for button in range(self.tableButton.numbuttons()):
            self.tableButton.button(button).config(state = 'disabled')
        self.analyzeButton.button(0).config(state='disabled')

    ### A set of functions for getting Larmord parameters
    def getLarmordBin(self):
        larmord_bin_fname = tkFileDialog.askopenfilename(title='Larmord Binary', initialdir='',filetypes=[('all','*')], parent=self.parent)
        if larmord_bin_fname: # if nonempty
            self.larmord_bin.set(larmord_bin_fname)
        return

    def getLarmordPara(self):
        larmord_para_fname = tkFileDialog.askopenfilename(title='Larmord Parameter', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_para_fname: # if nonempty
            self.larmord_para.set(larmord_para_fname)
        return

    def getLarmordRef(self):
        larmord_ref_fname = tkFileDialog.askopenfilename(title='Larmord Reference', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_ref_fname: # if nonempty
            self.larmord_ref.set(larmord_ref_fname)
        return

    def getLarmordAcc(self):
        larmord_acc_fname = tkFileDialog.askopenfilename(title='Accuracy File', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_acc_fname: # if nonempty
            self.larmord_acc.set(larmord_acc_fname)
        return

    def getLarmordCS(self):
        larmord_cs_fname = tkFileDialog.askopenfilename(title='Chemical Shift File', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_cs_fname: # if nonempty
            self.cs.set(larmord_cs_fname)
        return

    def getLarmordCS2(self):
        larmord_cs2_fname = tkFileDialog.askopenfilename(title='Predicted Chemical Shift File', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_cs2_fname: # if nonempty
            self.cs2.set(larmord_cs2_fname)
        return

    def disableNuclei(self):
        # Disable the selection of nuclei in 'table' tab
        # Will be called when metric selection is in 'state'
        self.sortby_nucleus.component('label').configure(state='disabled')
        for child in range(self.sortby_nucleus.numbuttons()):
            self.sortby_nucleus.button(child).configure(state='disabled')
        return True

    def enableNuclei(self):
        # Enable the selection of nuclei in 'table' tab
        # Will be called when metric selection is not in 'state'
        self.sortby_nucleus.component('label').configure(state='normal')
        for child in range(self.sortby_nucleus.numbuttons()):
            self.sortby_nucleus.button(child).configure(state='normal')
        return True

    def sortByMetricCallBack(self, tag):
        if tag in ['pearson']:
            self.enableNuclei()
            self.larmord_error_metric = 'pearson'
        if tag in ['kendall']:
            self.enableNuclei()
            self.larmord_error_metric = 'kendall'
        if tag in ['spearman']:
            self.enableNuclei()
            self.larmord_error_metric = 'spearman'
        if tag in ['MAE']:
            self.enableNuclei()
            self.larmord_error_metric = 'MAE'
        if tag in ['RMSE']:
            self.enableNuclei()
            self.larmord_error_metric = 'RMSE'
        if tag in ['PYM']:
            self.enableNuclei()
            self.larmord_error_metric = 'PYM'
        if tag in ['state']:
            self.larmord_error_metric = 'state'
            self.disableNuclei()
        if tag in ['weight']:
            self.larmord_error_metric = 'weight'
            self.disableNuclei()
        if tag in ['clusters']:
            self.larmord_error_metric = 'clusters'
            self.disableNuclei()

    def sortByNucleusCallBack(self, tag):
        self.larmord_error_sel = 'all'
        #self.tableButton.button(1).config(state = 'disabled')
        if tag in ['proton']:
            self.larmord_error_sel = 'proton'
        if tag in ['carbon']:
            self.larmord_error_sel = 'carbon'
        if tag in ['nitrogen']:
            self.larmord_error_sel = 'nitrogen'

    def modeCallBack(self, tag):
        # This is called whenever the user clicks on a button
        # during the selection of modes.
        Larmord_file = False
        if tag in ['Other']:
            self.get_shifts_from_larmord = False
            self.get_shifts_from_larmorca = False
            self.get_shifts_from_file = True
            self.cs2_but.configure(state = "normal")
            self.cs2_ent.component('entry').configure(state='normal')
            self.wterror_radio.component('label').configure(state='disabled')
            for child in range(self.wterror_radio.numbuttons()):
                 self.wterror_radio.button(child).configure(state='disabled')
        elif tag in ['LARMORCA']:
            self.get_shifts_from_larmord = False
            self.get_shifts_from_larmorca = True
            self.get_shifts_from_file = False
            self.wterror_radio.component('label').configure(state='normal')
            self.larmord_acc_but.configure(state = "normal")
            self.larmord_acc_ent.component('entry').configure(state='normal')
            self.cs2_but.configure(state = "disabled")
            self.cs2_ent.component('entry').configure(state='disabled')
            for child in range(self.wterror_radio.numbuttons()):
                 self.wterror_radio.button(child).configure(state='normal')
        elif tag in ['LARMORD']:
            self.get_shifts_from_larmord = True
            self.get_shifts_from_larmorca = False
            self.get_shifts_from_file = False
            if Larmord_file:
                # activate LARMORD relevant buttons and entry fields
                self.larmord_bin_but.configure(state = "normal")
                self.larmord_para_but.configure(state = "normal")
                self.larmord_ref_but.configure(state = "normal")
                self.larmord_bin_ent.component('entry').configure(state='normal')
                self.larmord_para_ent.component('entry').configure(state='normal')
                self.larmord_ref_ent.component('entry').configure(state='normal')
            self.wterror_radio.component('label').configure(state='normal')
            self.larmord_acc_but.configure(state = "normal")
            self.larmord_acc_ent.component('entry').configure(state='normal')
            self.cs2_but.configure(state = "disabled")
            self.cs2_ent.component('entry').configure(state='disabled')
            for child in range(self.wterror_radio.numbuttons()):
                 self.wterror_radio.button(child).configure(state='normal')
        else:
            self.get_shifts_from_larmord = True
            self.get_shifts_from_larmorca = True
            self.get_shifts_from_file = False
            self.wterror_radio.component('label').configure(state='normal')
            self.larmord_acc_but.configure(state = "normal")
            self.larmord_acc_ent.component('entry').configure(state='normal')
            self.cs2_but.configure(state = "disabled")
            self.cs2_ent.component('entry').configure(state='disabled')
            for child in range(self.wterror_radio.numbuttons()):
                 self.wterror_radio.button(child).configure(state='normal')

    def scaleCallBack(self, tag):
        # This is called whenever the user clicks on a button in the single selection of 'scale error' or 'not scale error' and pass the selction to
        # a global var self.weighted_errors.
        # If the user chooses to scale the error, self.weighted_errors = True; and vice versa.
        if tag in ['Yes']:
            self.weighted_errors = True
        else:
            self.weighted_errors = False

    def reset_predictedCS(self):
        # reset predicted chemical shift
        self.predictedCS = []

    def get_shifts_from_traj(self):
        """
        Run larmord in trajectory mode
        """

    def reset_MAE(self):
        self.mae = {}

    def load_MAE(self):
        # Load larmord accuracy(expected MAE for each residuals in larmord) from file.
        # If MAE file is not provided by the user, then load MAE from default list also given in this function.
        self.reset_MAE()
        if len(self.larmord_acc.get())>0:
            if self.check_file(self.larmord_acc.get()):
                mae_type = {'names': ('resname', 'nucleus','MAE'),'formats': ('S5','S5','float')}
                #data = np.loadtxt(self.larmord_acc.get(), dtype=mae_type)
                data = pd.read_csv(self.larmord_acc.get(), sep = '\s+', names = ['resname', 'nucleus','MAE'])
                self.larmord_acc_ext = data
                larmord_resname = data['resname'].values
                larmord_nucleus = data['nucleus'].values
                larmord_mae = data['MAE'].values
                print(data)
                for res in range(len(larmord_resname)):
                    keyMAE = str(larmord_resname[res]+":"+larmord_nucleus[res])
                    self.mae[keyMAE] = larmord_mae[res]
            else:
                #self.print_file_error(self.larmord_acc.get())
                return False

        # default MAE list
        else:
            print("loading default MAE...")
            # MAE for Larmord
            if (self.get_shifts_from_larmord) and not (self.get_shifts_from_larmord):
                self.mae["ADE:C1'"] = 0.700
                self.mae["ADE:C2"] = 0.764
                self.mae["ADE:C2'"] = 0.471
                self.mae["ADE:C3'"] = 1.017
                self.mae["ADE:C4'"] = 0.709
                self.mae["ADE:C5'"] = 0.914
                self.mae["ADE:C8"] = 0.684
                self.mae["ADE:H1'"] = 0.114
                self.mae["ADE:H2"] = 0.205
                self.mae["ADE:H2'"] = 0.107
                self.mae["ADE:H3'"] = 0.113
                self.mae["ADE:H4'"] = 0.074
                self.mae["ADE:H5'"] = 0.314
                self.mae["ADE:H5''"] = 0.132
                self.mae["ADE:H8"] = 0.161
                self.mae["GUA:C1'"] = 0.681
                self.mae["GUA:C2'"] = 0.546
                self.mae["GUA:C3'"] = 0.861
                self.mae["GUA:C4'"] = 0.817
                self.mae["GUA:C5'"] = 0.893
                self.mae["GUA:C8"] = 0.715
                self.mae["GUA:H1'"] = 0.168
                self.mae["GUA:H2'"] = 0.142
                self.mae["GUA:H3'"] = 0.124
                self.mae["GUA:H4'"] = 0.075
                self.mae["GUA:H5'"] = 0.271
                self.mae["GUA:H5''"] = 0.123
                self.mae["GUA:H8"] = 0.202
                self.mae["URA:C1'"] = 0.681
                self.mae["URA:C2'"] = 0.598
                self.mae["URA:C3'"] = 0.981
                self.mae["URA:C4'"] = 1.106
                self.mae["URA:C5"] = 0.562
                self.mae["URA:C5'"] = 0.735
                self.mae["URA:C6"] = 0.705
                self.mae["URA:H1'"] = 0.105
                self.mae["URA:H2'"] = 0.129
                self.mae["URA:H3'"] = 0.088
                self.mae["URA:H4'"] = 0.071
                self.mae["URA:H5"] = 0.120
                self.mae["URA:H5'"] = 0.303
                self.mae["URA:H5''"] = 0.121
                self.mae["URA:H6"] = 0.099
                self.mae["CYT:C1'"] = 0.642
                self.mae["CYT:C2'"] = 0.980
                self.mae["CYT:C3'"] = 1.147
                self.mae["CYT:C4'"] = 0.617
                self.mae["CYT:C5"] = 1.945
                self.mae["CYT:C5'"] = 0.938
                self.mae["CYT:C6"] = 0.584
                self.mae["CYT:H1'"] = 0.111
                self.mae["CYT:H2'"] = 0.101
                self.mae["CYT:H3'"] = 0.113
                self.mae["CYT:H4'"] = 0.094
                self.mae["CYT:H5"] = 0.114
                self.mae["CYT:H5'"] = 0.312
                self.mae["CYT:H5''"] = 0.194
                self.mae["CYT:H6"] = 0.117
                self.mae["URA:N3"] = 2.609
                self.mae["GUA:N1"] = 1.259
            # MAE for larmorca
            if not (self.get_shifts_from_larmord) and (self.get_shifts_from_larmorca):
                self.mae["ALA:HA"] = 0.39
                self.mae["ARG:HA"] = 0.39
                self.mae["ASN:HA"] = 0.39
                self.mae["ASP:HA"] = 0.39
                self.mae["CYS:HA"] = 0.39
                self.mae["GLN:HA"] = 0.39
                self.mae["GLU:HA"] = 0.39
                self.mae["GLY:HA"] = 0.39
                self.mae["HIS:HA"] = 0.39
                self.mae["ILE:HA"] = 0.39
                self.mae["LEU:HA"] = 0.39
                self.mae["LYS:HA"] = 0.39
                self.mae["MET:HA"] = 0.39
                self.mae["PHE:HA"] = 0.39
                self.mae["PRO:HA"] = 0.39
                self.mae["SER:HA"] = 0.39
                self.mae["THR:HA"] = 0.39
                self.mae["TRP:HA"] = 0.39
                self.mae["TYR:HA"] = 0.39
                self.mae["VAL:HA"] = 0.39
                self.mae["ALA:HN"] = 0.25
                self.mae["ARG:HN"] = 0.25
                self.mae["ASN:HN"] = 0.25
                self.mae["ASP:HN"] = 0.25
                self.mae["CYS:HN"] = 0.25
                self.mae["GLN:HN"] = 0.25
                self.mae["GLU:HN"] = 0.25
                self.mae["GLY:HN"] = 0.25
                self.mae["HIS:HN"] = 0.25
                self.mae["ILE:HN"] = 0.25
                self.mae["LEU:HN"] = 0.25
                self.mae["LYS:HN"] = 0.25
                self.mae["MET:HN"] = 0.25
                self.mae["PHE:HN"] = 0.25
                self.mae["PRO:HN"] = 0.25
                self.mae["SER:HN"] = 0.25
                self.mae["THR:HN"] = 0.25
                self.mae["TRP:HN"] = 0.25
                self.mae["TYR:HN"] = 0.25
                self.mae["VAL:HN"] = 0.25
                self.mae["ALA:CA"] = 0.9
                self.mae["ARG:CA"] = 0.9
                self.mae["ASN:CA"] = 0.9
                self.mae["ASP:CA"] = 0.9
                self.mae["CYS:CA"] = 0.9
                self.mae["GLN:CA"] = 0.9
                self.mae["GLU:CA"] = 0.9
                self.mae["GLY:CA"] = 0.9
                self.mae["HIS:CA"] = 0.9
                self.mae["ILE:CA"] = 0.9
                self.mae["LEU:CA"] = 0.9
                self.mae["LYS:CA"] = 0.9
                self.mae["MET:CA"] = 0.9
                self.mae["PHE:CA"] = 0.9
                self.mae["PRO:CA"] = 0.9
                self.mae["SER:CA"] = 0.9
                self.mae["THR:CA"] = 0.9
                self.mae["TRP:CA"] = 0.9
                self.mae["TYR:CA"] = 0.9
                self.mae["VAL:CA"] = 0.9
                self.mae["ALA:C"] = 1.03
                self.mae["ARG:C"] = 1.03
                self.mae["ASN:C"] = 1.03
                self.mae["ASP:C"] = 1.03
                self.mae["CYS:C"] = 1.03
                self.mae["GLN:C"] = 1.03
                self.mae["GLU:C"] = 1.03
                self.mae["GLY:C"] = 1.03
                self.mae["HIS:C"] = 1.03
                self.mae["ILE:C"] = 1.03
                self.mae["LEU:C"] = 1.03
                self.mae["LYS:C"] = 1.03
                self.mae["MET:C"] = 1.03
                self.mae["PHE:C"] = 1.03
                self.mae["PRO:C"] = 1.03
                self.mae["SER:C"] = 1.03
                self.mae["THR:C"] = 1.03
                self.mae["TRP:C"] = 1.03
                self.mae["TYR:C"] = 1.03
                self.mae["VAL:C"] = 1.03
                self.mae["ALA:CB"] = 1.55
                self.mae["ARG:CB"] = 1.55
                self.mae["ASN:CB"] = 1.55
                self.mae["ASP:CB"] = 1.55
                self.mae["CYS:CB"] = 1.55
                self.mae["GLN:CB"] = 1.55
                self.mae["GLU:CB"] = 1.55
                self.mae["GLY:CB"] = 1.55
                self.mae["HIS:CB"] = 1.55
                self.mae["ILE:CB"] = 1.55
                self.mae["LEU:CB"] = 1.55
                self.mae["LYS:CB"] = 1.55
                self.mae["MET:CB"] = 1.55
                self.mae["PHE:CB"] = 1.55
                self.mae["PRO:CB"] = 1.55
                self.mae["SER:CB"] = 1.55
                self.mae["THR:CB"] = 1.55
                self.mae["TRP:CB"] = 1.55
                self.mae["TYR:CB"] = 1.55
                self.mae["VAL:CB"] = 1.55
                self.mae["ALA:N"] = 2.45
                self.mae["ARG:N"] = 2.45
                self.mae["ASN:N"] = 2.45
                self.mae["ASP:N"] = 2.45
                self.mae["CYS:N"] = 2.45
                self.mae["GLN:N"] = 2.45
                self.mae["GLU:N"] = 2.45
                self.mae["GLY:N"] = 2.45
                self.mae["HIS:N"] = 2.45
                self.mae["ILE:N"] = 2.45
                self.mae["LEU:N"] = 2.45
                self.mae["LYS:N"] = 2.45
                self.mae["MET:N"] = 2.45
                self.mae["PHE:N"] = 2.45
                self.mae["PRO:N"] = 2.45
                self.mae["SER:N"] = 2.45
                self.mae["THR:N"] = 2.45
                self.mae["TRP:N"] = 2.45
                self.mae["TYR:N"] = 2.45
                self.mae["VAL:N"] = 2.45
            # MAE of the combined and other mode
            # Average of the two above
            else:
                self.mae["ADE:C1'"]=0.692
                self.mae["ADE:C2"]=0.624
                self.mae["ADE:C2'"]=0.467
                self.mae["ADE:C3'"]=0.976
                self.mae["ADE:C4'"]=0.672
                self.mae["ADE:C5'"]=0.866
                self.mae["ADE:C8"]=0.676
                self.mae["ADE:H1'"]=0.140
                self.mae["ADE:H2"]=0.195
                self.mae["ADE:H2'"]=0.123
                self.mae["ADE:H3'"]=0.110
                self.mae["ADE:H4'"]=0.073
                self.mae["ADE:H5'"]=0.240
                self.mae["ADE:H5''"]=0.113
                self.mae["ADE:H8"]=0.163
                self.mae["CYT:C1'"]=0.602
                self.mae["CYT:C2'"]=0.716
                self.mae["CYT:C3'"]=1.030
                self.mae["CYT:C4'"]=0.507
                self.mae["CYT:C5"]=1.199
                self.mae["CYT:C5'"]=0.814
                self.mae["CYT:C6"]=0.562
                self.mae["CYT:H1'"]=0.136
                self.mae["CYT:H2'"]=0.111
                self.mae["CYT:H3'"]=0.106
                self.mae["CYT:H4'"]=0.088
                self.mae["CYT:H5"]=0.121
                self.mae["CYT:H5'"]=0.223
                self.mae["CYT:H5''"]=0.150
                self.mae["CYT:H6"]=0.113
                self.mae["GUA:C1'"]=0.739
                self.mae["GUA:C2'"]=0.575
                self.mae["GUA:C3'"]=0.899
                self.mae["GUA:C4'"]=0.735
                self.mae["GUA:C5'"]=0.838
                self.mae["GUA:C8"]=0.801
                self.mae["GUA:H1'"]=0.178
                self.mae["GUA:H2'"]=0.131
                self.mae["GUA:H3'"]=0.122
                self.mae["GUA:H4'"]=0.076
                self.mae["GUA:H5'"]=0.198
                self.mae["GUA:H5''"]=0.114
                self.mae["GUA:H8"]=0.194
                self.mae["URA:C1'"]=0.676
                self.mae["URA:C2'"]=0.509
                self.mae["URA:C3'"]=0.964
                self.mae["URA:C4'"]=0.888
                self.mae["URA:C5"]=0.742
                self.mae["URA:C5'"]=0.813
                self.mae["URA:C6"]=0.768
                self.mae["URA:H1'"]=0.123
                self.mae["URA:H2'"]=0.135
                self.mae["URA:H3'"]=0.093
                self.mae["URA:H4'"]=0.078
                self.mae["URA:H5"]=0.141
                self.mae["URA:H5'"]=0.220
                self.mae["URA:H5''"]=0.107
                self.mae["URA:H6"]=0.108

    def reset_measuredCS(self):
        self.measuredCS = {}

    def load_measuredCS(self):
        # load measured Chemical shift data from file
        # If measurd CS file not given, an error box will pop up and return false

        self.reset_measuredCS()
        measuredCS_fn = self.cs.get()
        print('loading measured chemical shift from file: %s...'%measuredCS_fn)
        if self.check_file(measuredCS_fn):
            # read in from external file
            if measuredCS_fn.lower().endswith(('nmrstar', '.str')):
                print('Converting nmrstar chemical shift to measured chemical shift')
                from nmrstarConverter import nmrstar2measure
                tmp_measured_file =  "tmp_measured_shifts.dat"
                nmrstar_fn = "tmp_nmrstar.str"
                os.system("sed -n '/_Atom_chem_shift/,/stop_/p' %s | grep -v stop_ > %s" %(measuredCS_fn, nmrstar_fn))
                nmrstar2measure(nmrstar_fn, tmp_measured_file)
                expCS_data = pd.read_csv(tmp_measured_file, sep = '\s+', names = ['resname', 'resid', 'nucleus','expCS', 'junk'])
                os.remove(tmp_measured_file)
                os.remove(nmrstar_fn)
            else:
                expCS_data = pd.read_csv(measuredCS_fn, sep = '\s+', names = ['resname', 'resid', 'nucleus','expCS', 'junk'])
            self.expCS_data_ext = expCS_data
            measured_resname = expCS_data['resname'].values
            measured_resid = expCS_data['resid'].values
            measured_nucleus = expCS_data['nucleus'].values
            measured_csvalue = expCS_data['expCS'].values
            for res in range(len(measured_resid)):
                k2 = str(measured_resid[res])+":"+str(measured_resname[res])+":"+str(measured_nucleus[res])
                self.measuredCS[k2] = measured_csvalue[res]
        else:
            print(self.predictedCS[0].keys())
            self.reset_measuredCS()
            for key in self.predictedCS[0].keys():
                self.measuredCS[key] = self.predictedCS[0][key]

    ## Functions related to self.runAnalysis()
    def conv_resname_format(self):
        """ Convert the format of residual names to a consistent 3-capital-character form, e.g. GUA
            Dictionary of measured CS data already loaded as self.measuredCS
            GUA = G RG RG3 RG5 RGA; ADE = A RA RA3 RA5 RAD; CYT = C RC RC3 RC5 RCY; URA = U RU RU3 RU5
            Loop over the dictionary to convert each elements into desired form
        """
        for key in self.measuredCS.keys():
            resid, resname, nucleus = key.split(":")
            if resname.upper() in ['CYT','C','CYTOSINE','RC','RC3','RC5','RCY']:
                resname = 'CYT'
            elif resname.upper() in ['URA','U','URACIL','RU','RU3','RU5']:
                resname = 'URA'
            elif resname.upper() in ['GUA','G','GUANINE','RG','RG3','RG5','RGA']:
                resname = 'GUA'
            elif resname.upper() in ['ADE','A','ADENINE','RA','RA3','RA5','RAD']:
                resname = 'ADE'
            elif resname.upper() in ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]:
                resname = resname.upper()
            else:
                print('Error: residue name %s not found!' %resname)
            key_new = str(resid+":"+resname+":"+nucleus)
            if key_new != key:
                self.measuredCS[key_new] = self.measuredCS[key]
                del self.measuredCS[key]

    def parse_larmor_output(self, larmord_tmpout_fn):
        """ Parse Larmord output to self.predictedCS
            @param larmord_tmpout_fn: larmord output file
            @param type: string
            @output self.predictedCS: predicted chemical shift data provided by larmord output file
            @output type: a list of dictionaries. List index: state number; dictionary key: redidues (id + name + nucleus), dictionary value: predicted chemical shifts data
        """
        predCS = {}
        predCS_data = pd.read_csv(larmord_tmpout_fn, sep = '\s+', names = ['resid', 'resname', 'nucleus', 'junk', 'predCS'])
        larmord_resid = predCS_data['resid'].values
        larmord_resname = predCS_data['resname'].values
        larmord_nucleus = predCS_data['nucleus'].values
        larmord_predCS = predCS_data['predCS'].values
        # print larmord_predCS

        for res in range(len(larmord_resid)):
            #keypred = str(str(larmord_resid[res])+":"+larmord_resname[res]+":"+larmord_nucleus[res]).strip(sep=b'')
            keypred = str(str(larmord_resid[res])+":"+larmord_resname[res]+":"+larmord_nucleus[res])
            # generate keys to self.predictedCS
            predCS[keypred] = larmord_predCS[res]
        return predCS

    def parse_larmord_dataframe(self, data_ext, state):
        """ Parse Larmord output to self.predictedCS
            @param larmord_tmpout_fn: larmord output file
            @param type: string
            @output self.predictedCS: predicted chemical shift data provided by larmord output file
            @output type: a list of dictionaries. List index: state number; dictionary key: redidues (id + name + nucleus), dictionary value: predicted chemical shifts data
        """
        predCS = {}
        data_ext = data_ext[data_ext['state']==state]
        larmord_resid = data_ext['resid'].values
        larmord_resname = data_ext['resname'].values
        larmord_nucleus = data_ext['nucleus'].values
        larmord_predCS = data_ext['predCS'].values
        # print larmord_predCS
        for res in range(len(larmord_resid)):
            #keypred = str(str(larmord_resid[res])+":"+larmord_resname[res]+":"+larmord_nucleus[res]).strip(sep=b'')
            keypred = str(str(larmord_resid[res])+":"+larmord_resname[res]+":"+larmord_nucleus[res])
            # generate keys to self.predictedCS
            predCS[keypred] = larmord_predCS[res]
        return predCS

    def prepare_bme_files(self):
        states = [i for i in range(1,1+cmd.count_states(self.pymol_sel.get()))]
        predCS_data_ext = self.predCS_data_ext[self.predCS_data_ext['state'].isin(states)]

        # remove data flagged as outlier
        self.expCS_data_ext['keys'] = self.expCS_data_ext.agg('{0[resid]}:{0[resname]}:{0[nucleus]}'.format, axis=1)
        expCS_data_ext = self.expCS_data_ext[~self.expCS_data_ext['keys'].isin(self.outlier_keys)]

        # filter chemical shifts if needed
        if self.larmord_error_sel == 'proton':
            expCS_data_ext = expCS_data_ext[self.expCS_data_ext['nucleus'].isin(self.proton_list)]
        if self.larmord_error_sel == 'carbon':
            expCS_data_ext = expCS_data_ext[self.expCS_data_ext['nucleus'].isin(self.carbon_list)]
        if self.larmord_error_sel == 'nitrogen':
            expCS_data_ext = expCS_data_ext[self.expCS_data_ext['nucleus'].isin(self.nitrogen_list)]

        # merge predCS, expCS, and errors
        self.mergedCS = predCS_data_ext.merge(expCS_data_ext, on = ['resname', 'resid', 'nucleus'])
        self.mergedCS = self.mergedCS.merge(self.larmord_acc_ext, on = ['resname', 'nucleus'])
        self.mergedCS = self.mergedCS.sort_values(['nucleus', 'resid', 'state'], ascending=[None,True,True])
        print(self.mergedCS.head())

        # shape into matrices
        nrow = len(self.mergedCS.state.unique())
        self.cs_matrix = pd.DataFrame(self.mergedCS.predCS.values.reshape(-1, nrow).transpose())
        self.sim_matrix = pd.DataFrame(self.mergedCS.expCS.values.reshape(-1, nrow).transpose()).transpose()
        self.error_matrix = pd.DataFrame(self.mergedCS.MAE.values.reshape(-1, nrow).transpose()).transpose()

        self.sim_fn = None
        sim_os_fh, self.sim_fn = tempfile.mkstemp(suffix='.txt') # file os handle, file name
        os.close(sim_os_fh)

        self.exp_fn = None
        exp_os_fh, self.exp_fn = tempfile.mkstemp(suffix='.txt') # file os handle, file name
        os.close(exp_os_fh)

        # rewrite out files
        self.cs_matrix.to_csv(self.sim_fn, sep = ' ', header = None)
        pd.DataFrame([self.sim_matrix[0], self.error_matrix[0]]).transpose().to_csv(self.exp_fn, sep = ' ', index_label = '#', header = ['DATA=JCOUPLINGS', 'PRIOR=GAUSS'])

    def prepare_for_clustering(self):
        states = [i for i in range(1,1+cmd.count_states(self.pymol_sel.get()))]

        # get predicted chemical shifts
        predCS_data_ext = self.predCS_data_ext[self.predCS_data_ext['state'].isin(states)]

        # merge data
        self.clusterCS = predCS_data_ext.merge(self.larmord_acc_ext, on = ['resname', 'nucleus'])
        self.clusterCS = self.clusterCS.sort_values(['nucleus', 'resid', 'state'], ascending=[None,True,True])

        # shape into matrices
        nrow = len(self.clusterCS.state.unique())
        self.clusterCS = self.clusterCS.predCS.values.reshape(-1, nrow).transpose()

    def run_clustering(self):
        X = self.clusterCS
        X = StandardScaler().fit_transform(self.clusterCS)
        m = KMeans(n_clusters=self.min_size.get(), random_state=0).fit(X)
        m.fit(X)
        clusters = list(m.labels_)
        clusters.insert(0,0)
        return(clusters)

    def prepare_file_for_analysis(self, predictor, sel_name, objname):
        one_obj_sel = '%s and %s' % (sel_name, objname)
        pdb_fn = None
        pdb_os_fh, pdb_fn = tempfile.mkstemp(suffix='.pdb') # file os handle, file name
        os.close(pdb_os_fh)
        cmd.save(filename=pdb_fn, selection=one_obj_sel)
        if VERBOSE:
            print('Selection %s saved to %s.' % (one_obj_sel, pdb_fn))
        if pdb_fn is None:
            print('WARNING: %s has no pdb file to work on!' % predictor)
            return None
        print('Started Running %s for %s ...' % (predictor, one_obj_sel,))
        larmord_sse_dict = {}
        larmord_tmpout_os_fh, larmord_tmpout_fn = tempfile.mkstemp(suffix='.larmord')
        os.close(larmord_tmpout_os_fh)
        return pdb_fn, larmord_tmpout_fn


    def runAnalysisOneState(self, sel_name, objname):
        """
        For one state, compute chemical shift using Larmord or larmorca or load cs from file.
        Save the outcome (predited CS) in file: larmord_tmpout_fn, and call parse_larmor_output function to save all data in self.predictedCS.
        @param selname: selection name
        @param objname: state number
        @param type: string, int
        """
        if (self.get_shifts_from_larmord):
            pdb_fn, larmord_tmpout_fn = self.prepare_file_for_analysis('LARMORD', sel_name, objname)
            larmord_cmd = '%s/larmord -cutoff 15.0 -parmfile %s -reffile %s %s | awk \'{print $3, $4, $5, $7, $6}\' > %s' % (self.larmord_bin.get(), self.larmord_para.get(), self.larmord_ref.get(), pdb_fn, larmord_tmpout_fn)
            self.predCS_data_ext = larmord_tmpout_fn
            os.system(larmord_cmd)
        if (self.get_shifts_from_larmorca):
            pdb_fn, lamorca_tmpout_fn = self.prepare_file_for_analysis('LARMORCA', sel_name, objname)
            # ['resid', 'resname', 'nucleus', 'junk', 'predCS']
            # 1 164 N LEU 121.8 121.581 0 ID
            lamorca_cmd = '%s/larmorca %s | awk \'{print $2, $4, $3, $5, $6}\' > %s' % (self.larmorca_bin.get(), pdb_fn, lamorca_tmpout_fn)
            self.predCS_data_ext = lamorca_tmpout_fn
            os.system(lamorca_cmd)
        # Here call function self.parse_larmor_output and save predicted CS data to a global dictionary self.predCS
        if self.get_shifts_from_larmord:
            predCS = self.parse_larmor_output(larmord_tmpout_fn)
            self.predictedCS.append(predCS.copy())
        if self.get_shifts_from_larmorca:
            predCS = self.parse_larmor_output(lamorca_tmpout_fn)
            self.predictedCS.append(predCS.copy())
        else:
            if self.check_file(self.cs2.get()):
                predCS = self.parse_larmord_dataframe(data_ext = self.predCS_data_ext, state = cmd.get_state())
                self.predictedCS.append(predCS.copy())
            else:
                self.print_file_error(self.cs2.get())
                return False
        return True

    ## Functions related to self.runCompare()
    def runBME(self):
        bmea = bme.Reweight(verbose=VERBOSE)

        # load data
        bmea.load(self.exp_fn, self.sim_fn)

        # do optimization using theta=2
        bme_converged = False
        theta=1.0
        while(not bme_converged):
            try:
                chi2_before,chi2_after, srel = bmea.optimize(theta=theta)
                bme_converged = True
            except:
                theta+=1.0
        self.w_opt = list(np.round_(bmea.get_weights(), 5))
        self.w_opt.insert(0,0)
        return(theta)

    def root_mean_square_error(self,x,y,mae):
        # Calculate the root mean square error of two vectors x,y
        # Avoid dependence on python packages
        mse = 0.0 # mean square error
        N = len(x)
        if len(x)!=len(y) or len(x)!=len(mae):
            return False
        for k in range(len(x)):
            mse = mse + (x[k]-y[k])*(x[k]-y[k])/(mae[k]*mae[k])
            #print(mse, x[k], y[k], mae[k])
        try:
            rmse = np.sqrt(mse/N)
        except:
            return 0
        return rmse

    def confidence(self,x, y, mae):
        # Calculate the confidence PYM
        # see: A Bayesian approach to NMR crystal structure determination (10.1039/C9CP04489B)
        prob = 1.0
        N = len(x)
        if len(x)!=len(y) or len(x)!=len(mae):
            return False
        for k in range(len(x)):
            error = (x[k]-y[k])/mae[k]
            prefactor =  (1/np.sqrt(2*np.pi*mae[k]*mae[k]))
            prob_k =  prefactor * np.exp(-0.5*error*error)
            prob = prob*prob_k
        print("My confidence is: %s"%prob)
        return prob

    def getMAE(self, error, nOfAtoms):
      '''
      error: total error
      nOfAtoms: total number of atoms
      '''
      try:
        mae = error/nOfAtoms
      except:
        return 0.0
      return mae

    def computePearson(self, nucleus, output_nucleus, nnucleus, state_number):
        # Compute Pearson coefficients between exp. and pred. data
        # seperately for different nucleus within the same group
        # e.g. Compute Pearson coefficients for N1 and N3 seperately
        # and take a weighted average to obtain Pearson coefficients for Nitrogen
        nuclei_list = eval('self.'+nucleus+'_list')
        pearson = 0.0
        p_value = 0.0
        for nuclei in nuclei_list:
            nnuclei = 0
            list_expCS_nuclei = []
            list_predCS_nuclei = []
            for key in output_nucleus.keys():
                resid, resname, nucleus = key.split(":")
                predCS = self.predictedCS[state_number - 1][key]
                expCS = self.measuredCS[key]
                if nucleus == nuclei:
                    nnuclei += 1
                    list_expCS_nuclei.append(expCS)
                    list_predCS_nuclei.append(predCS)
            if nnuclei > 1:
                pears, p_value = stats.pearsonr(list_predCS_nuclei, list_expCS_nuclei)
                pearson += (pears * nnuclei) / nnucleus
        return pearson

    def runCompareOneState(self, state_number):
        """
        For one state, compute error and correlation coefficients between measured and predicted chemical shifts.
        @param state_number: specify the state to compare
        @param type: int
        @output total MAE, proton MAE, carbon MAE, nitrogen MAE, error dict
        @output type: float, float, float, dictionary
        """
        # reset b-factors
        cmd.alter("n. *", "b=0.0")

        total_error = 0.0
        proton_error = 0.0
        carbon_error = 0.0
        nitrogen_error = 0.0
        ntotal = 0
        nprotons = 0
        ncarbons = 0
        nnitrogens = 0
        predCS = 0.0
        output_total = {}
        output_carbon = {}
        output_proton = {}
        output_nitrogen = {}
        list_expCS = []
        list_expCS_carbon = []
        list_expCS_proton = []
        list_expCS_nitrogen = []
        list_predCS = []
        list_predCS_carbon = []
        list_predCS_proton = []
        list_predCS_nitrogen = []
        list_mae = []
        list_mae_carbon = []
        list_mae_proton = []
        list_mae_nitrogen = []

        # load measured Chemical Shifts data
        for key in self.predictedCS[state_number - 1].keys():
            ch = " "
            resid, resname, nucleus = key.split(":")
            predCS = self.predictedCS[state_number - 1][key]
            k1 = (ch, resid, resname, nucleus)
            k2 = str(resname+":"+nucleus).strip()
            # try to get mae (expected error)
            try:
                mae = self.mae[k2]
            except:
                mae = 1.0
                print(k2)

            # try to get measured chemical shifts for each predicted value
            try:
                expCS = self.measuredCS[key]

                # ignore predictions that exhibit errors that are greater than mae * self.larmord_outlier_threshold
                if nucleus in self.carbon_list:
                    if (np.abs(predCS - expCS - self.larmord_carbon_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        self.outlier_keys.append(key)
                        continue
                if nucleus in self.proton_list:
                    if (np.abs(predCS - expCS - self.larmord_proton_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        self.outlier_keys.append(key)
                        continue
                if nucleus in self.nitrogen_list:
                    if (np.abs(predCS - expCS - self.larmord_nitrogen_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        self.outlier_keys.append(key)
                        continue

                list_expCS.append(expCS)
                list_predCS.append(predCS)
            except:
                continue
            list_mae.append(mae)

            if nucleus in self.carbon_list:
                error = (predCS - expCS - self.larmord_carbon_offset.get())/mae
                carbon_error += np.abs(error)
                ncarbons += 1
                list_expCS_carbon.append(expCS)
                list_predCS_carbon.append(predCS)
                list_mae_carbon.append(mae)
                output_carbon[key] = np.abs(error)
            elif nucleus in self.proton_list:
                error = (predCS - expCS - self.larmord_proton_offset.get())/mae
                proton_error += np.abs(error)
                nprotons += 1
                list_expCS_proton.append(expCS)
                list_predCS_proton.append(predCS)
                list_mae_proton.append(mae)
                output_proton[key] = np.abs(error)
            elif nucleus in self.nitrogen_list:
                error = (predCS - expCS - self.larmord_nitrogen_offset.get())/mae
                nitrogen_error += np.abs(error)
                nnitrogens += 1
                list_expCS_nitrogen.append(expCS)
                list_predCS_nitrogen.append(predCS)
                list_mae_nitrogen.append(mae)
                output_nitrogen[key] = np.abs(error)
            else:
                continue
            output_total[key] = np.abs(error)
            ntotal +=1
            total_error += np.abs(error)
            cmd.alter("resi %s and n. %s"%(resid, nucleus), "b=%s"%error)

        # all shifts
        pearson = self.computePearson('total', output_total, ntotal, state_number)
        self.Pearson_coef.append(pearson)
        RMSE = self.root_mean_square_error(list_predCS, list_expCS, list_mae)
        self.RMSE_coef.append(RMSE)

        PYM = self.confidence(list_predCS, list_expCS, list_mae)
        self.PYM_coef.append(PYM)

        # carbon shifts
        pearson = self.computePearson('carbon', output_carbon, ncarbons, state_number)
        self.Pearson_carbon.append(pearson)
        RMSE = self.root_mean_square_error(list_predCS_carbon, list_expCS_carbon, list_mae_carbon)
        self.RMSE_carbon.append(RMSE)
        PYM = self.confidence(list_predCS_carbon, list_expCS_carbon, list_mae_carbon)
        self.PYM_carbon.append(PYM)

        # proton shifts
        pearson = self.computePearson('proton', output_proton, nprotons, state_number)
        self.Pearson_proton.append(pearson)
        RMSE = self.root_mean_square_error(list_predCS_proton, list_expCS_proton, list_mae_proton)
        self.RMSE_proton.append(RMSE)
        PYM = self.confidence(list_predCS_proton, list_expCS_proton, list_mae_proton)
        self.PYM_proton.append(PYM)

        # nitrogen shifts
        pearson = self.computePearson('nitrogen', output_nitrogen, nnitrogens, state_number)
        self.Pearson_nitrogen.append(pearson)
        RMSE = self.root_mean_square_error(list_predCS_nitrogen, list_expCS_nitrogen, list_mae_nitrogen)
        self.RMSE_nitrogen.append(RMSE)
        PYM = self.confidence(list_predCS_nitrogen, list_expCS_nitrogen, list_mae_nitrogen)
        self.PYM_nitrogen.append(PYM)


        print('Complete calculating error for state %d' %state_number)
        return self.getMAE(total_error, ntotal), self.getMAE(carbon_error, ncarbons), self.getMAE(proton_error, nprotons), self.getMAE(nitrogen_error, nnitrogens), output_total, output_carbon, output_proton, output_nitrogen

    ## Functions related to self.runSort()
    def render_larmord_errors(self, objname, type, scale):
        """
        """
        # set VDW for all RNA atoms to 0
        cmd.alter(objname+" and n. *", "vdw = 0")
        cmd.rebuild()
        # reinitialize b-factors
        stored.b = []
        # iterated over selection and store b-factors
        if(type == "proton"):
            sele = "("+objname+" and n. " + "+".join(self.proton_list)+ ")"
            stored.min = 0.0
        if(type == "carbon"):
            sele = "("+objname+" and n. " + "+".join(self.carbon_list)+ ")"
            stored.min = 0.0
        if(type == "nitrogen"):
            sele = "("+objname+" and n. " + "+".join(self.nitrogen_list)+ ")"
            stored.min = 0.0
        if(type == "all"):
            sele = "("+objname+" and n. " + "+".join(self.total_list)+ ")"
            stored.min = 0.0
        cmd.iterate (sele, "stored.b.append(b)")
        # alter vDW radius on b-factor
        cmd.alter(sele, "vdw = (%s * b)"%scale)
        # rebuild
        cmd.rebuild()

    def render_one_object(self, obj_name, error_sel, error_scale, error_color):
        # modify sphere radius in portion to error
        self.render_larmord_errors(obj_name, error_sel, error_scale)
        # display cartoon and spheres and then color based on b-factors
        cmd.set( "cartoon_ring_mode", 3)
        cmd.set("cartoon_tube_radius", 1)
        cmd.set("cartoon_ring_width", 0.3)
        cmd.show("sticks", "not backbone")
        cmd.set("stick_radius", 0.5)
        cmd.show("cartoon", obj_name)
        cmd.show("spheres", obj_name)
        cmd.spectrum("b", error_color, obj_name)

    def runRender(self):
        """
        Show all states in Pymol in sorted order given by self.best_model_indices
        """
        # delete old results
        self.sel_obj_list = []
        pdb_fn = None
        sel_name= None
        sel = self.pymol_sel.get()
        error_color = self.larmord_error_color.get()
        error_sel = self.larmord_error_sel
        error_scale = self.larmord_error_scale.get()
        error_height = self.larmord_error_height.get()
        error_lsize = self.larmord_error_lsize.get()

        # checking selection
        sel_name = self.check_selection(sel)
        if (not sel_name):
          return False

        # disable parent object
        cmd.enable("all")
        cmd.disable(sel_name)
        cmd.delete("best*")

        # turn grid mode ON
        cmd.set("grid_mode", 1)

        counter = 1
        for s in self.best_model_indices[:self.larmord_ndisplayed.get()]:
            # get indices of and load best models
            obj_new = "best_"+str(s)
            pdb_fn = sel_name+"_"+str(s)+".pdb"
            cmd.load(pdb_fn, obj_new)

            # render errors
            self.render_one_object(obj_new, error_sel, error_scale, error_color)
            cmd.set("grid_slot", counter, obj_new)
            counter += 1

        # Create groups
        cmd.group('best', 'best_'+'*')
        cmd.disable("best")
        return True

    def showModels(self):
        """
        Run render and enable group 'best' to show best models in a list
        """
        self.runRender()
        self.currentstate = 0
        cmd.disable("all")
        self.resetCSTable()
        cmd.enable("best*")

    ## Primary buttons ('Predict or Load Chemical Shifts', 'Compare', 'Sort') callback functions
    def runAnalysis(self):
        """
        Callback function for 'Run' button in 'Compute or Load Shifts' tab
        Run larmord over all states within one object
        Call self.runAnalysisOneState in each loop and loop over all states
        """
        # delete old results
        self.sel_obj_list = []
        pdb_fn = None
        sel_name= None
        sel = self.pymol_sel.get()
        self.reset_errorTable()

        # Disable all the buttons while running analysis
        self.disableAll()

        # reset
        self.get_shifts_from_file_larmor = False

        # reset predicted and measured chemical shifts(if ever loaded) before running analysis
        self.reset_predictedCS()

        # reset outliers keys
        self.outlier_keys = []

        # checking selection
        sel_name = self.check_selection(sel)
        if (not sel_name):
          return False

        # each object in the selection is treated as an independent struc
        cmd.enable(sel_name) # enable selection
        objlist = cmd.get_object_list(sel_name)
        for objname in objlist:
            self.m.set(0)
            self.sel_obj_list.append('%s and %s' % (sel_name, objname))
            if not self.get_shifts_from_file:
                #if number of states is greater than 1
                if self.PSICO and cmd.count_states(objname) > 1:
                    self.get_shifts_from_file_larmor = True

                    pdb_fn = None
                    pdb_os_fh, pdb_fn = tempfile.mkstemp(suffix='.pdb') # file os handle, file name
                    os.close(pdb_os_fh)

                    dcd_fn = None
                    pdb_os_fh, dcd_fn = tempfile.mkstemp(suffix='.dcd') # file os handle, file name
                    os.close(pdb_os_fh)

                    larmord_tmpout_fn = None
                    pdb_os_fh, larmord_tmpout_fn = tempfile.mkstemp(suffix='.larmord') # file os handle, file name
                    os.close(pdb_os_fh)

                    psico.exporting.save_traj(filename = dcd_fn, selection = objname)
                    cmd.save(filename = pdb_fn, selection = objname, state = 1)

                    if self.get_shifts_from_larmord:
                        larmor_cmd = '%s/larmord -cutoff 15.0 -parmfile %s -reffile %s -trj %s %s | awk \'{print $2+1, $3, $4, $5, $6, $9}\' > %s' % (self.larmord_bin.get(), self.larmord_para.get(), self.larmord_ref.get(), dcd_fn, pdb_fn, larmord_tmpout_fn)
                    if self.get_shifts_from_larmorca:
                        larmor_cmd = '%s/larmorca  -trj %s %s | awk \'{print $1, $2, $4, $3, $6, $8}\' > %s' % (self.larmorca_bin.get(), dcd_fn, pdb_fn, larmord_tmpout_fn)
                    print(larmor_cmd)
                    os.system(larmor_cmd)
                    self.cs2_internal = larmord_tmpout_fn
                    self.cs2.set(larmord_tmpout_fn)
                    self.predCS_data_ext = pd.read_csv(larmord_tmpout_fn, sep = '\s+', names = ['state', 'resid', 'resname', 'nucleus', 'predCS', 'id'])
                    #print(self.predCS_data_ext.shape)
            # load external file
            if self.get_shifts_from_file:
                self.predCS_data_ext = pd.read_csv(self.cs2.get(), sep = '\s+', names = ['state', 'resid', 'resname', 'nucleus', 'predCS', 'id'])
                print("from external file")
            if self.get_shifts_from_file_larmor:
                self.predCS_data_ext = pd.read_csv(self.cs2_internal, sep = '\s+', names = ['state', 'resid', 'resname', 'nucleus', 'predCS', 'id'])
                print("from internal file")
            if self.get_shifts_from_file or self.get_shifts_from_file_larmor:
                self.get_shifts_from_larmord = False
                self.get_shifts_from_larmorca = False
                # check for mismatch between  number of states in external file and number states in pymol object
                if self.check_states_mismatch():
                    self.print_mismatch_error()
                    return False

            # loop over states
            for a in range(1,1+cmd.count_states("(all)")):
              cmd.frame(a)
              self.runAnalysisOneState(sel_name, objname)
              progress = float(a)/float((cmd.count_states("(all)")))
              self.m.set(progress)
        self.analyzeButton.button(0).config(state = 'normal')
        self.tableButton.button(0).config(state = 'normal')
        self.notebook.nextpage()
        return True

    def runCompare(self):
        """
        Callback functions for 'Compare' button in 'table' tab
        Compare measured CS value with predicted CS values for all states within one object
        """
        # delete old results
        self.sel_obj_list = []
        pdb_fn = None
        sel_name= None
        sel = self.pymol_sel.get()
        self.reset_errorTable()

        # Cleanup
        cmd.delete("label*")
        cmd.delete("best*")
        cmd.delete("model*")
        cmd.enable(sel)

        self.disableAll()

        lowerLimit = 10 * cmd.count_states("(all)")
        # Initialize error coefficients
        self.Pearson_coef = [0]
        self.Kendall_coef = [0]
        self.Spearman_coef = [0]
        self.RMSE_coef = [0]
        self.PYM_coef = [0]

        self.Pearson_carbon = [0]
        self.Kendall_carbon = [0]
        self.Spearman_carbon = [0]
        self.RMSE_carbon = [0]
        self.PYM_carbon = [0]

        self.Pearson_proton = [0]
        self.Kendall_proton = [0]
        self.Spearman_proton = [0]
        self.RMSE_proton = [0]
        self.PYM_proton = [0]

        self.Pearson_nitrogen = [0]
        self.Kendall_nitrogen = [0]
        self.Spearman_nitrogen = [0]
        self.RMSE_nitrogen = [0]
        self.PYM_nitrogen = [0]

        self.load_measuredCS()
        self.conv_resname_format()

        # load MAEs and maybe prep for BME
        if self.weighted_errors:
            self.load_MAE()
        else:
            self.reset_MAE()

        # checking selection
        sel_name = self.check_selection(sel)
        if (not sel_name):
          return False

        # each object in the selection is treated as an independent struc
        objlist = cmd.get_object_list(sel_name)
        for objname in objlist:
            # reset progress bar
            self.m.set(0)
            temp = '%s and %s' % (sel_name, objname)
            self.null_w_opt = [1.0]
            self.null_clusters = [-1]
            self.sel_obj_list.append(temp)
            for a in range(1,1+cmd.count_states("(all)")):
                cmd.frame(a)
                print('now in state %d' %cmd.get_state())
                self.total_error[a], self.carbon_error[a], self.proton_error[a], self.nitrogen_error[a], self.larmord_error_all[a], self.larmord_error_carbon[a], self.larmord_error_proton[a], self.larmord_error_nitrogen[a] = self.runCompareOneState(a)
                # write output PDB with the b-factors in the
                obj_new = objname+"_"+str(a)
                pdb_fn = obj_new+".pdb"
                cmd.save(filename=pdb_fn, selection=temp)
                progress = float(a)/float((cmd.count_states("(all)")))
                self.null_w_opt.append(1.0)
                self.null_clusters.append(-1)
                self.m.set(progress)

        # initial best indices
        self.runSort()

        self.analyzeButton.button(0).config(state = 'normal')
        self.tableButton.button(1).config(state = 'normal')
        self.tableButton.button(0).config(state = 'normal')
        self.save_CStable.button(3).config(state = 'normal')
        return True

    def runSort(self):
        """
        Callback function of 'Sort' button in 'table' tab
        'Sort' -> self.runSort -> self.showModels -> self.runRender
        """
        self.resetCSTable()
        # setup files for clustering and weighting

        # finalize outlier key
        RUN_BME_ = self.BME and cmd.count_states(self.pymol_sel.get()) > 1 and self.check_file(self.cs.get())
        self.outlier_keys = list(set(self.outlier_keys))
        if self.SKLEARN and cmd.count_states(self.pymol_sel.get()) > 1:
            self.prepare_for_clustering()
        if RUN_BME_:
            self.prepare_bme_files()
        if self.SKLEARN  and cmd.count_states(self.pymol_sel.get()) > 1:
            self.clusters = self.run_clustering()
            print("SHOULD BE ABLE TO RUN SKLEARN: %s"%self.clusters)
        else:
            self.clusters = self.null_clusters
        if RUN_BME_:
            theta = self.runBME()
            print(self.w_opt)
            print("BME %s %s %s"%(len(self.w_opt), np.sum(self.w_opt), theta))
        else:
            self.w_opt = self.null_w_opt

        metric = self.larmord_error_metric
        if metric in ['MAE']:
            self.sort_error()
        else:
            self.sort_coef()
        self.printError(self.best_model_indices)
        self.showModels()
        return True

    ## Methods and callback functions related to error table
    def reset_errorTable(self):
        """
        """
        self.error_table.delete(2,'end')
        self.save_CStable.button(0).config(state = 'disabled')
        self.save_CStable.button(1).config(state = 'disabled')
        self.resetCSTable()

    def printError(self, orderList):
        """
        Print MAE and correlation coefficients between predicted and measured CS for each state in the following format:
        state || MAE || P coef || RMSE || PYM
        print error and coefs row by row in chosen order given by orderList
        @param orderList: control the order that the table to be printed in
        @param type: list
        """
        self.reset_errorTable()
        nucleus = self.larmord_error_sel
        if nucleus == 'proton':
            self.table_header = 'state_number proton_error Pearson_proton RMSE_proton w_opt clusters'
        if nucleus == 'carbon':
            self.table_header = 'state_number carbon_error Pearson_carbon RMSE_carbon w_opt clusters'
        if nucleus == 'nitrogen':
            self.table_header = 'state_number nitrogen_error Pearson_nitrogen RMSE_nitrogen w_opt clusters'
        if nucleus == 'all':
            self.table_header = 'state_number total_error Pearson_coef RMSE_coef w_opt clusters'
        self.table_header = str.split(self.table_header)

        #Initialize state_number array
        self.state_number = orderList

        # Print error data row by row
        for row in range(cmd.count_states("(all)")):
            state = orderList[row]
            data = '%-5s ' % (str(state),)
            dataline = ' ' + data + '   '
            for column in range(1, len(self.table_header)):
                value = eval('self.' + self.table_header[column] + '[' + str(state) + ']')
                data = str(value)[:5]
                data = '%-5s ' % (data,)
                dataline = dataline + data + '   '
            self.error_table.insert('end', dataline)
        self.save_CStable.button(0).config(state = 'normal')
        self.save_CStable.button(1).config(state = 'normal')

    def saveErrortable(self):
        objname = self.pymol_sel.get()
        filename = asksaveasfilename(initialdir = "saved_data/", initialfile = str('error_table_'+objname+'.txt'))
        try:
            errorTable = open(filename,'w')
        except:
            return False
        metric_temp = self.larmord_error_metric
        nucleus_temp = self.larmord_error_sel

        # Loop over all possible sorting methods and save error table
        # When metric = 'state'
        tag = 'Sorted by state:\n'
        errorTable.write(tag)
        self.sort_state_number()
        value = self.error_table.get(1,'end')
        for dataline in value:
            dataline = dataline + '\n'
            errorTable.write(dataline)

        # When metric in {everthing else}
        for self.larmord_error_metric in ['MAE', 'pearson', 'RMSE', 'weight', 'clusters',]:
            metric = self.larmord_error_metric
            for self.larmord_error_sel in ['proton', 'carbon', 'nitrogen', 'all']:
                nucleus = self.larmord_error_sel
                if metric in ['MAE']:
                    self.sort_error()
                else:
                    self.sort_coef()
                self.printError(self.best_model_indices)
                tag = str('Sorted by ' + nucleus+ ' ' + metric + ':\n')
                errorTable.write(tag)
                value = self.error_table.get(1,'end')
                for dataline in value:
                    dataline = dataline + '\n'
                    errorTable.write(dataline)
        errorTable.close()

        # Reset global parameters to initial condition
        self.larmord_error_metric = metric_temp
        self.larmord_error_sel = nucleus_temp
        if self.currentstate!=0:
            self.showCStable()
        else:
            self.runSort()

    ### Method(s) and callback function(s) related to selection in error table
    def sele_from_table(self, event):
        """
        Render the corresponding state when double click on a certain row in table
        """
        # clear old results in CS table and remove hightlight
        self.resetCSTable()

        error_color = self.larmord_error_color.get()
        error_scale = self.larmord_error_scale.get()
        error_height = self.larmord_error_height.get()
        error_lsize = self.larmord_error_lsize.get()
        sel_name = self.pymol_sel.get()
        error_sel = self.larmord_error_sel

        # disable parent object
        cmd.disable("all")
        #cmd.disable(sel_name)
        #cmd.disable("best*")
        # turn ON grid mode
        cmd.set("grid_mode", 0)

        s = self.error_table.curselection()
        selection = []
        cmd.delete("focus")
        for index in s:
            a = int(index) - 2 #offset for the two header rows
            if a >= 0:
                a = self.best_model_indices[a]
                selection.append(a)
                objname = sel_name +"_"+str(a)
                cmd.load(objname+".pdb", objname)
                # load error values for selected states
                total_error = self.total_error[a]
                proton_error = self.proton_error[a]
                carbon_error = self.carbon_error[a]
                nitrogen_error = self.nitrogen_error[a]
                self.render_one_object(objname, error_sel , error_scale, error_color)
                self.currentstate = a
                self.showCStable()
                cmd.enable(objname)
                cmd.orient(objname)
        cmd.group("focus", sel_name+"_*")

    ### Methods and callback functions related to sorting error table
    def sort_state_number(self):
        orderList = []
        for s in range(1, 1+cmd.count_states()):
            orderList.append(s)
        self.best_model_indices = orderList
        self.printError(orderList)
        return True

    def sort_error(self):
        nucleus = self.larmord_error_sel
        if nucleus == 'proton':
            self.best_model_indices = np.argsort(list(self.proton_error.values())) + 1
        if nucleus == 'carbon':
            self.best_model_indices = np.argsort(list(self.carbon_error.values())) + 1
        if nucleus == 'nitrogen':
            self.best_model_indices = np.argsort(list(self.nitrogen_error.values())) + 1
        if nucleus == 'all':
            self.best_model_indices = np.argsort(list(self.total_error.values())) + 1
        return True

    def sort_coef(self):
        metric = self.larmord_error_metric
        if metric in ['pearson']:
            self.sort_Pearson_coef()
        if metric in ['kendall']:
            self.sort_Kendall_coef()
        if metric in ['spearman']:
            self.sort_Spearman_coef()
        if metric in ['RMSE']:
            self.sort_RMSE_coef()
        if metric in ['PYM']:
            self.sort_PYM_coef()
        if metric in ['state']:
            self.sort_state_number()
        if metric in ['weight']:
            self.sort_weight()
        if metric in ['clusters']:
            self.sort_clusters()

    def sort_weight(self):
        self.best_model_indices = np.argsort(self.w_opt)[::-1]
        return True

    def sort_clusters(self):
        clusters = []
        for c in range(1, len(self.clusters)):
            clusters.append(self.clusters[c])
        self.best_model_indices = np.argsort(clusters)[::-1] + 1
        print(self.best_model_indices)
        return True

    def sort_Pearson_coef(self):
        nucleus = self.larmord_error_sel
        if nucleus == 'proton':
            self.best_model_indices = np.argsort(self.Pearson_proton)[::-1]
        if nucleus == 'carbon':
            self.best_model_indices = np.argsort(self.Pearson_carbon)[::-1]
        if nucleus == 'nitrogen':
            self.best_model_indices = np.argsort(self.Pearson_nitrogen)[::-1]
        if nucleus == 'all':
            self.best_model_indices = np.argsort(self.Pearson_coef)[::-1]
        return True

    def sort_Kendall_coef(self):
        nucleus = self.larmord_error_sel
        if nucleus == 'proton':
            self.best_model_indices = np.argsort(self.Kendall_proton)[::-1]
        if nucleus == 'carbon':
            self.best_model_indices = np.argsort(self.Kendall_carbon)[::-1]
        if nucleus == 'nitrogen':
            self.best_model_indices = np.argsort(self.Kendall_nitrogen)[::-1]
        if nucleus == 'all':
            self.best_model_indices = np.argsort(self.Kendall_coef)[::-1]
        return True

    def sort_Spearman_coef(self):
        nucleus = self.larmord_error_sel
        if nucleus == 'proton':
            self.best_model_indices = np.argsort(self.Spearman_proton)[::-1]
        if nucleus == 'carbon':
            self.best_model_indices = np.argsort(self.Spearman_carbon)[::-1]
        if nucleus == 'nitrogen':
            self.best_model_indices = np.argsort(self.Spearman_nitrogen)[::-1]
        if nucleus == 'all':
            self.best_model_indices = np.argsort(self.Spearman_coef)[::-1]
        return True

    def sort_RMSE_coef(self):
        nucleus = self.larmord_error_sel
        if nucleus == 'proton':
            self.best_model_indices = np.argsort(self.RMSE_proton)[1:]
        if nucleus == 'carbon':
            self.best_model_indices = np.argsort(self.RMSE_carbon)[1:]
        if nucleus == 'nitrogen':
            self.best_model_indices = np.argsort(self.RMSE_nitrogen)[1:]
        if nucleus == 'all':
            self.best_model_indices = np.argsort(self.RMSE_coef)[1:]
        return True

    def sort_PYM_coef(self):
        nucleus = self.larmord_error_sel
        if nucleus == 'proton':
            self.best_model_indices = np.argsort(self.PYM_proton)[1:]
        if nucleus == 'carbon':
            self.best_model_indices = np.argsort(self.PYM_carbon)[1:]
        if nucleus == 'nitrogen':
            self.best_model_indices = np.argsort(self.PYM_nitrogen)[1:]
        if nucleus == 'all':
            self.best_model_indices = np.argsort(self.PYM_coef)[1:]
        return True

    ## Methods and callback functions related to chemical shifts table
    def resetCSTable(self):
        """
        """
        self.CS_table.delete(2,'end')
        cmd.delete('highlight')
        for button in range(self.sort_CStable.numbuttons()):
            self.sort_CStable.button(button).config(state = 'disabled')
        self.save_CStable.button(2).config(state = 'disabled')

    def sortResidual(self):
        state = self.currentstate
        temp_dict = {}
        temp_dict = eval('self.larmord_error_'+self.larmord_error_sel+'['+str(state)+'].copy()')
        for key in temp_dict.keys():
            if key[1] == ':':
                new_key = '0' + key
                temp_dict[new_key] = temp_dict.pop(key)
        self.sorted_residual = sorted(temp_dict.keys())

    def showCStable(self):
        """
        Sort chemical shift data table in ascending res id order
        """
        temp_dict = {}
        self.resetCSTable()
        state = self.currentstate
        self.sortResidual()
        self.printCStable(state)

    def sort_CStable_error(self):
        """
        Sort chemical shift data table of each residue for selected state in error descending order
        """
        temp_dict = {}
        self.resetCSTable()
        state = self.currentstate
        temp_dict = eval('self.larmord_error_'+self.larmord_error_sel+'['+str(state)+'].copy()')
        self.sorted_residual = sorted(temp_dict, key=temp_dict.get, reverse = True)
        self.printCStable(state)

    def saveCStableOnestate(self):
        """
        save chemical shift table for selected one state
        """
        state = self.currentstate
        sel_name = self.pymol_sel.get()
        filename = asksaveasfilename(initialdir = "saved_data/", initialfile = str('predicted_CS_table_' + sel_name + '_' + str(state)+'.dat'))
        try:
            CStable = open(filename,'w')
        except:
            return False
        for key in self.sorted_residual:
            if key[0] == '0':
                key = key[1:]
            predCS = self.predictedCS[state - 1][key]
            resid, resname, nucleus = key.split(':')

            # try to get mae (expected error)
            k2 = str(resname+":"+nucleus).strip()
            try:
                mae = self.mae[k2]
            except:
                mae = 1.0
            try:
                error = self.larmord_error_all[state][key]
                expCS = self.measuredCS[key]
                # ignore predictions that exhibit errors that are greater than mae * self.larmord_outlier_threshold
                if nucleus in self.carbon_list:
                    if (np.abs(predCS - expCS - self.larmord_carbon_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        continue
                if nucleus in self.proton_list:
                    if (np.abs(predCS - expCS - self.larmord_proton_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        continue
                if nucleus in self.nitrogen_list:
                    if (np.abs(predCS - expCS - self.larmord_nitrogen_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        continue
            except:
                continue
            dataline = "{:<7} {:<7} {:<7} {:<7,.3f}\n".format(resname, resid, nucleus, predCS)
            CStable.write(dataline)
        CStable.close()
        msg = 'Predicted chemical shifts and error has been saved to %s' %filename
        return True

    def saveCStable(self):
        flag = 'bound'
        objname = self.pymol_sel.get()
        filename = asksaveasfilename(initialdir = "saved_data/", initialfile = str('predicted_CS_table_' + objname + '.txt'))
        try:
            CStable = open(filename,'w')
        except:
            return False
        try:
            self.sorted_residual
        except:
            temp_dict = {}
            temp_dict = eval('self.larmord_error_'+self.larmord_error_sel+'[1].copy()')
            self.sorted_residual = sorted(temp_dict, key=temp_dict.get, reverse = True)
        for state in range(1,1+cmd.count_states("(all)")):
            for key in self.sorted_residual:
                if key[0] == '0':
                    key = key[1:]
                predCS = self.predictedCS[state - 1][key]
                resid, resname, nucleus = key.split(':')

                # try to get mae (expected error)
                k2 = str(resname+":"+nucleus).strip()
                try:
                    mae = self.mae[k2]
                except:
                    mae = 1.0

                try:
                    expCS = self.measuredCS[key]
                    # ignore predictions that exhibit errors that are greater than mae * self.larmord_outlier_threshold
                    if nucleus in self.carbon_list:
                        if (np.abs(predCS - expCS - self.larmord_carbon_offset.get())) > mae * self.larmord_outlier_threshold.get():
                            continue
                    if nucleus in self.proton_list:
                        if (np.abs(predCS - expCS - self.larmord_proton_offset.get())) > mae * self.larmord_outlier_threshold.get():
                            continue
                    if nucleus in self.nitrogen_list:
                        if (np.abs(predCS - expCS - self.larmord_nitrogen_offset.get())) > mae * self.larmord_outlier_threshold.get():
                            continue
                    expCS = "%.3f" %expCS
                    predCS = "%.3f" %predCS
                except:
                    continue
                try:
                    error = self.larmord_error_all[state][key]
                except:
                    continue
                dataline = "{:<3} {:<5} {:<3} {:<5} {:<7} {:<7} {:<7}\n".format(state, resname, resid, nucleus, predCS, expCS, error)
                CStable.write(dataline)
        CStable.close()
        return True
	
    def saveSession(self):
        """
        Save session
        """
        objname = self.pymol_sel.get()
        filename = asksaveasfilename(initialdir = "saved_data/", initialfile = str('save_session_'+objname+".pse"))
        cmd.save(filename, format = 'pse')

    def printCStable(self,state):
        """
        Print chemical shift data table in the order of self.sorted_residual
        """
        for button in range(self.sort_CStable.numbuttons()):
            self.sort_CStable.button(button).config(state = 'normal')
        self.save_CStable.button(2).config(state = 'normal')

        for key in self.sorted_residual:
            if key[0] == '0':
                key = key[1:]
            predCS = self.predictedCS[state - 1][key]
            resid, resname, nucleus = key.split(':')

            # try to get mae (expected error)
            k2 = str(resname+":"+nucleus).strip()
            try:
                mae = self.mae[k2]
            except:
                mae = 1.0
            try:
                error = self.larmord_error_all[state][key]
                expCS = self.measuredCS[key]

                # ignore predictions that exhibit errors that are greater than mae * self.larmord_outlier_threshold
                if nucleus in self.carbon_list:
                    if (np.abs(predCS - expCS - self.larmord_carbon_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        continue
                if nucleus in self.proton_list:
                    if (np.abs(predCS - expCS - self.larmord_proton_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        continue
                if nucleus in self.nitrogen_list:
                    if (np.abs(predCS - expCS - self.larmord_nitrogen_offset.get())) > mae * self.larmord_outlier_threshold.get():
                        continue
            except:
                continue
            dataline = "{:<7} {:<6} {:<6} {:<6,.1f} {:<8,.1f} {:<6,.1f}".format(resname, resid, nucleus, expCS, predCS, error)
            dataline = ' ' + dataline
            self.CS_table.insert('end', dataline)

    ### Methods and callback functions related to selection in CS table
    def highlight_selection(self, sele, color="yellow", scale = 0.5, transparency=0.5):
        """ function that highlights a selection
            e.g. highlight_selection("resi 1 and name H1' and test_1", scale=0.5, color="yellow")
            Taken from pymol.py by Dr.Frank
        """
        scale = 1.2 * self.larmord_error_scale.get()
        # remove old hightlight
        cmd.delete("highlight")
        # get original view
        v = cmd.get_view()
        # make copy of sele
        cmd.create("highlight", sele)
        # increase size
        cmd.alter("highlight", "vdw = (%s * b)"%(scale))
        # color
        cmd.color(color, "highlight")
        cmd.show("spheres", "highlight")
        cmd.set("sphere_transparency", transparency, "highlight")
        # reset view
        cmd.set_view(v)

    def seleRes(self, event):
        # Double-click callback function in chemical shift table
        state = self.currentstate
        if state == 0:
            return False
        sel = self.pymol_sel.get()
        res = self.CS_table.curselection()
        list = self.sorted_residual

        for index in res:
            a = int(index) - 2 #offset for the two header rows
            if a >= 0:
                key = list[a]
                if key[0] == '0':
                    key = key[1:]
                resid, resname, nuclei = key.split(':')
                sele = "resi %s and name %s and %s_%s" %(str(int(resid)), nuclei, sel, str(state))
                print("You select %s of state %d" %(key, state))
                self.highlight_selection(sele)
        return True

    ## File check, selection check and others
    # Check the validation of user provided file and other user inputs before every running
    # to avoid waste of computing resources
    def check_file(self, file):
        """
        Check files
        @param file: the file to be checked
        @param type: string
        """
        try:
            if os.path.getsize(file) > 0:
              return True
            else:
              return False
        except:
            return False

    def print_mismatch_error(self):
        err_msg = 'Mismatch # states in %s !=  # states in %s.' %(self.cs2.get(), self.pymol_sel.get())
        print('ERROR: %s' % (err_msg,))
        tkMessageBox.showinfo(title='ERROR', message=err_msg)

    def check_states_mismatch(self):
        return(len(self.predCS_data_ext.state.unique()) != cmd.count_states(self.pymol_sel.get()))

    def print_file_error(self, file):
        err_msg = 'The %s does not exist or is empty.' %(file,)
        print('ERROR: %s' % (err_msg,))
        tkMessageBox.showinfo(title='ERROR', message=err_msg)

    def check_selection(self, sel):
        if len(sel) > 0: # if any pymol selection/object is specified
            all_sel_names = cmd.get_names('all') # get names of all selections
            if sel in all_sel_names:
                if cmd.count_atoms(sel) == 0:
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print('ERROR: %s' % (err_msg,))
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = sel
                    return sel_name
            # no selection/object with the input name is found
            # we assume either a single-word selector or
            # some other selection-expression is uesd
            else:
                print('The selection/object you specified is not found.')
                print('Your input will be interpreted as a selection-expression.')
                tmpsel = self.randomSeleName(prefix='your_sele_')
                cmd.select(tmpsel,sel)
                if cmd.count_atoms(tmpsel) == 0:
                    cmd.delete(tmpsel)
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print('ERROR: %s' % (err_msg,))
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = tmpsel
                    return sel_name

        else:   # what structure do you want Larmord to work on?
            err_msg = 'No PyMOL selection/object specified!'
            print('ERROR: %s' % (err_msg,))
            tkMessageBox.showinfo(title='ERROR', message=err_msg)
            return False

    def randomSeleName(self, prefix='sele', suffix=''):
        """ generate a random selection name.
        """
        sel_list = cmd.get_names('all')
        sel_dict = dict(zip(sel_list, range(len(sel_list))))
        sel_name = '%s%d%s' % (prefix, random.randint(1000,9999), suffix)
        while(sel_name in sel_dict):
            sel_name = '%s%d%s' % (prefix, random.randint(1000,9999), suffix)
        return sel_name

    ## Primary Execute function
    def execute(self, butcmd):
        """ Run the cmd represented by the botton clicked by user.
        """
        if butcmd == 'OK':
            print('is everything OK?')
        elif butcmd == 'Run':
            rtn = self.runAnalysis()
            if rtn and VERBOSE:
                 print('Done with Larmord!')
        elif butcmd == 'Compare Chemcial Shifts':
            rtn = self.runCompare()
            if rtn and VERBOSE:
                 print('Done comparing chemical shifts!')
        elif butcmd == 'Sort':
            rtn = self.runSort()
            if rtn and VERBOSE:
                 print('Done rendering chemical shift errors!')
        elif butcmd == 'Exit':
            print('Exiting PyShifts Plugin ...')
            if __name__ == '__main__':
                self.parent.destroy()
            else:
                self.dialog.withdraw()
            print('Done.')
        else:
            print('Exiting PyShifts Plugin because of unknown button click ...')
            self.dialog.withdraw()
            print('Done.')

    def quit(self):
        self.dialog.destroy()

#############################################
#
# Create demo in root window for testing.
#
##############################################
if __name__ == '__main__':

    class App:
        def my_show(self,*args,**kwargs):
            pass

    app = App()
    app.root = tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('It seems to work!')

    widget = PyShiftsPlugin(app)
    app.root.mainloop()
