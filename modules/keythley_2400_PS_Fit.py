'''
 MIT License

Copyright (c) 2024 Tanay Dey
... [License text preserved] ...
Author: Dr. Tanay Dey
Present Affiliation: SAHA INSTITUTE OF NUCLEAR PHYSICS (SINP), KOLKATA, INDIA
'''

import pyvisa as visa
import time
import matplotlib.pyplot as plt
from tkinter import *
import tkinter as Tk
from tkinter import ttk
from tkinter import messagebox as msg
import os
from tkinter import filedialog
import sys
import numpy as np
import numpy.linalg as npl
import matplotlib as matplotlib
import matplotlib.patches as patches
import copy
import glob as glob
from matplotlib.patches import Polygon
import requests
import io
from datetime import date, datetime
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from io import StringIO
from contextlib import redirect_stdout
import logging
from PIL import Image, ImageTk
from matplotlib.widgets import Slider
import screeninfo
from tkinter.font import Font
import re
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator
from matplotlib import style
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
from serial.tools import list_ports
import matplotlib.path as mpath

import math
import serial
import socket
import pickle
import tempfile
import serial.tools.list_ports
import pandas as pd
from scipy.optimize import curve_fit

# --- DPI AWARENESS FOR WINDOWS ---
try:
    from ctypes import windll
    windll.shcore.SetProcessDpiAwareness(1)
except:
    pass

matplotlib.use("TkAgg")

# ==========================================
#  BREAKDOWN ANALYSIS HELPER FUNCTIONS
# ==========================================
def dinu_eq8_model(V, V_bd, V_cr, p, A, leak_a, leak_b):
    I_leak = np.exp(leak_a * V + leak_b)
    I_aval = np.zeros_like(V)
    mask = (V > V_bd) & (V < V_cr)
    
    if np.any(mask):
        dV = V[mask] - V_bd 
        gain = dV
        prob = 1 - np.exp(-p * dV)
        afterpulse = (V_cr - V_bd) / (V_cr - V[mask])
        I_aval[mask] = A * gain * prob * afterpulse
        
    return I_leak + I_aval

def fit_wrapper(V, *args):
    model_I = dinu_eq8_model(V, *args)
    model_I = np.where(model_I <= 1e-13, 1e-13, model_I)
    return np.log(model_I)

def find_vbd_derivative(voltage, current):
    current_safe = np.where(current <= 1e-9, 1e-9, current)
    log_I = np.log(current_safe)
    log_derivative = np.gradient(log_I, voltage)
    
    search_mask = voltage > 20
    if np.sum(search_mask) > 5:
        masked_deriv = log_derivative[search_mask]
        masked_volt = voltage[search_mask]
        peak_idx = np.argmax(masked_deriv)
        v_bd = masked_volt[peak_idx]
    else:
        peak_idx = np.argmax(log_derivative)
        v_bd = voltage[peak_idx]
        
    return v_bd
        
def optimize_fit(voltage, current, v_bd_guess, user_params=None):
    if len(voltage) < 5:
        return np.zeros(6), False
    print("max(d(log (I)/dV :: )",v_bd_guess)
    mask = voltage < v_bd_guess
    V_cut = voltage[mask]
    I_cut = current[mask]

# Step 2: midpoint indices
    mid = len(V_cut) // 2
    V1, V2 = V_cut[mid-1], V_cut[mid]
    I1, I2 = I_cut[mid-1], I_cut[mid]

# Step 3: compute a and b
    a = (np.log(I2) - np.log(I1)) / (V2 - V1)
    b = np.log(I1) - a * V1
    if np.max(voltage)> v_bd_guess:
       V_cr=np.max(voltage)
    else:
       V_cr=v_bd_guess+10
    if user_params:
        p0 = [
            user_params.get('v_bd', v_bd_guess),
            user_params.get('v_cr', V_cr),
            user_params.get('p', 10.0),
            user_params.get('A', np.max(current)),
            user_params.get('leak_a', a),
            user_params.get('leak_b', b)
        ]
    else:
        pre_bd_mask = (voltage < (v_bd_guess - 2)) & (current > 0)
        leak_a, leak_b = 0.1, -5.0 
        
        if np.sum(pre_bd_mask) > 3:
            try:
                log_I_pre = np.log(current[pre_bd_mask])
                coeffs = np.polyfit(voltage[pre_bd_mask], log_I_pre, 1)
                leak_a = coeffs[0]
                leak_b = coeffs[1]
            except:
                pass

        v_cr_guess = max(max(voltage) + 5.0, v_bd_guess + 10.0)
        p0 = [v_bd_guess, v_cr_guess, 1.0, 0.5, leak_a, leak_b]

    vbd_center = p0[0]
    bounds = (
        [vbd_center - 10.0, vbd_center + 0.1, 0.01, 0.0, -10.0, -50.0],
        [vbd_center + 5.0,  200.0,            10.0, 1e5, 10.0,  10.0] 
    )

    current_safe = np.where(current <= 0, 1e-13, current)
    log_current_data = np.log(current_safe)
    
    try:
        popt, pcov = curve_fit(fit_wrapper, voltage, log_current_data, p0=p0, bounds=bounds, maxfev=10000)
        return popt, True
    except Exception as e:
        print(f"Fit failed: {e}")


        return np.zeros(6), False        



class DraggableAnnotation:
    def __init__(self, annotation):
        self.annotation = annotation
        self.got_artist = False
        self.canvas = self.annotation.figure.canvas
        # Connect to matplotlib events
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.annotation.axes: return
        contains, _ = self.annotation.contains(event)
        if not contains: return
        self.got_artist = True

    def on_motion(self, event):
        if not self.got_artist or event.inaxes != self.annotation.axes: return
        # Update text position to mouse location
        self.annotation.xytext = (event.xdata, event.ydata)
        self.canvas.draw_idle()

    def on_release(self, event):
        self.got_artist = False


# ==========================================
#  MAIN GUI CLASS
# ==========================================

class KeithleyGUI:
    def __init__(self):
        # Initialize Main Window
        self.window = Tk.Tk()
        self.window.title('SINP')
        self.window.configure(bg="red")
        self.analysis_artists = []      
        
        # --- Initialize Global State Variables ---
        self.counter = 0
        self.current_threshold = 1
        self.max_voltage = 1500
        self.plt_flag = 1
        
        # Plotting Lists
        self.xp = []
        self.yp = []
        self.ypp = []
        self.xp_ap = []
        self.temp_arr = []
        self.humid_arr = []
        self.time_arr = []
       
        # Simulation Data
        #self.voltage_array_sim = np.array([ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0, 25.5, 26.0, 26.5, 27.0, 27.5, 28.0, 28.5, 29.0, 29.5, 30.0])
        
        #self.current_array_sim = np.array([ 0.7104839000000001, 1.2014749999999998, 0.8774472999999999, 1.940415, 1.705872, 1.274704, 2.186259, 1.96661, 1.894825, 2.018717, 2.064218, 1.941652, 2.573901, 1.9488300000000005, 2.58141, 2.818302, 2.691768, 2.866244, 2.385318, 2.806337, 3.081027, 2.648411, 2.857701, 2.644859, 2.8332390000000003, 2.913645, 3.332809, 3.310944, 3.468049, 3.617178, 3.733407, 3.737988, 3.736461, 4.008461, 3.881432, 3.770413, 4.475728, 3.956766, 4.240734, 4.340819, 4.632298, 4.449694999999999, 4.246217, 5.70826, 3.521805, 3.3010390000000003, 7.359289, 5.850054, 5.308605, 25.1172, 192.9801, 517.5822000000001, 1007.169, 1585.119, 2397.39, 3382.047, 4665.905, 6149.853, 8115.172, 10512.3])

        self.voltage_array_sim = np.array([
     0.4999859, 0.999952, 1.499947, 1.999928, 2.4998, 2.999674, 
    3.499631, 3.999846, 4.499805, 4.999751, 5.499687, 5.999937, 6.499873, 
    6.999689, 7.499645, 7.999514, 8.499705, 8.999655, 9.499553, 9.999434, 
    10.49975, 10.9996, 11.49948, 11.99941, 12.49961, 12.99953, 13.49946, 
    13.99931, 14.49956, 14.99948, 15.49936, 15.99926, 16.49949, 16.9995, 
    17.49935, 17.99931, 18.4996, 18.99953, 19.4994, 19.99928, 20.49948, 
    20.99946, 21.47763, 21.986, 22.48145, 22.97448, 23.43093, 23.97558, 
    24.47595, 24.93115, 25.48596, 25.98343, 26.4994, 26.98108, 27.47599, 
    27.99158, 28.48512, 28.98525, 29.49876, 29.99419
])

# Current (nA) - Assuming units are nA based on previous context
        self.current_array_sim = np.array([
    0.7972846, 1.071266, 0.9152755, 0.3773804, 3.883199, 1.607349, 
    1.948123, 1.853812, 2.26815, 2.152795, 1.635359, 2.737103, 2.327063, 
    2.674965, 2.541881, 2.638132, 2.742049, 2.592083, 2.4066, 2.865247, 
    2.800307, 2.840492, 3.457065, 3.302141, 3.093801, 3.626013, 4.004491, 
    3.907786, 3.897441, 4.115551, 4.529222, 4.330207, 4.619445, 4.696694, 
    4.699408, 4.793971, 4.82329, 5.057756, 4.620898, 4.868884, 4.642936, 
    4.796926, 6.430016, 4.66504, 5.922447, 8.494838, 4.520731, 2.77246, 
    5.683654, 21.71631, 178.5566, 440.3047, 823.2871, 1340.467, 1975.572, 
    2828.87, 3826.657, 5109.14, 6757.934, 8820.594
])
        # Instrument & Control Flags
        self.C_ucell=0
        self.rm = None
        self.instrument = None
        self.search_flag = 0
        self.run_time_flag = 0
        self.pause_plot = 0
        self.figure_canvas = None
        self.canvas_analysis = None
        self.plot1 = None
        self.plot2 = None
        self.plot3 = None
        self.plot4 = None 
        self.plot5 = None 
        self.plot6 = None 
        self.warn_flag = 0
        self.legn_flag = 0
        self.start_vol = 0
        self.end_vol = 0
        self.step_vol = 0
        self.time_delay = 0
        self.curr_th = 0
        self.legend1 = None
        self.sim_flag = 0
        self.run_flag = 0
        self.stop_flag = 0
        self.run_index = 0
        self.baud_rate = 9600
        self.ard_flag = 0
        self.rmp_dwn_flag = 0
        self.all_ports = {''}
        self.ii = 0
        self.run_init_flg = 0 
        self.polarinit = 0
        
        # Arduino specific
        self.ser = None
        self.label8 = None
        self.arduino_port_list = None
        
        # --- Initialize Tkinter Variables ---
        self.datapath = StringVar()
        self.p_address = StringVar()
        self.module_name = StringVar(value="Sipm V-I Characteristic Test")
        self.current_th = StringVar(value="10000")
        self.Nmeas = StringVar(value="5")
        self.start_voltage = StringVar()
        self.end_voltage = StringVar()
        self.step_voltage = StringVar(value='0.5')
        self.down_step_voltage = StringVar(value='1')
        self.delay_time = StringVar(value='1')
        self.current_datetimes = StringVar()
        self.arduino_ports = StringVar(value="Choose Option")
        self.p_reading = StringVar(value='VOLTAGE::  V\nCURRENT:: μA')
        self.single_voltage = StringVar()
        self.user_answer = StringVar()
        self.var = IntVar()
        self.scale_var = StringVar(value="log")
        
        self.calc_vbd_var = Tk.BooleanVar(value=True)
        self.show_geiger_var = Tk.BooleanVar(value=True)
        self.show_dcr_var = Tk.BooleanVar(value=False)
        self.user_fit_params = {}

        # --- Load Resources (Images) ---
        self.photo1 = self.light_images('r1.png')
        self.photo2 = self.light_images('y1.png')
        self.photo3 = self.light_images('g1.png')
        
         # Post Process Variables
        self.selected_log_file = StringVar()
        self.selected_log_file.set("")
        self.x_start_var = DoubleVar(value=0)
        self.x_end_var = DoubleVar(value=10)
        self.log_scale_var = Tk.BooleanVar(value=False)
        self.show_temp_hum_var = Tk.BooleanVar(value=False)
        self.current_unit_var = Tk.StringVar(value="µA")
        self.breakdown_voltage_var = Tk.BooleanVar(value=False)
        self.giger_prob_var = Tk.BooleanVar(value =False)
        self.voltage_min = 0
        self.voltage_max = 10
        self.set_title = Tk.StringVar()
        self.set_title.set("V-I Characteristic Post Processed Plot")

        self.set_ovv = Tk.StringVar()
        self.set_ovv.set("2.5")
        
        self.CURRENT_SCALE = {"A": 1, "mA": 1e3, "µA": 1e6, "nA": 1e9}
        self.post_canvas = None

        # --- Build GUI ---
        self.setup_gui()
        
        try:
            self.window.state('zoomed')
        except:
            w, h = self.window.winfo_screenwidth(), self.window.winfo_screenheight()
            self.window.geometry(f"{w}x{h}")

        self.window.protocol("WM_DELETE_WINDOW", self.exits)

    def _configure_styles(self):
        self.style = ttk.Style()
        self.style.theme_use('clam')  
        self.colors = {
            'bg_main': '#F4F6F9', 'bg_sidebar': '#2C3E50', 'fg_sidebar': '#ECF0F1',
            'accent': '#2980B9', 'accent_hover': '#3498DB', 'warning': '#E74C3C',
            'success': '#27AE60', 'text_dark': '#2C3E50'
        }
        self.style.configure('Main.TFrame', background=self.colors['bg_main'])
        self.style.configure('Sidebar.TFrame', background=self.colors['bg_sidebar'])
        self.style.configure('Sidebar.TLabel', background=self.colors['bg_sidebar'], foreground=self.colors['fg_sidebar'], font=('Segoe UI', 10))
        self.style.configure('Header.TLabel', background=self.colors['bg_sidebar'], foreground=self.colors['fg_sidebar'], font=('Segoe UI', 12, 'bold'))
        self.style.configure('Panel.TLabel', background=self.colors['bg_main'], foreground=self.colors['text_dark'], font=('Segoe UI', 10))
        self.style.configure('Group.TLabelframe', background=self.colors['bg_sidebar'], foreground=self.colors['fg_sidebar'])
        self.style.configure('Group.TLabelframe.Label', background=self.colors['bg_sidebar'], foreground='#BDC3C7', font=('Segoe UI', 9, 'bold'))
        self.style.configure('Action.TButton', font=('Segoe UI', 10, 'bold'), padding=6)
        self.style.map('Action.TButton', background=[('active', self.colors['accent_hover'])], foreground=[('active', 'black')])
        self.style.configure('TNotebook', background=self.colors['bg_main'], tabposition='n')
        self.style.configure('TNotebook.Tab', padding=[12, 4], font=('Segoe UI', 10))

    def run(self):
        self.window.mainloop()

    # ==========================================
    #  GUI SETUP (SCROLLBAR + TABS)
    # ==========================================
    def setup_gui(self):
        """Builds the GUI layout ensuring it fits on screen without forced scrolling."""
        self._configure_styles()
        self.window.configure(bg=self.colors['bg_main'])
        
        # 1. Main Container
        self.master_frame = ttk.Frame(self.window, style='Main.TFrame')
        self.master_frame.pack(fill=Tk.BOTH, expand=True)

        # Tabs
        self.tab_control = ttk.Notebook(self.master_frame)
        self.tab1 = ttk.Frame(self.tab_control, style='Main.TFrame')
        self.tab3 = ttk.Frame(self.tab_control, style='Main.TFrame')  # Post Process Tab
        
        self.tab_control.add(self.tab1, text=' IV / HV Characterization ')
        self.tab_control.add(self.tab3, text=' Post Process ')
        self.tab_control.pack(fill=Tk.BOTH, expand=True, padx=5, pady=5)

        # =========================================================================
        # TAB 1: SIDEBAR (SCROLLABLE & COMPACT)
        # =========================================================================
        sidebar_frame = Tk.Frame(self.tab1, bg=self.colors['bg_sidebar'], width=300)
        sidebar_frame.pack(side=Tk.LEFT, fill=Tk.Y, expand=False)
        
        # Scrollbar Configuration
        sidebar_canvas = Tk.Canvas(sidebar_frame, bg=self.colors['bg_sidebar'], width=350, highlightthickness=0)
        sidebar_scroll = ttk.Scrollbar(sidebar_frame, orient="vertical", command=sidebar_canvas.yview)
        sidebar_canvas.configure(yscrollcommand=sidebar_scroll.set)
        
        sidebar_scroll.pack(side=Tk.RIGHT, fill=Tk.Y)
        sidebar_canvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True)
        
        # Inner Frame
        self.button_frame = ttk.Frame(sidebar_canvas, style='Sidebar.TFrame')
        canvas_window = sidebar_canvas.create_window((0, 0), window=self.button_frame, anchor="nw")

        # Bindings
        def configure_scroll_region(event):
            sidebar_canvas.configure(scrollregion=sidebar_canvas.bbox("all"))
        def configure_window_width(event):
            sidebar_canvas.itemconfig(canvas_window, width=event.width)
        self.button_frame.bind("<Configure>", configure_scroll_region)
        sidebar_canvas.bind("<Configure>", configure_window_width)

        def _on_mousewheel(event):
            sidebar_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        sidebar_canvas.bind("<Enter>", lambda _: sidebar_canvas.bind_all("<MouseWheel>", _on_mousewheel))
        sidebar_canvas.bind("<Leave>", lambda _: sidebar_canvas.unbind_all("<MouseWheel>"))

        # --- SIDEBAR WIDGETS ---
        ttk.Label(self.button_frame, text="KEITHLEY 2410 CONTROL", style='Header.TLabel', padding=(5, 5)).pack(fill=Tk.X)

        conn_group = ttk.LabelFrame(self.button_frame, text="CONNECTION", style='Group.TLabelframe', padding=5)
        conn_group.pack(fill=Tk.X, padx=5, pady=2, expand=True) 
        conn_inner = ttk.Frame(conn_group, style='Sidebar.TFrame')
        conn_inner.pack(fill=Tk.X)
        ttk.Label(conn_inner, text="Address:", style='Sidebar.TLabel').pack(side=Tk.LEFT)
        self.ps_address_screen = ttk.Entry(conn_inner, textvariable=self.p_address, width=15)
        self.ps_address_screen.pack(side=Tk.LEFT, fill=Tk.X, expand=True, padx=5)
        ttk.Button(conn_group, text="Connect", style='Action.TButton', command=self.search_or_set).pack(fill=Tk.X, pady=(2,0))

        exp_group = ttk.LabelFrame(self.button_frame, text="CONFIG", style='Group.TLabelframe', padding=5)
        exp_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        ttk.Label(exp_group, text="Module Name:", style='Sidebar.TLabel').pack(anchor='w')
        ttk.Entry(exp_group, textvariable=self.module_name).pack(fill=Tk.X, pady=(0, 5))
        mode_frame = ttk.Frame(exp_group, style='Sidebar.TFrame')
        mode_frame.pack(fill=Tk.X)
        Radiobutton(mode_frame, text='Positive IV', variable=self.user_answer, value='HV', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self.HVTEST).pack(side=Tk.LEFT, expand=True)
        Radiobutton(mode_frame, text='Negative IV', variable=self.user_answer, value='IV', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self.IVTEST).pack(side=Tk.LEFT, expand=True)

        volt_group = ttk.LabelFrame(self.button_frame, text="PARAMETERS (V)", style='Group.TLabelframe', padding=5)
        volt_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        for i in range(3): volt_group.columnconfigure(i, weight=1)
        ttk.Label(volt_group, text="Start (V):", style='Sidebar.TLabel').grid(row=0, column=0, sticky='e')
        ttk.Entry(volt_group, textvariable=self.start_voltage, width=6).grid(row=0, column=1, sticky='ew', padx=2)
        ttk.Label(volt_group, text="End (V):", style='Sidebar.TLabel').grid(row=0, column=2, sticky='e')
        ttk.Entry(volt_group, textvariable=self.end_voltage, width=6).grid(row=0, column=3, sticky='ew', padx=2)
        ttk.Label(volt_group, text="Up (V) :", style='Sidebar.TLabel').grid(row=1, column=0, sticky='e')
        ttk.Entry(volt_group, textvariable=self.step_voltage, width=6).grid(row=1, column=1, sticky='ew', padx=2)
        ttk.Label(volt_group, text="Down (V) :", style='Sidebar.TLabel').grid(row=1, column=2, sticky='e')
        ttk.Entry(volt_group, textvariable=self.down_step_voltage, width=6).grid(row=1, column=3, sticky='ew', padx=2)
        ttk.Label(volt_group, text="Delay (Sec) :", style='Sidebar.TLabel').grid(row=2, column=0, sticky='e')
        ttk.Entry(volt_group, textvariable=self.delay_time, width=6).grid(row=2, column=1, sticky='ew', padx=2)
        ttk.Label(volt_group, text="Curr Lim (uA):", style='Sidebar.TLabel').grid(row=2, column=2, sticky='e')
        ttk.Entry(volt_group, textvariable=self.current_th, width=6).grid(row=2, column=3, sticky='ew', padx=2)
        ttk.Label(volt_group, text="No. Meas Per Step:", style='Sidebar.TLabel').grid(row=3,columnspan=2, column=0, sticky='e')
        #ttk.Entry(volt_group, textvariable=self.Nmeas, width=6).grid(row=3, column=2,columnspan=2, sticky='ew', padx=2)
        #ttk.Entry(volt_group, textvariable=self.Nmeas, width=6).tap(lambda w: w.bind("<FocusOut>", lambda e: self.validate_and_run())).grid(row=3, column=2,columnspan=2, sticky='ew', padx=2)
        e = ttk.Entry(volt_group, textvariable=self.Nmeas, width=6); e.bind("<FocusOut>", lambda _: self.validate_and_run()); e.grid(row=3, column=2,columnspan=2, sticky='ew', padx=2)
        plot_group = ttk.LabelFrame(self.button_frame, text="PLOT", style='Group.TLabelframe', padding=5)
        plot_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        scale_frame = ttk.Frame(plot_group, style='Sidebar.TFrame')
        scale_frame.pack(fill=Tk.X)
        Radiobutton(scale_frame, text='Log', variable=self.scale_var, value='log', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self.change_scale).pack(side=Tk.LEFT, expand=True)
        Radiobutton(scale_frame, text='Linear', variable=self.scale_var, value='linear', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self.change_scale).pack(side=Tk.LEFT, expand=True)

        analysis_group = ttk.LabelFrame(self.button_frame, text="ANALYSIS", style='Group.TLabelframe', padding=5)
        analysis_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        Checkbutton(analysis_group, text="Breakdown V", variable=self.calc_vbd_var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar']).pack(anchor='w')
        Checkbutton(analysis_group, text="Geiger Prob", variable=self.show_geiger_var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar']).pack(anchor='w')
        Checkbutton(analysis_group, text="DCR", variable=self.show_dcr_var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'],command=self.open_dcr_window ).pack(anchor='w')

        ttk.Button(analysis_group, text="Fit Params", style='Action.TButton', command=self.open_param_window).pack(fill=Tk.X, pady=(2,0))

        manual_group = ttk.LabelFrame(self.button_frame, text="MANUAL", style='Group.TLabelframe', padding=5)
        manual_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        man_row = ttk.Frame(manual_group, style='Sidebar.TFrame')
        man_row.pack(fill=Tk.X)
        ttk.Entry(man_row, textvariable=self.single_voltage, width=8).pack(side=Tk.LEFT, fill=Tk.X, expand=True, padx=(0,5))
        ttk.Button(man_row, text="Set V", width=5, command=self.set_single_voltage).pack(side=Tk.LEFT, padx=1)
        ttk.Button(man_row, text="Zero", width=5, command=self.ramp_down_single_voltage).pack(side=Tk.LEFT, padx=1)

        self.env_group = ttk.LabelFrame(self.button_frame, text="ENV", style='Group.TLabelframe', padding=5)
        self.env_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        Checkbutton(self.env_group, text="Connect Arduino", variable=self.var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=lambda: self.check_button_clicked(self.var)).pack(anchor='w')

        act_frame = ttk.Frame(self.button_frame, style='Sidebar.TFrame')
        act_frame.pack(fill=Tk.X, padx=5, pady=5, side=Tk.BOTTOM)
        Button(act_frame, text="START TEST", bg=self.colors['success'], fg='white', font=('Segoe UI', 10, 'bold'), relief=Tk.FLAT, pady=2, command=self.start_process).pack(fill=Tk.X, pady=1)
        Button(act_frame, text="EXIT", bg=self.colors['warning'], fg='white', font=('Segoe UI', 9), relief=Tk.FLAT, pady=0, command=self.exits).pack(fill=Tk.X, pady=1)
        ttk.Label(act_frame, text="Designed By: Dr. Tanay Dey", style='Sidebar.TLabel', font=('Segoe UI', 7)).pack(pady=0)

        # =========================================================================
        # TAB 1: MAIN AREA (RIGHT SIDE)
        # =========================================================================
        right_panel = Tk.Frame(self.tab1, bg='white') 
        right_panel.pack(side=Tk.RIGHT, fill=Tk.BOTH, expand=True)

        # 1. Monitor Strip
        monitor_container = Tk.Frame(right_panel, bg='white')
        monitor_container.pack(side=Tk.TOP, fill=Tk.X, padx=10, pady=5)

        readout_frame = Tk.LabelFrame(monitor_container, text="Real-time Readings", bg='white', fg='#7f8c8d', font=("arial", 10))
        readout_frame.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True, padx=(0, 20))
        self.labels1 = Label(readout_frame, textvariable=self.p_reading, bg='white', fg='#2980B9', font=("Noto Sans", 16, 'bold'), justify=Tk.LEFT, anchor="w")
        self.labels1.pack(fill=Tk.BOTH, expand=True, padx=10, pady=5)

        ctrl_container = Tk.Frame(monitor_container, bg='white')
        ctrl_container.pack(side=Tk.RIGHT)
        self.image_label2 = ttk.Label(ctrl_container, image=self.photo1, background='white')
        self.image_label2.pack(side=Tk.LEFT, padx=10)
        self.pause = Button(ctrl_container, text='PAUSE', bg='#E0E0E0', relief=GROOVE, command=self.pause_plots)
        self.pause.pack(side=Tk.LEFT, padx=5)
        Button(ctrl_container, text='STOP', bg='#E74C3C', fg='white', relief=GROOVE, font=("arial", 10, "bold"), command=self.stop_run).pack(side=Tk.LEFT, padx=5)
        Button(ctrl_container, text='Simulate', bg='#E0E0E0', relief=GROOVE, command=self.simulation_run).pack(side=Tk.LEFT, padx=5)

        # ---------------------------------------------------------------------
        # 2. PLOT AREA WITH TABS (NEW)
        # ---------------------------------------------------------------------
        self.canvas_frame = Tk.Frame(right_panel, bg='white')
        self.canvas_frame.pack(side=Tk.TOP, fill=Tk.BOTH, expand=True, padx=10, pady=(0, 10))

        self.plot_notebook = ttk.Notebook(self.canvas_frame)
        self.plot_notebook.pack(fill=Tk.BOTH, expand=True)

        # Tab 1: Measurement
        self.tab_measure = ttk.Frame(self.plot_notebook, style='Main.TFrame')
        self.plot_notebook.add(self.tab_measure, text='  Measurement  ')

        # Tab 2: Analysis Result
        self.tab_analysis = ttk.Frame(self.plot_notebook, style='Main.TFrame')
        self.plot_notebook.add(self.tab_analysis, text='  Analysis Result  ')

        # Measurement Setup
        self.keithley_img_frame = ttk.Frame(self.tab_measure, style='Main.TFrame')
        self.keithley_img_frame.pack(fill=Tk.BOTH, expand=True)

        self.figure = plt.Figure(figsize=(5, 4), dpi=100) 
        self.figure.patch.set_facecolor(self.colors['bg_main'])
        self.ax = self.figure.add_subplot(111)
        self.plot1, = self.ax.plot([], [], 'o', color='#3498DB', markersize=4, label="Measured I-V")
        self.plot2, = self.ax.plot([], [], 'x', color='#E74C3C', markersize=4, label="Set I-V")
        self.plot3, = self.ax.plot([], [], 'b', linestyle='None', label="Limit")
        self.ax2 = self.ax.twinx()
        self.plot4, = self.ax2.plot([], [], 'ro', linestyle='None', label="Set voltage")
        self.plot5, = self.ax2.plot([], [], 'b:', label="Temp")
        self.plot6, = self.ax2.plot([], [], 'g-.', label="Humidity")
        self.ax.set_yscale(self.scale_var.get())
        self.ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        self.ax.set_facecolor('white')

        self.figure_canvas = FigureCanvasTkAgg(self.figure, master=self.tab_measure) 

        # Analysis Setup
        self.fig_analysis = plt.Figure(figsize=(5, 4), dpi=100)
        self.fig_analysis.patch.set_facecolor('white')
        self.canvas_analysis = FigureCanvasTkAgg(self.fig_analysis, master=self.tab_analysis)
        self.canvas_analysis.get_tk_widget().pack(fill=Tk.BOTH, expand=True)


        # Setup Tab 2 (Post Process)
        self._setup_post_process_tab()

    def open_param_window(self):
        top = Tk.Toplevel(self.window)
        top.title("Fit Initial Parameters")
        top.geometry("350x380")
        
        lbl_info = ttk.Label(top, text="Leave blank to use Auto-Guess", foreground="blue")
        lbl_info.pack(pady=5)
        
        fields = [
            ("V_bd Guess (V)", "v_bd"), ("V_cr Guess (V)", "v_cr"),
            ("Geiger Shape (p)", "p"), ("Amplitude (A)", "A"),
            ("Leakage Coeff (a)", "leak_a"), ("Leakage Offset (b)", "leak_b")
        ]
        
        self.entries = {}
        for text, key in fields:
            row = ttk.Frame(top)
            row.pack(fill=X, padx=10, pady=4)
            ttk.Label(row, text=text, width=18, anchor='w').pack(side=LEFT)
            ent = ttk.Entry(row)
            ent.pack(side=RIGHT, expand=True, fill=X)
            if key in self.user_fit_params: ent.insert(0, str(self.user_fit_params[key]))
            self.entries[key] = ent
            
        def save_params():
            for key, ent in self.entries.items():
                val = ent.get().strip()
                if val:
                    try:
                        self.user_fit_params[key] = float(val)
                    except ValueError:
                        msg.showerror("Error", f"Invalid value for {key}")
                        return
                else:
                    if key in self.user_fit_params: del self.user_fit_params[key]
            top.destroy()
            msg.showinfo("Success", "Parameters updated for next run.")
        ttk.Button(top, text="Save Parameters", command=save_params).pack(pady=15)

    def open_dcr_window(self):
        #if self.show_dcr_var.get() == 0:
        #    return

        win = Tk.Toplevel(self.window)
        win.title("DCR Parameters")
        win.geometry("420x300")
        win.configure(bg=self.colors['bg_main'])

        # Track whether Save was pressed
        self.dcr_saved = False

        def on_close():
            if not self.dcr_saved:
                self.show_dcr_var.set(0)
            win.destroy()

        win.protocol("WM_DELETE_WINDOW", on_close)

        # Title
        Tk.Label(
            win,
            text="DCR Calculation Parameters",
            font=("Arial", 12, "bold"),
            bg=self.colors['bg_main']
        ).pack(pady=10)

        # Input field
        frame = Tk.Frame(win, bg=self.colors['bg_main'])
        frame.pack(pady=10)

        Tk.Label(
            frame,
            text="Microcell Capacitance  Cμcell (F):",
            bg=self.colors['bg_main']
        ).grid(row=0, column=0, sticky="w", padx=5)

        self.c_ucell_entry = Tk.Entry(frame, width=20)
        self.c_ucell_entry.grid(row=0, column=1, padx=5)
        self.c_ucell_entry.insert(0, '1.79E-13')

        # Datasheet note
        note_text = (
            "NOTE:\n"
            "Cμcell = Ctotal / Ncells\n\n"
            "For SensL MicroFC-60035:\n"
            "Ctotal = 3400 pF\n"
            "Ncells = 18980\n\n"
            "⇒ Cμcell ≈ 1.79×10⁻13 F"
        )

        Tk.Label(
            win,
            text=note_text,
            justify="left",
            bg=self.colors['bg_main'],
            fg="gray20",
            wraplength=380
        ).pack(padx=10, pady=10)

        # Save button
        Tk.Button(
            win,
            text="Save",
            command=lambda: self.save_c_ucell(win)
        ).pack(pady=10)
    def save_c_ucell(self, win):
        try:
            self.C_ucell = float(self.c_ucell_entry.get())
            print("microcell value is :: ",self.C_ucell)
            self.dcr_saved = True
            win.destroy()
        except ValueError:
            Tk.messagebox.showerror(
                "Input Error",
                "Please enter Cμcell in Farads (e.g. 1.79e-13)"
            )
            self.show_dcr_var = Tk.BooleanVar(value=False)
            win.destroy()
                          
    def run_breakdown_analysis(self):
        if not self.xp or len(self.xp) < 5: return
        volts = np.array(self.xp)
        if self.sim_flag == 1:
        	currents_nA = np.array(self.yp)*1000000000.
        else:
                currents_nA = np.array(self.yp)	
        #currents_nA = currents_nA / 1000.0
        
        v_bd_deriv = find_vbd_derivative(volts, currents_nA)
        popt, success = optimize_fit(volts, currents_nA, v_bd_deriv, user_params=self.user_fit_params)
        self.plot_analysis_results(volts, currents_nA, v_bd_deriv, popt, success)

    '''def plot_analysis_results(self, volts, currents_nA, v_bd_deriv, popt, success):
        self.fig_analysis.clf() 
        self.fig_analysis.subplots_adjust(left=0.10, right=0.95, top=0.88, bottom=0.10, hspace=0.1)

        if hasattr(self, 'show_geiger_var'): show_geiger = self.show_geiger_var.get()
        else: show_geiger = True 
        
        if show_geiger:
            gs = GridSpec(2, 1, height_ratios=[3, 1]) 
            ax_iv = self.fig_analysis.add_subplot(gs[0])
            ax_prob = self.fig_analysis.add_subplot(gs[1], sharex=ax_iv)
        else:
            ax_iv = self.fig_analysis.add_subplot(111)
            ax_prob = None

        star = mpath.Path.unit_regular_star(6)
        circle = mpath.Path.unit_circle()
        cut_star = mpath.Path(vertices=np.concatenate([circle.vertices, star.vertices[::-1, ...]]), codes=np.concatenate([circle.codes, star.codes]))

        ax_iv.plot(volts, currents_nA, marker=cut_star, color='indigo', markersize=10, alpha=0.6, label="Measured Data", linestyle='None')
        
        if success:
            v_bd_fit = popt[0]
            y_val_nA = dinu_eq8_model(v_bd_fit, *popt) #* 1000.0
            idx = (np.abs(volts - v_bd_fit)).argmin()
            ax_iv.plot(v_bd_fit, y_val_nA, 'rx', markersize=10, markeredgewidth=2, label="Breakdown Point")

            ax_iv.annotate(f"Breakdown Point: {v_bd_fit:.2f}V", xy=(v_bd_fit, y_val_nA), xytext=(v_bd_fit - (max(volts)*0.15), currents_nA[idx]-y_val_nA/2), color='red', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), bbox=dict(boxstyle="round", fc="white", alpha=0.7), fontsize=13)
            #################################################################################
            y_val_nA_ov = dinu_eq8_model(v_bd_fit+2.7, *popt)
            
            ax_iv.plot(v_bd_fit+2.7, y_val_nA_ov, 'mP', markersize=10, markeredgewidth=2,label="Current at 2.7 Overvoltage")
            ax_iv.annotate(f"$V_{{bd}}+2.7$: {v_bd_fit+2.7:.2f} V\n I: {y_val_nA_ov:0.0f} nA", xy=(v_bd_fit+2.7, y_val_nA_ov), xytext=(v_bd_fit+2.7 -7, y_val_nA_ov-0.7*y_val_nA_ov), color='m', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), bbox=dict(boxstyle="round", fc="white", alpha=0.7), fontsize=13)
            #################################################################################
            v_smooth = np.linspace(min(volts), min(max(volts), popt[1]-0.1), 1000)
            i_fit_nA = dinu_eq8_model(v_smooth, *popt)
            #i_fit_nA = i_fit_uA #* 1000.0
            ax_iv.plot(v_smooth, i_fit_nA, 'g--', linewidth=2, label=f"Fit Model")
            ax_iv.axvline(v_bd_fit, color='blue', linestyle='--', alpha=0.5) 
            if ax_prob: ax_prob.axvline(v_bd_fit, color='blue', linestyle='--', alpha=0.5)
            
            if show_geiger and ax_prob:
                p_factor = popt[2]
                p_geiger = np.zeros_like(v_smooth)
                mask_aval = v_smooth > v_bd_fit
                if np.any(mask_aval): p_geiger[mask_aval] = 1 - np.exp(-p_factor * (v_smooth[mask_aval] - v_bd_fit))
                ax_prob.plot(v_smooth, p_geiger, 'b-', linewidth=2, label="Geiger Prob.")
                ax_prob.fill_between(v_smooth, p_geiger, color='blue', alpha=0.1)
                ax_prob.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
                ax_prob.set_ylabel("Geiger Prob.", fontweight='bold', color='blue', fontsize=13)
                ax_prob.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)
                ax_prob.set_ylim(-0.05, 1.1)
                ax_prob.grid(True, which='both', linestyle='--', alpha=0.5)
                formula_txt = r"$P_{Geiger} = 1 - e^{-p(V - V_{bd})}$"
                ax_prob.text(0.02, 0.6, formula_txt, transform=ax_prob.transAxes, 
                             fontsize=11, color='darkblue', fontweight='bold',
                             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85, ec="blue"))		
            if self.show_dcr_var.get()==0:
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}")
            else:
                if abs(self.C_ucell)>0:
                    DCR=popt[3]*1e-9/(self.C_ucell*1e3)
                else: DCR=0   
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}\n"+f"DCR : {DCR:0.3f} kHz")
    
                
                
                
            equation_latex = (r"$I_{tot} = I_{leak} + I_{aval}$" + "\n" + r"$I_{leak} = \exp(aV + b)$" + "\n" + r"$I_{aval} = A \cdot \Delta V \cdot (1 - e^{-p \Delta V}) \cdot \frac{V_{cr}-V_{bd}}{V_{cr}-V}$" + "\n")
            
            ax_iv.text(0.33, 0.96, equation_para, transform=ax_iv.transAxes, verticalalignment='top', fontsize=15, bbox=dict(boxstyle="round", fc="white", alpha=0.90, ec="#27AE60"), color="#2C3E50")
            ax_iv.text(0.01, 0.92, equation_latex, transform=ax_iv.transAxes, verticalalignment='top', fontsize=15, bbox=dict(boxstyle="round", fc="white", alpha=0.90, ec="green"), color="black")
                         
        ax_iv.set_ylabel("Current (nA)", fontweight='bold',fontsize=14)
        ax_iv.set_yscale(self.scale_var.get())
        if not show_geiger: ax_iv.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)   
        ax_iv.grid(True, which='both', linestyle='--', alpha=0.5)
        ax_iv.set_title(self.module_name.get(), pad=35) 
        ax_iv.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=4, mode="expand", borderaxespad=0., frameon=False, fontsize=14)
        
        self.canvas_analysis.draw()
        self.plot_notebook.select(self.tab_analysis)
        print(f"Analysis Complete. Vbd: {v_bd_deriv:.2f}V")'''
        
    def plot_analysis_results(self, volts, currents_nA, v_bd_deriv, popt, success):
        self.fig_analysis.clf() 
        self.fig_analysis.subplots_adjust(left=0.10, right=0.95, top=0.88, bottom=0.10, hspace=0.1)

        if hasattr(self, 'show_geiger_var'): show_geiger = self.show_geiger_var.get()
        else: show_geiger = True 
        
        if show_geiger:
            gs = GridSpec(2, 1, height_ratios=[3, 1]) 
            ax_iv = self.fig_analysis.add_subplot(gs[0])
            ax_prob = self.fig_analysis.add_subplot(gs[1], sharex=ax_iv)
        else:
            ax_iv = self.fig_analysis.add_subplot(111)
            ax_prob = None

        star = mpath.Path.unit_regular_star(6)
        circle = mpath.Path.unit_circle()
        cut_star = mpath.Path(vertices=np.concatenate([circle.vertices, star.vertices[::-1, ...]]), codes=np.concatenate([circle.codes, star.codes]))

        ax_iv.plot(volts, currents_nA, marker=cut_star, color='indigo', markersize=10, alpha=0.6, label="Measured Data", linestyle='None')
        
        if success:
            v_bd_fit = popt[0]
            y_val_nA = dinu_eq8_model(v_bd_fit, *popt) #* 1000.0
            idx = (np.abs(volts - v_bd_fit)).argmin()
            ax_iv.plot(v_bd_fit, y_val_nA, 'rx', markersize=10, markeredgewidth=2, label="Breakdown Point")

            ax_iv.annotate(f"Breakdown Point: {v_bd_fit:.2f}V", xy=(v_bd_fit, y_val_nA), xytext=(v_bd_fit - (max(volts)*0.15), currents_nA[idx]-y_val_nA/2), color='red', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), bbox=dict(boxstyle="round", fc="white", alpha=0.7), fontsize=13)
            
            #################################################################################
            overvol=2.5
            y_val_nA_ov = dinu_eq8_model(v_bd_fit+overvol, *popt)
            
            ax_iv.plot(v_bd_fit+overvol, y_val_nA_ov, 'mP', markersize=10, markeredgewidth=2,label=f"Current at {overvol:0.2f} Overvoltage")
            ax_iv.annotate(f"$V_{{bd}}+overvol$: {v_bd_fit+overvol:.2f} V\n I: {y_val_nA_ov:0.0f} nA", xy=(v_bd_fit+overvol, y_val_nA_ov), xytext=(v_bd_fit+overvol -7, y_val_nA_ov-0.7*y_val_nA_ov), color='m', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), bbox=dict(boxstyle="round", fc="white", alpha=0.7), fontsize=13)
            #################################################################################
            
            v_smooth = np.linspace(min(volts), min(max(volts), popt[1]-0.1), 1000)
            i_fit_nA = dinu_eq8_model(v_smooth, *popt)
            #i_fit_nA = i_fit_uA #* 1000.0
            ax_iv.plot(v_smooth, i_fit_nA, 'g--', linewidth=2, label=f"Fit Model")
            ax_iv.axvline(v_bd_fit, color='blue', linestyle='--', alpha=0.5) 
            if ax_prob: ax_prob.axvline(v_bd_fit, color='blue', linestyle='--', alpha=0.5)
            
            if show_geiger and ax_prob:
                p_factor = popt[2]
                p_geiger = np.zeros_like(v_smooth)
                mask_aval = v_smooth > v_bd_fit
                if np.any(mask_aval): p_geiger[mask_aval] = 1 - np.exp(-p_factor * (v_smooth[mask_aval] - v_bd_fit))
                ax_prob.plot(v_smooth, p_geiger, 'b-', linewidth=2, label="Geiger Prob.")
                ax_prob.fill_between(v_smooth, p_geiger, color='blue', alpha=0.1)
                ax_prob.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
                ax_prob.set_ylabel("Geiger Prob.", fontweight='bold', color='blue', fontsize=13)
                ax_prob.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)
                ax_prob.set_ylim(-0.05, 1.1)
                ax_prob.grid(True, which='both', linestyle='--', alpha=0.5)
                formula_txt = r"$P_{Geiger} = 1 - e^{-p(V - V_{bd})}$"
                ax_prob.text(0.02, 0.6, formula_txt, transform=ax_prob.transAxes, 
                             fontsize=11, color='darkblue', fontweight='bold',
                             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85, ec="blue"))		
            if self.show_dcr_var.get()==0:
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}")
            else:
                if abs(self.C_ucell)>0:
                    DCR=popt[3]*1e-9/(self.C_ucell*1e3)
                else: DCR=0   
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}\n"+f"DCR : {DCR:0.3f} kHz")
    
                
                
                
            equation_latex = (r"$I_{tot} = I_{leak} + I_{aval}$" + "\n" + r"$I_{leak} = \exp(aV + b)$" + "\n" + r"$I_{aval} = A \cdot \Delta V \cdot (1 - e^{-p \Delta V}) \cdot \frac{V_{cr}-V_{bd}}{V_{cr}-V}$" + "\n")
            
            ax_iv.text(0.33, 0.96, equation_para, transform=ax_iv.transAxes, verticalalignment='top', fontsize=15, bbox=dict(boxstyle="round", fc="white", alpha=0.90, ec="#27AE60"), color="#2C3E50")
            ax_iv.text(0.01, 0.92, equation_latex, transform=ax_iv.transAxes, verticalalignment='top', fontsize=15, bbox=dict(boxstyle="round", fc="white", alpha=0.90, ec="green"), color="black")
                         
        ax_iv.set_ylabel("Current (nA)", fontweight='bold',fontsize=14)
        ax_iv.set_yscale(self.scale_var.get())
        if not show_geiger: ax_iv.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)   
        ax_iv.grid(True, which='both', linestyle='--', alpha=0.5)
        ax_iv.set_title(self.module_name.get(), pad=35) 
        ax_iv.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=4, mode="expand", borderaxespad=0., frameon=False, fontsize=14)
        
        self.canvas_analysis.draw()
        self.plot_notebook.select(self.tab_analysis)
        print(f"Analysis Complete. Vbd: {v_bd_deriv:.2f}V")
        

    def change_scale(self):
        if hasattr(self, 'ax') and self.ax:
            current_scale = self.scale_var.get()
            self.ax.set_yscale(current_scale)
            self.ax.relim()
            self.ax.autoscale_view()
            if self.figure_canvas: self.figure_canvas.draw()

    def set_plot_on_or_off(self, val=1):
        self.plt_flag = val

    def measure_current(self):

        self.instrument.write(":SENS:FUNC 'CURR'")
        self.instrument.write(":CURR:RANG:AUTO ON")
        data = self.instrument.query("READ?")
        vals = data.split(',')
        return float(vals[1])

    def measure_voltage(self):
        self.instrument.write("*CLS")
        self.instrument.write(":OUTP ON")
        self.instrument.write("*WAI")
        data = self.instrument.query("READ?")
        vals = data.split(',')
        return float(vals[0])

    '''def measure_all(self):
        vol=0
        curr=0
	for i from 0 to 5:
        	self.instrument.write("*CLS")
        	self.instrument.write(":OUTP ON")
        	self.instrument.write("*WAI")
        	data =(self.instrument.query("READ?"))
        	
        
        	vals = data.split(',')
        	vol=vol+float(vals[0])
        	curr=curr+float(vals[1])
        	
        return vol/5,curr/5'''
    def validate_and_run(self):
        try:
        # Attempt to convert the string to an integer
            iterations = int(self.Nmeas.get())
            print(iterations)
        # Optional: Check if the number is positive
            if iterations <= 0:
                self.Nmeas.set("5") 
                raise ValueError("Number must be greater than zero, reseting it to 5")
                      
                 # If successful, continue with your logic
            print(f"Starting with {iterations} measurements.")

        except ValueError:
        # This triggers if int() fails or if our custom 'positive' check fails
            self.Nmeas.set("5")
            msg.showwarning("Input Error", "Please enter a valid whole number for 'No. Meas Per Step'.,  reseting it to 5")   

    def measure_all(self):
       total_vol = 0.0
       total_curr = 0.0
       iterations = int(self.Nmeas.get())

    # 1. Clear status and turn Output ON once
       self.instrument.write("*CLS")
       self.instrument.write(":OUTP ON")
       self.instrument.write("*WAI") 
    
       try:
           for i in range(iterations):
               data = self.instrument.query("READ?")
               vals = data.split(',')
               total_vol += float(vals[0])
               total_curr += float(vals[1])
  
       except Exception as e:
           print(f"Error during measurement: {e}")
       return total_vol / iterations, total_curr / iterations        

    def check_output_state(self):
        output_state = self.instrument.query(":OUTPUT:STATE?")
        if output_state.strip() == "1": return 1
        else: return 0

    def set_current_threshold(self, threshold):
        self.instrument.write(":SOUR:VOLT:ILIM %f" % threshold)

    def clr_n_reset_powersupply(self, vol_step):
        self.instrument.write(":OUTP ON")
        self.instrument.write("*WAI")
        voltage_r1 = self.measure_voltage()
        while voltage_r1 > 1e-10:
            voltage_r1 = voltage_r1 - vol_step
            voltage_r2, current_r2 = self.setVoltage(voltage_r1)
            self.instrument.write("*WAI")
            voltage_r1 = voltage_r2
        voltage_r2, current_r2 = self.setVoltage(0.0)
        self.instrument.write("*CLS")
        self.instrument.write("*RST")

    def ramp_down_zero(self, v_step=1.0, delay_t=0.01):
        self.instrument.write("*WAI")
        voltage_r1 = self.measure_voltage()
        curr_r1 = self.measure_current()
        print('Threshod crossed Curr:: ', curr_r1*1e9, ' VOlTAGE:: ', voltage_r1)
        self.rmp_dwn_flag = 1
        end_volt = 0
        diff = abs(self.end_vol - voltage_r1)
        if (diff <= v_step): v_step = v_step / 2.0
        self.ramp_up(end_volt, v_step, delay_t)

    def chk_polarity(self, voltage, pol_voltage):
        if voltage > pol_voltage: return 1
        else: return 0

    def ramp_up_run(self, voltage_r1, voltage, vol_step, polar1, sec_t):
        if abs(voltage_r1 - voltage) > 1e-2:
            voltage_r1 = voltage_r1 + vol_step
            polar2 = self.chk_polarity(voltage_r1, voltage)
            if polar1 != polar2:
                voltage_r2, current_r2 = self.setVoltage(voltage)
                current_r2 = current_r2 * 1000000000.0
                if self.var.get() == 1:
                    temp, humid = self.run_arduino()
                    self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 4)) + ' nA' + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                    self.labels1.config(text=self.p_reading.get())
                else:
                    self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 4)) + ' nA')
                    self.labels1.config(text=self.p_reading.get())
                if (self.rmp_dwn_flag == 1):
                    self.instrument.write("OUTP OFF")
                    self.rmp_dwn_flag = 0
                return
            voltage_r2, current_r2 = self.setVoltage(voltage_r1)
            current_r2 = current_r2 * 1000000000.0
            time.sleep(sec_t)
            if self.var.get() == 1:
                temp, humid = self.run_arduino()
                self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 2)) + ' nA' + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                self.labels1.config(text=self.p_reading.get())
            else:
                self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 2)) + ' nA')
                self.labels1.config(text=self.p_reading.get())
            self.instrument.write("*WAI")
            if (abs(voltage_r2 - voltage) <= 1e-2):
                voltage_r2, current_r2 = self.setVoltage(voltage)
                current_r2 = current_r2 * 1000000000.0
                if self.var.get() == 1:
                    temp, humid = self.run_arduino()
                    self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 2)) + ' nA' + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                    self.labels1.config(text=self.p_reading.get())
                else:
                    self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 2)) + ' nA')
                    self.labels1.config(text=self.p_reading.get())
                if (self.rmp_dwn_flag == 1):
                    self.instrument.write("OUTP OFF")
                    self.rmp_dwn_flag = 0
                    return
                return
            if (self.rmp_dwn_flag == 1): print('Ramping Down\nCurr:: ', current_r2, ' VOlTAGE:: ', voltage_r2)
            self.window.after(1000, lambda: self.ramp_up_run(voltage_r1, voltage, vol_step, polar1, sec_t)) # Tanay Dey
        else:
            voltage_r2, current_r2 = self.setVoltage(voltage)
            current_r2 = current_r2 * 1000000000.0
            if self.var.get() == 1:
                temp, humid = self.run_arduino()
                self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 4)) + ' nA' + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                self.labels1.config(text=self.p_reading.get())
            else:
                self.p_reading.set('VOLTAGE:: ' + str(round(voltage_r2, 4)) + ' V\n' + 'CURRENT::' + str(round(current_r2, 4)) + ' nA')
                self.labels1.config(text=self.p_reading.get())
            if (self.rmp_dwn_flag == 1):
                self.instrument.write("OUTP OFF")
                self.rmp_dwn_flag = 0
            return

    def ramp_up(self, voltage, vol_step=.50, sec_t=0.01):
        self.instrument.write("*WAI")
        voltage_r1 = self.measure_voltage()
        current_r1 = self.measure_current() * 1000000000.0
        if self.var.get() == 1:
            temp, humid = self.run_arduino()
            self.p_reading.set('VOLTAGE:: ' + str(round(self.measure_voltage(), 4)) + ' V\n' + 'CURRENT::' + str(round(current_r1, 2)) + ' nA' + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
            self.labels1.config(text=self.p_reading.get())
        else:
            self.p_reading.set('VOLTAGE:: ' + str(round(self.measure_voltage(), 4)) + ' V\n' + 'CURRENT::' + str(round(current_r1, 2)) + ' nA')
            self.labels1.config(text=self.p_reading.get())

        indx = 0
        polar1 = self.chk_polarity(voltage_r1, voltage)
        if voltage_r1 > voltage: vol_step = -1.0 * vol_step
        if abs(voltage) < 0.5 and abs(voltage_r1) < 1:
            self.setVoltage(voltage)
            if (self.rmp_dwn_flag == 1):
                self.instrument.write("OUTP OFF")
                self.rmp_dwn_flag = 0
        else:
            self.ramp_up_run(voltage_r1, voltage, round(vol_step, 1), polar1, sec_t)

    def set_single_voltage(self):
        if (self.is_blank_string(self.p_address.get()) == True or self.search_flag == 0):
            msg.showwarning('warning', 'Power supply is not detected \n SEARCH OR SET SOURCE ADDRESS')
            return 0
        step_voltage_num = 0
        delay_time_num = 1
        if (self.down_step_voltage.get() == '' or self.delay_time.get() == ''):
            step_voltage_num = 5
            delay_time_num = 1
        else:
            flag4, step_voltage_num = self.is_number(self.step_voltage.get())
            flag5, delay_time_num = self.is_number(self.delay_time.get())
        self.ramp_up(float(self.single_voltage.get()), step_voltage_num, delay_time_num)

    def ramp_down_single_voltage(self):
        if (self.is_blank_string(self.p_address.get()) == True or self.search_flag == 0):
            msg.showwarning('warning', 'Power supply is not detected \n SEARCH OR SET SOURCE ADDRESS')
            return 0
        step_voltage_num = 0
        delay_time_num = 1
        if (self.down_step_voltage.get() == '' or self.delay_time.get() == ''):
            step_voltage_num = 5
            delay_time_num = 1
        else:
            flag4, step_voltage_num = self.is_number(self.down_step_voltage.get())
            flag5, delay_time_num = self.is_number(self.delay_time.get())
        self.ramp_down_zero(step_voltage_num, delay_time_num)

    def setVoltage(self, voltage):
        try:
            self.instrument.write("SOUR:FUNC VOLT")
            self.instrument.write("*WAI")
            self.instrument.write("SOUR:VOLT:RANG:AUTO ON")
            self.instrument.write("SOUR:VOLT %f" % voltage)
            self.instrument.write("*WAI")
            self.instrument.write("OUTP ON")
            self.instrument.write("*WAI")
            voltage_r = self.measure_voltage()
            current = self.measure_current()
            return voltage_r, current
        except visa.VisaIOError as error:
            print("Error:", error)
            return -9999, -9999

    def set_output_off(self):
        self.instrument.write("OUTP OFF")

    '''def find_powersupply(self, location):
        flag = 0
        address_powersupply = ''
        for loc in location:
            index = loc.find('USB')
            if (index != -1):
                flag = 1
                address_powersupply = loc
                break
        return flag, address_powersupply'''
    def manual_port_selection(self, port_list):
        selected_port = {"value": None}

        top = Tk.Toplevel(self.window)
        top.title("Select Power Supply Port")
        top.geometry("600x400")

        Tk.Label(
            top,
            text="Auto-detection failed.\nSelect port manually:",
            font=("Arial", 11)
        ).pack(pady=10)

        listbox = Tk.Listbox(top, width=50)
        listbox.pack(pady=10, fill=Tk.BOTH, expand=True)

        for port in port_list:
            listbox.insert(Tk.END, port)

        def confirm_selection():
            try:
                selected_port["value"] = listbox.get(Tk.ACTIVE)
            except:
                selected_port["value"] = None
            top.destroy()

        Tk.Button(top, text="Select", command=confirm_selection).pack(pady=10)

        # Wait until window closes
        self.window.wait_window(top)

        return selected_port["value"]        

    def find_powersupply(self, location):
        flag = 0
        address_powersupply = ''

        # Automatic detection
        for loc in location:
            if 'USB' in loc:
                flag = 1
                address_powersupply = loc
                return flag, address_powersupply

        # Manual selection if not found
        selected = self.manual_port_selection(location)

        if selected:
            flag = 1
            address_powersupply = selected

        return flag, address_powersupply


    def find_powersupply1(self, location):
        flag = 0
        index = location.find('USB')
        if (index != -1): flag = 1
        return flag, location

    def search(self):
        self.rm = visa.ResourceManager()
        self.plot_VI_graph(-1, 1)
        try:
            self.location = self.rm.list_resources()
            print(self.location)
            value, address = self.find_powersupply(self.location)
            if (value == 1):
                self.instrument = self.rm.open_resource(address)
                print('Power supply Detected at address:: ', address)
                self.p_address.set(address)
                self.window.after(0, self.show_yellow_light)
                self.search_flag = 1
            else:
                msg.showwarning("warning", "Power supply is not detected.")
                self.p_address.set('')
                self.window.after(0, self.show_red_light)
        except pyvisa.Error as e:
            self.window.after(0, self.show_red_light)
            self.p_address.set('')
            msg.showerror("Error", f"PyVISA error: {e}")

    def is_blank_string(self, s):
        return not s or s.isspace()

    def set_address(self):
        self.rm = visa.ResourceManager()
        try:
            if (self.is_blank_string(self.p_address.get()) == False):
                flag, address = self.find_powersupply1(self.p_address.get())
                if flag == 1:
                    self.instrument = self.rm.open_resource(address)
                    msg.showinfo("Information", "Power supply address set::\n" + self.p_address.get())
                    self.p_address.set(address)
                    self.window.after(0, self.show_yellow_light)
                    self.search_flag = 1
                    return
                else:
                    msg.showwarning("warning", "Set a valid address")
                    self.window.after(0, self.show_red_light)
                    self.p_address.set('')
                    return
            else:
                msg.showwarning("warning", "First Search or Set Power supply address")
                self.p_address.set('')
                return
        except visa.Error as e:
            self.window.after(0, self.show_red_light)
            self.p_address.set('')
            msg.showerror("Error", f"PyVISA error: {e}")

    def search_or_set(self):
        self.plot_VI_graph(0, 0)
        self.sim_flag = 0
        if (self.is_blank_string(self.p_address.get()) == False): self.set_address()
        else: self.search()

    def get_temp_dir(self):
        temp_dir = tempfile.gettempdir()
        return os.path.join(temp_dir, '_MEI<some_random_string>')

    def light_images(self, s1):
        if getattr(sys, 'frozen', False): base_path = sys._MEIPASS
        else: base_path = os.path.abspath(".")
        image_path = os.path.join(base_path, 'light_files', str(s1))
        try:
            image = Image.open(image_path)
            resized_image1 = image.resize((200, 85))
            photo1 = ImageTk.PhotoImage(resized_image1)
            return photo1
        except Exception: return None

    def exits(self, event=None):
        self.window.quit()

    def HVTEST(self):
        user_response = msg.askquestion("Positive IV TEST", "Positive IV test selected \n Do you want to continue?").lower()
        if user_response in ('no', 'n'):
            self.user_answer.set('')
            return
        else:
            self.start_voltage.set('0')
            self.end_voltage.set('30')
            self.current_th.set('10')
            self.step_voltage.set('0.5')
            self.down_step_voltage.set('5')
            self.delay_time.set('0.5')

    def IVTEST(self):
        user_response = msg.askquestion("Negative IV TEST", "Negative IV test selected \n Do you want to continue?").lower()
        if user_response in ('no', 'n'):
            self.user_answer.set('')
            return
        else:
            self.start_voltage.set('0')
            self.end_voltage.set('30')
            self.current_th.set('10')
            self.step_voltage.set('0.5')
            self.down_step_voltage.set('5')
            self.delay_time.set('0.5')
            return

    def is_number(self, num):
        try: return True, float(num)
        except ValueError: return False

    def RUN_IV_HV(self):
        self.run_time_flag = 0
        if (self.is_blank_string(self.p_address.get()) == True or self.search_flag == 0):
            msg.showwarning('warning', 'Power supply is not detected \n SEARCH OR SET SOURCE ADDRESS')
            return 0
        if (self.user_answer.get() == ''):
            msg.showwarning('warning', 'Please choose any option from \nTEST TYPE HV/IV')
            return 0

        try:
            flag1, current_th_num = self.is_number(self.current_th.get())
            flag2, start_voltage_num = self.is_number(self.start_voltage.get())
            flag3, end_voltage_num = self.is_number(self.end_voltage.get())
            flag4, step_voltage_num = self.is_number(self.step_voltage.get())
            flag5, delay_time_num = self.is_number(self.delay_time.get())
            
            if (not flag1 or not flag2 or not flag3 or not flag4 or not flag5):
                 msg.showwarning('warning', 'Please provide numbers to the parameters')
                 return 0
            
            self.run_time_flag = 1
            return 1
        except Exception as e:
            msg.showerror("Error", f"{e}")
            return 0

    def auto_run_process(self):
        temp, humid = 0.0, 0.0
        temp1, humid1 = '', ''
        if self.var.get() == 1:
            temp1, humid1 = self.run_arduino()
            if temp1 == '-999': self.stop_run()
            elif temp1 == '-998':
                while temp1 == '-998':
                    temp1, humid1 = self.run_arduino()
                    time.sleep(0.3)
            temp = float(temp1)
            humid = float(humid1)
        try:
            if abs(self.end_vol - self.start_vol) > 1e-3:
                polarrun = self.chk_polarity(self.end_vol, self.start_vol)
                if self.polarinit != polarrun: self.start_vol = self.end_vol
                if abs(self.start_vol) > 1: self.ramp_up(self.start_vol, self.step_vol, self.time_delay)
                else: self.setVoltage(self.start_vol)
                time.sleep(self.time_delay)
                voltage_tmp = self.measure_voltage()
                current_tmp = self.measure_current()
                current_tmp_store = current_tmp * 1000000000.0
                current_tmp_1 = abs(current_tmp)
                diff_I = abs(self.curr_th - current_tmp_1) * 1000000000.0
                diff_V = abs(self.start_vol - voltage_tmp)
                #if abs(voltage_tmp<0.0001 and current_tmp<0:
                #	self.window.after(100, self.auto_run_process)                                
                if (diff_V >= self.step_vol or diff_I <= 20) and self.warn_flag == 0:
                    current_tmp1 = current_tmp * 1000000000.0
                    warning_message = 'WARNING: Limit reached \n ' + str(round(current_tmp1, 1)) + ' nA'
                    print(warning_message)
                    self.warn_flag = 1
                    self.xp.append(voltage_tmp)
                    self.xp_ap.append(self.start_vol)
                    self.yp.append(current_tmp_store)
                    self.temp_arr.append(float(temp))
                    self.humid_arr.append(float(humid))
                    self.time_arr.append(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
                    self.plot1.set_data(self.xp, self.yp)
                    self.plot2.set_data(self.xp_ap, self.yp)
                    self.ax.relim()
                    self.ax.autoscale_view()
                    self.figure_canvas.draw()
                    self.ramp_down_zero(self.down_step_vol, self.time_delay)
                    self.save_results()
                    if self.calc_vbd_var.get(): self.window.after(100, self.run_breakdown_analysis)
                    self.window.after(0, self.show_yellow_light)
                    return

                if (abs(voltage_tmp - self.end_vol) < 1e-3 or abs(self.start_vol - self.end_vol) < 1e-3):
                    self.run_flag = 0
                    self.ramp_down_zero(self.down_step_vol, self.time_delay)
                    self.save_results()
                    if self.calc_vbd_var.get(): self.window.after(100, self.run_breakdown_analysis)
                    self.window.after(0, self.show_yellow_light)
                    return
                if abs(voltage_tmp)>0.0001 and current_tmp>0:
                 voll_avg,curr_avg=self.measure_all()
                 #self.xp.append(voltage_tmp)
                 self.xp.append(voll_avg) #changed tanay
                 self.xp_ap.append(self.start_vol)
                 #self.yp.append(current_tmp * 1e9)
                 self.yp.append(curr_avg * 1e9) #changed by tanay
                 self.temp_arr.append(float(temp))
                 self.humid_arr.append(float(humid))
                 self.time_arr.append(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
                 self.plot1.set_data(self.xp, self.yp)
                 self.plot2.set_data(self.xp_ap, self.yp)
                 self.plot5.set_data(self.xp, self.temp_arr)
                 self.plot6.set_data(self.xp, self.humid_arr)
                 self.ax.relim()
                 self.ax.autoscale_view()
                 self.figure_canvas.draw()
                
                if self.check_output_state() == 1:
                    if self.pause_plot == 0 and self.stop_flag == 0:
                        self.window.after(100, self.auto_run_process)
                    elif self.stop_flag == 1:
                        self.ramp_down_zero(self.down_step_vol, self.time_delay)
                        if self.calc_vbd_var.get(): self.window.after(100, self.run_breakdown_analysis)
                        self.save_results()    
                        self.window.after(0, self.show_yellow_light)
                        return
                    else:
                        self.window.after(0, self.show_yellow_light)
                        return

                if True:
                    if self.start_vol <= self.end_vol: self.start_vol = self.start_vol + self.step_vol
                    else: self.start_vol = self.start_vol - self.step_vol

            else:
                self.ramp_down_zero(self.down_step_vol, self.time_delay)
                if self.calc_vbd_var.get(): self.window.after(100, self.run_breakdown_analysis)
                self.save_results()    
                self.window.after(1000, self.show_yellow_light)
                return

        except Exception as e:
            msg.showerror("Error", f"{e}")

    def start_process(self, event=None):
        self.sim_flag = 0
        self.run_flag = self.RUN_IV_HV()
        if (not self.run_flag): return
        self.run_index = 0
        self.warn_flag = 0
        self.legn_flag = 0
        self.run_flag = 1
        self.pause_plot = 0
        self.stop_flag = 0
        self.run_init_flg = 0
        self.xp = []
        self.yp = []
        self.ypp = []
        self.xp_ap = []
        self.temp_arr = []
        self.humid_arr = []
        self.time_arr = []

        flag1, current_th_num = self.is_number(self.current_th.get())
        flag2, start_voltage_num = self.is_number(self.start_voltage.get())
        flag3, end_voltage_num = self.is_number(self.end_voltage.get())
        flag4, step_voltage_num = self.is_number(self.step_voltage.get())
        flag6, down_step_vol_num = self.is_number(self.down_step_voltage.get())
        flag5, delay_time_num = self.is_number(self.delay_time.get())

        self.start_vol = start_voltage_num
        self.end_vol = end_voltage_num
        self.step_vol = step_voltage_num
        self.down_step_vol = down_step_vol_num
        self.time_delay = delay_time_num

        self.plot_VI_graph(-1, 1)
        if (self.user_answer.get() == 'HV'): self.ax.set_ylim(0.001, current_th_num * 1e3 + 10)
        
        self.polarinit = self.chk_polarity(self.end_vol, self.start_vol)
        current_th_num = current_th_num * 0.000001
        self.curr_th = current_th_num
        self.set_current_threshold(current_th_num)
        self.show_green_light()
        if self.var.get()==1:
           self.ax2.set_ylim(0, 80)
        
        self.auto_run_process()

    def sensel_current(self,indx):
        return self.current_array_sim[indx]*1e-9
        
    def simulation_run(self, event=None):
        self.module_name.set(f"IV Characteristic of a SenSL SiPM")
        self.sim_flag = 1
        self.plot_VI_graph(-1, 1)
        self.pause_plot = 0
        self.warn_flag = 0
        self.stop_flag = 0
        if self.var.get()==1:
           self.ax2.set_ylim(0, 80)

        self.xp = []
        self.yp = []
        self.ypp = []
        self.xp_ap = []
        self.temp_arr = []
        self.humid_arr = []
        self.time_arr = []
        self.run_index = 0
        self.simulation()
            
    def simulation(self):
        self.window.after(0, self.show_green_light)
        temp = '25'
        humid = '30'
        if self.var.get() == 1:
            temp, humid = self.run_arduino()
            if temp == '-999': self.stop_run()

        voltage = self.voltage_array_sim[self.run_index]
        self.xp.append(voltage)
        self.temp_arr.append(float(temp))
        self.humid_arr.append(float(humid))
        self.time_arr.append(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        cur = self.sensel_current(self.run_index)
        self.yp.append(cur)
        self.plot1.set_data(self.xp, self.yp)
        self.plot5.set_data(self.xp, self.temp_arr)
        self.plot6.set_data(self.xp, self.humid_arr)

        if self.var.get() == 1:
            self.p_reading.set('VOLTAGE:: ' + str(round(voltage, 4)) + ' V\n' + 'CURRENT::' + str(round(cur, 8)) + ' μA' + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
            self.labels1.config(text=self.p_reading.get())
        else:
            self.p_reading.set('VOLTAGE:: ' + str(round(voltage, 4)) + ' V\n' + 'CURRENT::' + str(round(cur, 8)) + ' μA')
            self.labels1.config(text=self.p_reading.get())

        self.ax.relim()
        self.ax.autoscale_view()
        self.figure_canvas.draw()
        time.sleep(self.time_delay)
        self.run_index = self.run_index + 1
        
        if self.run_index < 60 and self.warn_flag == 0:
            if self.pause_plot == 0 and self.stop_flag == 0: self.window.after(100, self.simulation)
            else: self.window.after(0, self.show_yellow_light)
        else:
            self.sim_flag = 0
            self.window.after(0, self.show_red_light)
            if self.calc_vbd_var.get(): self.window.after(100, self.run_breakdown_analysis)
            self.save_results()
            
    def pause_plots(self, event=None):
        if (self.stop_flag == 0):
            if self.pause_plot == 0:
                self.pause_plot = 1
                self.pause.config(text='RESUME')
            else:
                self.pause_plot = 0
                self.pause.config(text='PAUSE')
                if self.sim_flag == 1: self.simulation()
                elif self.run_flag == 1: self.auto_run_process()
        else:
            msg.showwarning("warning", 'Run is stopped.Can\'t resume. Please start again.')

    def stop_run(self, event=None):
        should_save = (self.run_flag == 1) or (self.sim_flag == 1)
        self.pause_plot = 1
        self.stop_flag = 1
        self.pause.config(text='PAUSE')
        if self.run_flag == 1:
            try:
                flag4, step_volt = self.is_number(self.step_voltage.get())
                flag5, delay_t = self.is_number(self.delay_time.get())
                if self.instrument: self.instrument.write("*WAI")
                self.ramp_down_zero(step_volt, delay_t)
            except Exception as e: print(f"Ramp Down Error: {e}")
        self.run_flag = 0
        self.sim_flag = 0
        if should_save:
            self.window.after(500, self.save_results)
            if self.calc_vbd_var.get(): self.window.after(600, self.run_breakdown_analysis)

    def get_sub(self, x):
        normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
        sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
        res = x.maketrans(''.join(normal), ''.join(sub_s))
        return x.translate(res)

    def get_super(self, x):
        normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
        super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
        res = x.maketrans(''.join(normal), ''.join(super_s))
        return x.translate(res)

    def multicolor_ylabel(self, axs, list_of_strings, list_of_colors, axis='y', anchorpad=0, xx=0.0, yy=0.0, **kw):
        if axis == 'x' or axis == 'both':
            boxes = [TextArea(text, textprops=dict(color=color, ha='left', va='bottom', **kw)) for text, color in zip(list_of_strings, list_of_colors)]
            xbox = HPacker(children=boxes, align="center", pad=0, sep=5)
            anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad, frameon=False, bbox_to_anchor=(0.2, -0.09), bbox_transform=axs.transAxes, borderpad=0.)
            axs.add_artist(anchored_xbox)
        if axis == 'y' or axis == 'both':
            boxes = [TextArea(text, textprops=dict(color=color, ha='left', va='bottom', rotation=90, **kw)) for text, color in zip(list_of_strings[::-1], list_of_colors)]
            ybox = VPacker(children=boxes, align="center", pad=0, sep=5)
            anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(xx, yy), bbox_transform=axs.transAxes, borderpad=0.)
            axs.add_artist(anchored_ybox)

    def plot_VI_graph(self, voltage_start, voltage_end):
        if hasattr(self, 'keithley_img_frame'): self.keithley_img_frame.pack_forget()
        if hasattr(self, 'analysis_artists'): self.analysis_artists = [] 

        self.figure.clf() 
        self.figure.subplots_adjust(left=0.12, right=0.90, top=0.92, bottom=0.12)

        self.ax = self.figure.add_subplot(111)
        self.ax2 = self.ax.twinx()

        self.plot1, = self.ax.plot([], [], 'o-', color='#3498DB', markersize=4, label="Measured I-V")
        self.plot2, = self.ax.plot([], [], 'x', color='#E74C3C', markersize=4, label="Set I-V")
        self.plot3, = self.ax.plot([], [], 'b', linestyle='None', label="Limit")
        self.plot4, = self.ax2.plot([], [], 'ro', linestyle='None', label="Set voltage")
        self.plot5, = self.ax2.plot([], [], 'bd', label="Temp")
        self.plot6, = self.ax2.plot([], [], 'ms', label="Humidity")

        if (self.sim_flag == 0):
            if self.var.get() == 1:
                self.ax2.set_visible(True)
                color_map = {'TEMP': 'red', 'Humidity': 'green', 'C': 'blue'}
                label = 'TEMP in {}C'.format(self.get_super('o'))
                self.multicolor_ylabel(self.ax2, (label, 'AND', 'Humidity in %'), ('m', 'k', 'b'), axis='y', size=15, xx=1.05, yy=0.2, weight='bold')
            else:
                self.ax2.set_visible(False)
            self.multicolor_ylabel(self.ax, ('Current', 'in nA '), ('r', 'r'), axis='y', size=15, weight='bold', xx=-0.08)
        else:
            if self.var.get() == 0:
                self.ax2.set_visible(False)
            else:
                self.ax2.set_visible(True)
                color_map = {'TEMP': 'red', 'Humidity': 'green', 'C': 'blue'}
                label = 'TEMP in {}C'.format(self.get_super('o'))
                self.multicolor_ylabel(self.ax2, (label, 'AND', 'Humidity in %'), ('m', 'k', 'b'), axis='y', size=15, xx=1.05, yy=0.2, weight='bold')
                
            #self.multicolor_ylabel(self.ax, ('Current', 'in μA'), ('r', 'r'), axis='y', size=15, weight='bold', xx=-0.08)
            self.multicolor_ylabel(self.ax, ('Current', 'in μA'), ('#3498DB', '#3498DB'), axis='y', size=15, weight='bold', xx=-0.08, yy=0.5)
            #self.ax2.set_ylabel('Current in μA', color='b')

        self.ax.set_xlabel('Voltage in V', color='green')
        self.ax.set_title(self.module_name.get())
        self.ax.tick_params(colors='#3498DB', axis='y')
        self.ax.tick_params(colors='green', axis='x')
        self.ax2.tick_params(colors='blue', axis='y')
        
        self.ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5, color='gray')
        self.ax.set_facecolor('white')
        self.ax2.grid(False) 

        current_scale = self.scale_var.get()
        self.ax.set_yscale(current_scale)

        self.ax.set_ylim(auto=True)
        self.ax2.set_ylim(auto=True)
        self.ax2.set_xlim(auto=True)
        self.ax.set_xlim(auto=True)
        
        if self.figure_canvas: self.figure_canvas.get_tk_widget().pack_forget()
        
        self.figure_canvas = FigureCanvasTkAgg(self.figure, master=self.tab_measure)
        self.figure_canvas.get_tk_widget().pack(anchor="center", fill=Tk.BOTH, expand=True) 
        self.figure_canvas.draw()
        
        self.plot_notebook.select(self.tab_measure)

    def search_all_words(self, my_string, words):
        for word in words:
            if word in my_string: return True
        return False

    def init_arduino(self):
        ports = serial.tools.list_ports.comports()
        self.all_ports = []
        words_to_find = ["ACM", "VID", "PID", "SER", "LOCATION"]#,"USB"]
        for port in ports:
            if self.search_all_words(port.device, words_to_find) or self.search_all_words(port.description, words_to_find):
                self.all_ports.append(port.device)
        if len(self.all_ports) < 1:
            msg.showwarning("warning", "Arduino is not found")
            self.var.set(0)
            self.stop_run()
        elif self.var.get() == 1:
            self.arduino_ports.set(self.all_ports[0])
            if self.ard_flag == 0:
                self.ser = serial.Serial(self.arduino_ports.get(), self.baud_rate)

    def arduino_port_on_select(self, event):
        self.arduino_ports.set(self.arduino_port_list.get())
        print("Selected:", self.arduino_ports.get())

    def check_button_clicked(self, var):
        if hasattr(self, 'label8') and self.label8: self.label8.destroy()
        if hasattr(self, 'arduino_port_list') and self.arduino_port_list: self.arduino_port_list.destroy()

        selected = var.get()
        if selected:
            self.all_ports = []
            self.init_arduino()
            if not self.all_ports: self.all_ports = ["No Device Found"]
            self.label8 = Label(self.env_group, text='SELECT ARDUINO PORT', fg='red', font=("arial", 9, "bold"), bg=self.colors['bg_sidebar'])
            self.label8.pack(anchor='w', padx=5, pady=(5, 0))
            self.arduino_port_list = ttk.Combobox(self.env_group, values=self.all_ports)
            if self.all_ports: self.arduino_port_list.set(self.all_ports[0])
            self.arduino_port_list.pack(fill=Tk.X, padx=5, pady=(0, 5))
            self.arduino_port_list.bind("<<ComboboxSelected>>", self.arduino_port_on_select)
        else:
            self.arduino_ports.set('')            

    def run_arduino(self):
        try:
            self.ser.write(b"all\n")
            l1 = self.ser.readline().decode('utf-8', errors='ignore').strip()
            temp = '-998'
            humid = '-998'
            numbers = re.findall(r'\d+\.\d+', str(l1))
            if len(numbers) >= 2: temp, humid = numbers[0], numbers[1]
            return str(temp), str(humid)
        except Exception:
            return -999, -999

    def save_results(self):
        user_response = msg.askquestion("Save results", "Do You Want to save results?").lower()
        if user_response in ('no', 'n'):
            alldata = pd.DataFrame({"VOLTS": self.xp, "CURRNT_NAMP": self.yp, "TEMP_DEGC": self.temp_arr, "RH_PRCNT": self.humid_arr, "TIME": self.time_arr})
            alldata.to_csv('temp.csv', index=False) 
            return
        
        outfile = self.module_name.get()
        outfile = re.sub(r'\s+', '', outfile)
        outfile = "".join(outfile.split())
        outfile = outfile.replace(":", "")
        self.current_datetimes.set(datetime.now().strftime("%d-%m-%Y-%H-%M"))
        directory = './Results/' + str(self.current_datetimes.get()) + '_' + outfile

        if os.path.exists(directory) == True:
            user_response = msg.askquestion("Path Clashes", "Same Module \nDo You Want to Continue?").lower()
            if user_response in ('yes', 'y'):
                while os.path.exists(directory) == True:
                    directory = directory + '_clone'
                os.makedirs(directory)
            else: return
        else:
            os.makedirs(directory)
            
        outfile = directory + '/' + str(self.current_datetimes.get()) + '_' + outfile + '_Result'
        log_file = outfile + '_Log.csv'
        
        # Save CSV
        alldata = pd.DataFrame({"VOLTS": self.xp, "CURRNT_NAMP": self.yp, "TEMP_DEGC": self.temp_arr, "RH_PRCNT": self.humid_arr, "TIME": self.time_arr})
        alldata.to_csv(log_file, index=False)
        
        # Save Measurement Graph
        meas_plot = outfile + '_IV_Graph.png'
        self.figure.savefig(meas_plot)
        
        # Save Analysis Graph (if available)
        try:
            analysis_plot = outfile + '_Analysis_Graph.png'
            self.fig_analysis.savefig(analysis_plot)
        except Exception:
            pass

    def show_yellow_light(self):
        self.image_label2.config(image=self.photo2)

    def show_red_light(self):
        self.image_label2.config(image=self.photo1)

    def show_green_light(self):
        self.image_label2.config(image=self.photo3)

     # ==========================================
    # POST PROCESS TAB METHODS
    # Add these methods to your KeithleyGUI class
    # ==========================================

    def _setup_post_process_tab(self):
        """Setup the Post Process tab"""
        # Get screen dimensions
        try:
            import screeninfo
            screen = screeninfo.get_monitors()[0]
            width = screen.width
            height = screen.height
        except:
            width = 1920
            height = 1080
        
        # Colors
        COLORS = {
            "bg_left": "#2C3E50",
            "bg_right": "#ECF0F1",
            "panel": "#34495E",
            "text": "#070707",
            "muted": "#95A5A6",
            "accent": "#3498DB",
            "success": "#2ECC71",
            "danger": "#E74C3C",
            "warning": "#F39C12",
            "header": "#1ABC9C"
        }
        
        # Left Panel
       # Left Panel (Scrollable)
        post_left = Frame(self.tab3, width=220, bg=COLORS["bg_left"])
        post_left.pack(side=Tk.LEFT, fill=Tk.Y)
        post_left.pack_propagate(False)

        post_left_canvas = Tk.Canvas(
            post_left,
            bg=COLORS["bg_left"],
            highlightthickness=0,
            width=220
        )
        post_left_scrollbar = ttk.Scrollbar(
            post_left,
            orient="vertical",
            command=post_left_canvas.yview
        )

        post_left_canvas.configure(yscrollcommand=post_left_scrollbar.set)

        post_left_scrollbar.pack(side=Tk.RIGHT, fill=Tk.Y)
        post_left_canvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True)

        post_left_inner = Frame(post_left_canvas, bg=COLORS["bg_left"])
        canvas_window = post_left_canvas.create_window(
            (0, 0),
            window=post_left_inner,
            anchor="nw"
        )

        def _resize_canvas(event):
            post_left_canvas.itemconfig(canvas_window, width=event.width)

        post_left_canvas.bind("<Configure>", _resize_canvas)

        post_left_inner.bind(
            "<Configure>",
            lambda e: post_left_canvas.configure(
                scrollregion=post_left_canvas.bbox("all")
            )
        )

        def _on_mousewheel(event):
            post_left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        post_left_canvas.bind_all("<MouseWheel>", _on_mousewheel)


        
        # Right Panel
        post_right = Frame(self.tab3, bg=COLORS["bg_right"])
        post_right.pack(side=Tk.RIGHT, fill=Tk.BOTH, expand=True)
        
        # Header
        header_frame = Frame(post_left_inner, bg=COLORS["header"], height=30)
        header_frame.pack(side=Tk.TOP, fill=Tk.X)
        header_frame.pack_propagate(False)
        
        Label(header_frame, text="Controls Panel", 
            bg=COLORS["header"], fg="white", 
            font=("Arial", 10, "bold")).pack(pady=6)
        
        # Main controls frame
        Monitor_frame3 = Frame(post_left_inner)
        Monitor_frame3.pack(side=Tk.TOP, anchor="nw", fill=Tk.X)
        Monitor_frame3.grid_rowconfigure(0, weight=1)
        
        # Section: File Selection
        file_section = Frame(Monitor_frame3, bg=COLORS["panel"], bd=1, relief="groove")
        file_section.pack(fill=Tk.X, pady=(0, 4), padx=1)
        
        Label(file_section, text="File", bg=COLORS["panel"], fg=COLORS["header"], 
            font=("Arial", 12, "bold")).pack(anchor="w", padx=4, pady=(2, 1))
        
        file_entry = Tk.Entry(file_section, textvariable=self.selected_log_file, 
                            font=('Arial', 12), bg="#ECF0F1", relief="flat")
        file_entry.pack(fill=Tk.X, pady=(0, 2), padx=4)
        
        select_file_btn = Button(file_section, text="Select CSV", 
                                command=self.select_log_file,
                                bg=COLORS["accent"], fg="white", 
                                font=("Arial", 12, "bold"),
                                relief="flat", cursor="hand2",
                                activebackground="#2980B9",
                                height=1, padx=2, pady=2)
        select_file_btn.pack(fill=Tk.X, pady=(0, 2), padx=4)
        
        # Section: Voltage Range
        voltage_section = Frame(Monitor_frame3, bg=COLORS["panel"], bd=1, relief="groove")
        voltage_section.pack(fill=Tk.X, pady=(0, 4), padx=1)
        
        Label(voltage_section, text="Voltage", bg=COLORS["panel"], fg=COLORS["warning"], 
            font=("Arial", 10, "bold")).pack(anchor="w", padx=4, pady=(2, 1))
        
        # Voltage Start
        Label(voltage_section, text="Start:", bg=COLORS["panel"], 
            fg=COLORS["text"], font=("Arial", 10,"bold")).pack(anchor="w", padx=4, pady=(1, 0))
        
        v_start_frame = Frame(voltage_section, bg=COLORS["panel"])
        v_start_frame.pack(fill=Tk.X, padx=4, pady=(0, 2))
        
        self.voltage_start_slider = Scale(v_start_frame, from_=self.voltage_min, to=self.voltage_max, 
                                    resolution=0.1, orient=Tk.HORIZONTAL,
                                    variable=self.x_start_var, bg=COLORS["panel"], 
                                    fg=COLORS["text"], highlightthickness=0,
                                    troughcolor="#B226DD", activebackground=COLORS["accent"],
                                    command=self.update_voltage_range_from_sliders, 
                                    length=80, width=6, font=("Arial", 10,"bold"))
        self.voltage_start_slider.pack(side=Tk.LEFT, fill=Tk.X, expand=True)
        
        Label(v_start_frame, textvariable=self.x_start_var, bg=COLORS["panel"], 
            fg=COLORS["warning"], font=("Arial", 10, "bold"), width=4).pack(side=Tk.RIGHT, padx=(1, 0))
        
        # Voltage End
        Label(voltage_section, text="End:", bg=COLORS["panel"], 
            fg=COLORS["text"], font=("Arial", 10, "bold")).pack(anchor="w", padx=4, pady=(1, 0))
        
        v_end_frame = Frame(voltage_section, bg=COLORS["panel"])
        v_end_frame.pack(fill=Tk.X, padx=4, pady=(0, 2))
        
        self.voltage_end_slider = Scale(v_end_frame, from_=self.voltage_min, to=self.voltage_max, 
                                resolution=0.1, orient=Tk.HORIZONTAL,
                                variable=self.x_end_var, bg=COLORS["panel"], 
                                fg=COLORS["text"], highlightthickness=0,
                                troughcolor="#B226DD", activebackground=COLORS["accent"],
                                command=self.update_voltage_range_from_sliders, 
                                length=80, width=6, font=("Arial", 10,"bold"))
        self.voltage_end_slider.pack(side=Tk.LEFT, fill=Tk.X, expand=True)
        
        Label(v_end_frame, textvariable=self.x_end_var, bg=COLORS["panel"], 
            fg=COLORS["warning"], font=("Arial", 10, "bold"), width=4).pack(side=Tk.RIGHT, padx=(1, 0))
        
        # Voltage Range Text Box
        Label(voltage_section, text="Range:", bg=COLORS["panel"],
            fg=COLORS["muted"], font=("Arial", 10, 'bold')).pack(anchor="w", padx=4, pady=(1, 0))
        
        self.voltage_range_text = Tk.Text(
            voltage_section,
            height=1,
            wrap="none",
            font=("Arial", 10, 'bold'),
            bg="#ECF0F1",
            relief="flat"
        )
        self.voltage_range_text.pack(fill=Tk.X, padx=4, pady=(0, 2))
        
        Label(voltage_section, text="Title box:", bg=COLORS["panel"],
            fg=COLORS["muted"], font=("Arial", 10, 'bold')).pack(anchor="w", padx=4, pady=(1, 0))
        
        title_box = Tk.Entry(voltage_section, textvariable=self.set_title, 
                            font=('Arial', 10), bg="#ECF0F1", relief="flat")
        title_box.pack(fill=Tk.X, pady=(0, 2), padx=4)
        
        Button(voltage_section, text="Update Plot Range and Title", 
            command=self.apply_voltage_range_from_text,
            bg=COLORS["accent"], fg="white", font=("Arial", 10, "bold"),
            relief="flat", cursor="hand2", activebackground="#2980B9",
            height=1, padx=2, pady=2).pack(fill=Tk.X, pady=(0, 2), padx=4)
        
        # Section: Display Options
        display_section = Frame(Monitor_frame3, bg=COLORS["panel"], bd=1, relief="groove")
        display_section.pack(fill=Tk.X, pady=(0, 4), padx=1)
        
        Label(display_section, text="Display", bg=COLORS["panel"], fg=COLORS["success"], 
            font=("Arial",12, "bold")).pack(anchor="w", padx=4, pady=(2, 1))
        
        Checkbutton(display_section, text="Log Scale", variable=self.log_scale_var, 
                    bg=COLORS["panel"], fg=COLORS["text"], selectcolor="#34495E",
                    font=("Arial", 12, "bold"), activebackground=COLORS["accent"],
                    command=lambda: self.post_plot(self.selected_log_file.get(), 
                                                self.x_start_var.get(), 
                                                self.x_end_var.get())).pack(anchor="w", padx=4, pady=1)
        
        Checkbutton(display_section, text="Temp & Humidity", variable=self.show_temp_hum_var, 
                    bg=COLORS["panel"], fg=COLORS["text"], selectcolor="#34495E",
                    font=("Arial", 12, "bold"), activebackground=COLORS["panel"],
                    command=lambda: self.post_plot(self.selected_log_file.get(), 
                                                self.x_start_var.get(), 
                                                self.x_end_var.get())).pack(anchor="w", padx=4, pady=1)
        
        Checkbutton(display_section, text="Breakdown Voltage(V)", variable=self.breakdown_voltage_var, 
                    bg=COLORS["panel"], fg=COLORS["text"], selectcolor="#34495E",
                    font=("Arial", 12, "bold"), activebackground=COLORS["panel"],
                    command=lambda: self.post_plot(self.selected_log_file.get(), 
                                                self.x_start_var.get(), 
                                                self.x_end_var.get())).pack(anchor="w", padx=4, pady=1)

        Checkbutton(display_section, text="Giger Probability", variable=self.giger_prob_var, 
                    bg=COLORS["panel"], fg=COLORS["text"], selectcolor="#34495E",
                    font=("Arial", 12, "bold"), activebackground=COLORS["panel"],
                    command=lambda: self.post_plot(self.selected_log_file.get(), 
                                                self.x_start_var.get(), 
                                                self.x_end_var.get())).pack(anchor="w", padx=4, pady=1)
        Checkbutton(display_section, text="DCR", variable=self.show_dcr_var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'],command=self.open_dcr_window ).pack(anchor='w')
        Label(display_section, text="Curr at OVV", bg=COLORS["panel"], fg=COLORS["success"], 
            font=("Arial",12, "bold")).pack(anchor="w", padx=4, pady=(2, 1))
        curr_ov = Tk.Entry(display_section, textvariable=self.set_ovv, 
                            font=('Arial', 10), bg="#ECF0F1", relief="flat")
        curr_ov.pack(fill=Tk.X, pady=(0, 2), padx=4)
        
        Label(display_section, text="Unit:", bg=COLORS["panel"], 
            fg=COLORS["text"], font=("Arial",12,"bold")).pack(anchor="w", padx=4, pady=(1, 0))
        
        current_unit_box = ttk.Combobox(display_section, textvariable=self.current_unit_var,
                                        values=["A", "mA", "µA", "nA"], state="readonly", 
                                        width=6, font=("Arial",12,"bold"))
        current_unit_box.pack(fill=Tk.X, pady=(0, 2), padx=4)
        current_unit_box.bind("<<ComboboxSelected>>", lambda e: self.post_plot(
            self.selected_log_file.get(), self.x_start_var.get(), self.x_end_var.get()))
        
        # Section: Actions
        action_section = Frame(Monitor_frame3, bg=COLORS["panel"], bd=1, relief="groove")
        action_section.pack(fill=Tk.X, pady=(0, 4), padx=1)
        
        Label(action_section, text="Actions", bg=COLORS["panel"], fg=COLORS["accent"], 
            font=("Arial", 12, "bold")).pack(anchor="w", padx=4, pady=(2, 1))
        
        Button(action_section, text="Save", bg=COLORS["danger"], fg="white",
            font=("Arial", 12, "bold"), relief="flat", cursor="hand2",
            activebackground="#C0392B", command=self.save_plot,
            height=1, padx=2, pady=2).pack(fill=Tk.X, pady=(0, 2), padx=4)
        
        Button(action_section, text="Live Data", bg=COLORS["success"], fg="white",
            font=("Arial", 12, "bold"), relief="flat", cursor="hand2",
            activebackground="#27AE60", command=lambda: self.live_data(),
            height=1, padx=2, pady=2).pack(fill=Tk.X, pady=(0, 2), padx=4)
        Button(action_section, text="EXIT", bg=self.colors['warning'], fg='white', font=('Segoe UI', 9), relief=Tk.FLAT, pady=0, command=self.exits).pack(fill=Tk.X, pady=1)        
        # -------- Placeholder Image --------
        self.image_label_vi = None
        
        if getattr(sys, 'frozen', False):
            base_path = sys._MEIPASS
        else:
            base_path = os.path.abspath(".")
        
        image_path_vi = os.path.join(base_path, "light_files", "keithley.png")
        
        try:
            img = Image.open(image_path_vi)
            img = img.resize((1000, 600))
            self.photo_vi = ImageTk.PhotoImage(img)
            
            self.image_label_vi = Label(
                post_right,
                image=self.photo_vi,
                bg=COLORS["bg_right"]
            )
            self.image_label_vi.place(relx=0.5, rely=0.5, anchor="center")
        
        except Exception as e:
            print("Placeholder image error:", e)
        
        self.plot_frame = Frame(post_right, bg="white")
        self.plot_frame.place(relx=0, rely=0, relwidth=1, relheight=1)
        self.plot_frame.lower()   # keep plot behind image initially


    def post_plot(self, log_file, voltage_start=None, voltage_end=None):
        """Enhanced post-processing plot with breakdown voltage and Geiger probability"""
        
        # ---------------- Screen size ----------------
        try:
            import screeninfo
            screen = screeninfo.get_monitors()[0]
            width, height = screen.width, screen.height
        except:
            width, height = 1920, 1080

        if not log_file:
            self.show_placeholder()
            return

        self.hide_placeholder()

        # Clear previous canvas
        if self.post_canvas:
            self.post_canvas.get_tk_widget().destroy()
            self.post_canvas = None

        # ---------------- Load data ----------------
        try:
            data = pd.read_csv(log_file)
            
            voltage = data['VOLTS'].to_numpy()
            current_nA = data['CURRNT_NAMP'].to_numpy()
            temperature = data['TEMP_DEGC'].to_numpy()
            humidity = data['RH_PRCNT'].to_numpy()
            voltage_fit = voltage
            current_fit = current_nA
        except Exception as e:
            msg.showerror("Error", f"Failed to load CSV file: {e}")
            return

        # GUI unit scaling
        scale_factor = self.CURRENT_SCALE[self.current_unit_var.get()]
        current = (current_nA * 1e-9) * scale_factor

        # ---------------- Voltage range ----------------
        self.voltage_min, self.voltage_max = voltage.min(), voltage.max()

        if voltage_start is not None and voltage_end is not None:
            mask = (voltage >= voltage_start) & (voltage <= voltage_end)
            voltage = voltage[mask]
            current = current[mask]
            current_nA = current_nA[mask]
            temperature = temperature[mask]
            humidity = humidity[mask]
            #voltage_fit=voltage[mask]
            #current_fit = current_nA[mask]
        else:
            self.x_start_var.set(self.voltage_min)
            self.x_end_var.set(self.voltage_max)

        # Check which features to show
        show_bd = self.breakdown_voltage_var.get()
        show_geiger = self.giger_prob_var.get()

        # ---------------- Breakdown voltage analysis ----------------
        v_bd_fit = None
        popt = None
        fit_success = False
        
        if show_bd or show_geiger:
            try:
                # Calculate breakdown voltage
                v_bd_deriv = find_vbd_derivative(voltage, current_nA)
                
                # Fit the model
                popt, fit_success = optimize_fit(voltage, current_nA, v_bd_deriv, 
                                                user_params=self.user_fit_params if hasattr(self, 'user_fit_params') else None)
                
                if fit_success:
                    v_bd_fit = popt[0]
            except Exception as e:
                print(f"Breakdown analysis failed: {e}")
                fit_success = False

        # ---------------- Figure layout ----------------
        if show_geiger and fit_success:
            fig = plt.figure(figsize=(width*0.0078, height*0.008))
            fig.subplots_adjust(top=0.82)
            gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.1)
            ax1 = fig.add_subplot(gs[0])
            ax_prob = fig.add_subplot(gs[1], sharex=ax1)
            plt.setp(ax1.get_xticklabels(), visible=False)
        else:
            fig, ax1 = plt.subplots(figsize=(width*0.0078, height*0.008))
            fig.subplots_adjust(top=0.82)
            ax_prob = None

        # ---------------- Main V–I plot ----------------
        # Plot measured data
        star = mpath.Path.unit_regular_star(6)
        circle = mpath.Path.unit_circle()
        cut_star = mpath.Path(vertices=np.concatenate([circle.vertices, star.vertices[::-1, ...]]), 
                            codes=np.concatenate([circle.codes, star.codes]))
        
        ax1.plot(voltage, current, marker=cut_star, color='indigo', markersize=8, 
                alpha=0.6, label='Measured Data', linestyle='None')

        # Plot breakdown voltage if enabled
        if show_bd and fit_success and v_bd_fit is not None:
            # Plot vertical line at breakdown voltage
            ax1.axvline(
                v_bd_fit,
                color='red',
                linestyle='--',
                linewidth=2,
                alpha=0.7,
                label=rf'$Breakdown\,Voltage \, V_{{bd}} = {v_bd_fit:.2f}\,\mathrm{{V}}$'
            )

            
            # Plot the fitted curve
            v_smooth = np.linspace(min(voltage), min(max(voltage), popt[1]-0.1), 500)
            i_fit_nA = dinu_eq8_model(v_smooth, *popt)
            i_fit_scaled = (i_fit_nA * 1e-9) * scale_factor
            ax1.plot(v_smooth, i_fit_scaled, 'g--', linewidth=2, label='Fit Model')
            
            # Mark breakdown point
            y_val_nA = dinu_eq8_model(v_bd_fit, *popt)
            y_val_scaled = (y_val_nA * 1e-9) * scale_factor
            ax1.plot(v_bd_fit, y_val_scaled, 'rx', markersize=12, 
                    markeredgewidth=3, label='Breakdown Point')


            
            #################################################################################
            overvol=float(self.set_ovv.get())#2.5
            y_val_nA_ov = dinu_eq8_model(v_bd_fit+overvol, *popt)* 1e-9* scale_factor
            text_pos=max(voltage)-3 #v_bd_fit+overvol +5
            ax1.plot(v_bd_fit+overvol, y_val_nA_ov, 'mP', markersize=10, markeredgewidth=2,label=f"Current at {overvol:0.2f} Overvoltage")
            ann=ax1.annotate(f"$V_{{bd}}+overvol$: {v_bd_fit+overvol:.2f} V\n I: {y_val_nA_ov:0.2f} {self.current_unit_var.get()}", xy=(v_bd_fit+overvol, y_val_nA_ov),     xytext=(text_pos, 0.001*y_val_nA_ov), color='m', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), bbox=dict(boxstyle="round", fc="white", alpha=0.7), fontsize=13)
            ann.draggable()
            #################################################################################
            
            
            
            # Add annotation
            idx = (np.abs(voltage - v_bd_fit)).argmin()
            if idx < len(current):
                ax1.annotate(f'V_bd: {v_bd_fit:.2f}V', 
                            xy=(v_bd_fit, y_val_scaled),
                            xytext=(v_bd_fit + (max(voltage)-min(voltage))*0.1, 
                                y_val_scaled*1.5),
                            color='red', fontweight='bold', fontsize=12,
                            arrowprops=dict(arrowstyle='->', color='red', lw=2),
                            bbox=dict(boxstyle="round,pad=0.5", fc="white", 
                                    alpha=0.8, ec="red"))

        # ---------------- Axis labels and title ----------------
        if not show_geiger:
            ax1.set_xlabel('Voltage (V)', fontsize=16, fontweight='bold')
        
        if self.log_scale_var.get():
            ax1.set_yscale('log')
            ax1.set_ylabel(f"Current ({self.current_unit_var.get()})", 
                        fontsize=16, fontweight='bold')
        else:
            ax1.set_ylabel(f"Current ({self.current_unit_var.get()})", 
                        fontsize=16, fontweight='bold')
        
        ax1.text(
            0.5, 1.18,   # ⬅ higher than legend
            self.set_title.get(),
            transform=ax1.transAxes,
            ha="center",
            va="bottom",
            fontsize=16,
            fontweight="bold",
            clip_on=False
        )


        ax1.grid(True, alpha=0.3, linestyle='--')

        # ---------------- Temp / Humidity on secondary axis ----------------
        if self.show_temp_hum_var.get():
            ax3 = ax1.twinx()
            ax3.plot(voltage, temperature, 'rs--', label='Temp (°C)', 
                    markersize=4, linewidth=1.5)
            ax3.plot(voltage, humidity, 'g^-.', label='RH (%)', 
                    markersize=4, linewidth=1.5)
            ax3.set_ylabel('Temperature (°C) / Humidity (%)', 
                        fontsize=14, fontweight='bold')
            ax3.set_ylim(0, max(max(temperature)*1.2, max(humidity)*1.2, 100))
            ax3.legend(loc='upper right', fontsize=10)
        
        # ---------------- Geiger Probability subplot ----------------
        if show_geiger and ax_prob is not None and fit_success:
            p_factor = popt[2]
            v_bd_fit = popt[0]
            
            v_smooth = np.linspace(min(voltage), max(voltage), 500)
            p_geiger = np.zeros_like(v_smooth)
            mask_aval = v_smooth > v_bd_fit
            if np.any(mask_aval):
                p_geiger[mask_aval] = 1 - np.exp(-p_factor * (v_smooth[mask_aval] - v_bd_fit))
            
            ax_prob.plot(v_smooth, p_geiger, 'b-', linewidth=2.5, label='Geiger Probability')
            ax_prob.fill_between(v_smooth, p_geiger, color='blue', alpha=0.15)
            ax_prob.axvline(v_bd_fit, color='red', linestyle='--', alpha=0.5, linewidth=2)
            ax_prob.axhline(1.0, color='gray', linestyle='--', alpha=0.4)
            
            ax_prob.set_ylabel('Geiger Prob.', fontweight='bold', fontsize=14, color='blue')
            ax_prob.set_xlabel('Bias Voltage (V)', fontweight='bold', fontsize=16)
            ax_prob.set_ylim(-0.05, 1.15)
            ax_prob.grid(True, which='both', linestyle='--', alpha=0.3)
            ax_prob.tick_params(axis='y', labelcolor='blue')
            
            # Add formula
            formula_txt = r"$P_{Geiger} = 1 - e^{-p(V - V_{bd})}$"
            greiger=ax_prob.text(0.3, 0.8, formula_txt, transform=ax_prob.transAxes, 
                        fontsize=12, color='darkblue', fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.9, ec="blue"),
                        ha='right', va='top')
            #greiger.draggable()
            if self.show_dcr_var.get()==0:
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}")
            else:
                if abs(self.C_ucell)>0:
                    DCR=popt[3]*1e-9/(self.C_ucell*1e3)
                else: DCR=0   
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}\n"+f"DCR : {DCR:0.3f} kHz")
    
                               
                
            equation_latex = (r"$I_{tot} = I_{leak} + I_{aval}$" + "\n" + r"$I_{leak} = \exp(aV + b)$" + "\n" + r"$I_{aval} = A \cdot \Delta V \cdot (1 - e^{-p \Delta V}) \cdot \frac{V_{cr}-V_{bd}}{V_{cr}-V}$" + "\n")
            
            eqn=ax1.text(0.01, 0.95,equation_latex,transform=ax1.transAxes,va='top',fontsize=13,zorder=5,bbox=dict(boxstyle="round", fc="white", alpha=0.92, ec="green"))

            para=ax1.text(0.35, 0.95, equation_para,transform=ax1.transAxes,va='top',fontsize=13,zorder=4,bbox=dict(boxstyle="round", fc="white", alpha=0.85, ec="#27AE60"))

            #eqn.draggable()
            #para.draggable()
        if not show_geiger: ax1.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)   
        ax1.grid(True, which='both', linestyle='--', alpha=0.5)
        ax1.legend(
            bbox_to_anchor=(0., 1.00, 1., .102),  # ⬅ lower
            loc='lower left',
            ncol=4,
            mode="expand",
            borderaxespad=0.,
            frameon=False,
            fontsize=14
        )

        # ---------------- Embed in Tkinter ----------------
        self.post_canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.post_canvas.draw()
        self.post_canvas.get_tk_widget().pack(fill='both', expand=True)

    def update_voltage_range_from_sliders(self, val):
        """Update plot and text box when sliders are moved"""
        self.voltage_range_text.delete("1.0", Tk.END)
        self.voltage_range_text.insert(Tk.END, f"{self.x_start_var.get():.2f}, {self.x_end_var.get():.2f}")
        self.post_plot(self.selected_log_file.get(), self.x_start_var.get(), self.x_end_var.get())

    def select_log_file(self):
        file_path = filedialog.askopenfilename(
            title="Select Result CSV File",
            initialdir="./Results",
            filetypes=[("CSV Files", "*.csv")]
        )
        
        if file_path:
            self.selected_log_file.set(file_path)
            
            # Read file to get voltage range
            data = pd.read_csv(file_path)
            v_min = float(data['VOLTS'].min())
            v_max = float(data['VOLTS'].max())
            
            # Store global voltage range
            self.voltage_min = v_min
            self.voltage_max = v_max
            
            # Update slider ranges dynamically
            self.voltage_start_slider.config(from_=v_min, to=v_max)
            self.voltage_end_slider.config(from_=v_min, to=v_max)
            
            # Set initial values to full range
            self.x_start_var.set(v_min)
            self.x_end_var.set(v_max)
            
            # Update text box
            self.voltage_range_text.delete("1.0", Tk.END)
            self.voltage_range_text.insert(Tk.END, f"{v_min:.2f}, {v_max:.2f}")
            
            # Plot with full range
            self.post_plot(file_path, v_min, v_max)

    def save_plot(self):
        if self.post_canvas is None:
            print("No plot to save!")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save Plot",
            defaultextension=".png",
            filetypes=[("JPEG Image", "*.jpeg"),("PNG Image", "*.png"), 
                    ("PDF File", "*.pdf"), ("All Files", "*.*")]
        )
        
        if file_path:
            self.post_canvas.figure.savefig(file_path, dpi=300, bbox_inches="tight")
            print("Plot saved to:", file_path)

    def live_data(self):
        file_path = 'temp.csv'
        self.selected_log_file.set(file_path)
        
        # Read file to get voltage range
        data = pd.read_csv(file_path)
        v_min = float(data['VOLTS'].min())
        v_max = float(data['VOLTS'].max())
        
        # Store global voltage range
        self.voltage_min = v_min
        self.voltage_max = v_max
        
        # Update slider ranges dynamically
        self.voltage_start_slider.config(from_=v_min, to=v_max)
        self.voltage_end_slider.config(from_=v_min, to=v_max)
        
        # Set initial values to full range
        self.x_start_var.set(v_min)
        self.x_end_var.set(v_max)
        
        # Update text box
        self.voltage_range_text.delete("1.0", Tk.END)
        self.voltage_range_text.insert(Tk.END, f"{v_min:.2f}, {v_max:.2f}")
        
        # Plot with full range
        self.post_plot(file_path, v_min, v_max)

    def apply_voltage_range_from_text(self, event=None):
        """Parse voltage range from text box and update sliders + plot"""
        try:
            text = self.voltage_range_text.get("1.0", Tk.END).strip()
            text = text.replace(",", " ").replace("FROM", "").replace("TO", "")
            values = [float(v) for v in text.split() if v.replace('.', '', 1).replace('-', '', 1).isdigit()]
            
            if len(values) >= 2:
                v_start, v_end = values[0], values[1]
                
                if v_start < self.voltage_min:
                    v_start = self.voltage_min
                if v_end > self.voltage_max:
                    v_end = self.voltage_max
                if v_start >= v_end:
                    return
                
                self.x_start_var.set(v_start)
                self.x_end_var.set(v_end)
                
                self.post_plot(self.selected_log_file.get(), v_start, v_end)
        
        except Exception as e:
            print("Invalid voltage range text:", e)

    def show_placeholder(self):
        if self.image_label_vi:
            self.image_label_vi.lift()

    def hide_placeholder(self):
        if self.plot_frame:
            self.plot_frame.lift()

if __name__ == "__main__":
    app = KeithleyGUI()
    app.run()
