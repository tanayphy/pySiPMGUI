# ==========================================
# 1. IMPORTS & SYSTEM CONFIGURATION
# ==========================================
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
# 2. CUSTOM EXCEPTIONS
# ==========================================
class BatchStopRequested(Exception):
    """Raised to cleanly unwind out of run_batch_sequence when the user
    hits the main-GUI Batch Stop button, regardless of which nested loop
    (sweep wait, ramp-down wait, inter-channel delay, inter-board delay)
    is currently executing."""
    pass

# ==========================================
# 3. PHYSICS & DATA ANALYSIS HELPERS
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

def optimize_fit(voltage, current, v_bd_guess, user_params=None, current_std=None):
    """Fit the Dinu Eq.8 breakdown model to (voltage, current).

    current_std, if provided, is the per-point standard deviation of the
    current readings (same units/array-length as `current`). It is
    propagated into log-space (sigma_log = std/I, since the fit is done
    on log(I)) and passed to curve_fit as `sigma` with absolute_sigma=True,
    so that noisier points are down-weighted and the returned covariance
    matrix carries genuine physical units. Returns (popt, success, perr)
    where perr is the 1-sigma parameter uncertainty array (zeros if the
    fit failed or no uncertainty could be computed).
    """
    if len(voltage) < 5:
        return np.zeros(6), False, np.zeros(6)
    print("max(d(log (I)/dV :: )",v_bd_guess)
    mask = voltage < v_bd_guess
    V_cut = voltage[mask]
    I_cut = current[mask]

    if len(V_cut) < 2:
        a, b = 0.1, -5.0
    else:
        mid = len(V_cut) // 2
        if mid == 0:
            mid = 1
        V1, V2 = V_cut[mid-1], V_cut[mid]
        I1, I2 = I_cut[mid-1], I_cut[mid]
        a = (np.log(max(I2, 1e-13)) - np.log(max(I1, 1e-13))) / (V2 - V1) if (V2 - V1) != 0 else 0.1
        b = np.log(max(I1, 1e-13)) - a * V1

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

    # Propagate per-point current uncertainty (std-dev, same units as
    # `current`) into log-space: since the fit is performed on log(I),
    # d(log I) = dI / I, so sigma_log = current_std / current. Points
    # with zero/unknown std fall back to a small floor relative to the
    # signal scale so they don't get an artificial infinite weight, and
    # any non-finite values are floored too so a single bad std doesn't
    # break the whole fit.
    sigma_log = None
    if current_std is not None:
        current_std = np.asarray(current_std, dtype=float)
        if current_std.shape == current.shape:
            floor = 0.05 * np.median(np.abs(current_safe))  # 5% relative floor
            if not np.isfinite(floor) or floor <= 0:
                floor = 1e-13
            std_safe = np.where((current_std <= 0) | ~np.isfinite(current_std), floor, current_std)
            sigma_log = std_safe / current_safe
            sigma_log = np.where(~np.isfinite(sigma_log) | (sigma_log <= 0), floor / current_safe, sigma_log)

    try:
        if sigma_log is not None:
            popt, pcov = curve_fit(fit_wrapper, voltage, log_current_data, p0=p0, bounds=bounds,
                                    sigma=sigma_log, absolute_sigma=True, maxfev=10000)
        else:
            popt, pcov = curve_fit(fit_wrapper, voltage, log_current_data, p0=p0, bounds=bounds, maxfev=10000)
        perr = np.sqrt(np.diag(pcov)) if pcov is not None and np.all(np.isfinite(pcov)) else np.zeros(6)
        return popt, True, perr
    except Exception as e:
        print(f"Fit failed: {e}")
        return np.zeros(6), False, np.zeros(6)

# ==========================================
# 4. MATPLOTLIB HELPER CLASSES
# ==========================================
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
# 5. MAIN GUI & CONTROLLER CLASS
# ==========================================
class KeithleyGUI:
    # ------------------------------------------
    # Class-level constants
    # ------------------------------------------
    # Maximum |voltage| allowed for a Forward-bias sweep (V). Forward
    # characterization of a SiPM (quenching resistance) only needs a small
    # excursion past the diode turn-on knee -- going much higher risks
    # over-driving the device, so anything past this is flagged.
    FORWARD_MODE_MAX_V = 2.0

    # ------------------------------------------
    # 5.1 INITIALIZATION & MACRO LOADING
    # ------------------------------------------
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
        self.address_powersupply=''
        # Plotting Lists
        self.xp = []
        self.yp = []
        self.ypp = []
        self.xp_ap = []
        self.temp_arr = []
        self.humid_arr = []
        self.time_arr = []
        self.curr_std_arr = []   # per-point std-dev of current readings (for error bars)

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
]) #SenSL

        # Instrument & Control Flags
        self.C_ucell=0
        self.rm = None
        self.instrument = None
        self.search_flag = 0
        self.run_time_flag = 0
        self.pause_plot = 0
        # Set while an I-V sweep (single or batch) is frozen because the
        # power supply appeared to disconnect mid-sweep. While True, the
        # main PAUSE/RESUME button (labelled RESUME) and the BATCH
        # PAUSE/RESUME button both route through attempt_reconnect_and_resume()
        # instead of just re-arming the sweep loop, and the batch wait-loop
        # holds in place instead of tearing the channel down. Cleared once
        # a reconnect succeeds or the user presses Stop/Batch Stop.
        self.awaiting_reconnect = False
        self.disconnect_resume_v = None
        self.figure_canvas = None
        self.canvas_analysis = None
        self.plot1 = None
        self.plot2 = None
        self.plot3 = None
        self.plot4 = None
        self.plot5 = None
        self.plot6 = None
        self.errbar_container = None   # I-V error-bar artist (Measured I-V), rebuilt on each redraw
        self.warn_flag = 0
        self.batch_mode_active = False
        self.batch_pause_flag = 0   # batch-level pause, survives per-channel resets
        self.batch_stop_flag  = 0   # batch-level stop, survives per-channel resets
        self.legn_flag = 0
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
        self.ramp_down_complete = False
        # Generation token for ramp_up()/ramp_up_run(). Each call to
        # ramp_up() bumps this counter and stamps its own recursive
        # window.after() chain with the new value. If a second ramp is
        # started (e.g. ramp_down_zero() during a Batch Stop) while an
        # earlier ramp's chain is still pending, the earlier chain's next
        # tick sees a stale token and quietly stops instead of continuing
        # to call setVoltage()/toggle rmp_dwn_flag in parallel with the
        # new ramp. Without this, two concurrent chains fight over the
        # same instrument and the same shared rmp_dwn_flag/
        # ramp_down_complete flags -- which is what caused voltage to
        # climb back up after a ramp-down had already started, and
        # ramp_down_complete to fire before the output was actually at 0V.
        self._ramp_gen = 0

        # Add these to __init__
        self.rq_rbias_var = StringVar(value="0")
        self.rq_ncells_var = StringVar(value="1")
        self.rq_unit_var = StringVar(value="auto")
        self.rq_precision_var = StringVar(value="3")
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
        self.auto_yscale_var = Tk.BooleanVar(value=True)
        self.ymin_var = StringVar(value="0.001")
        self.ymax_var = StringVar(value="10010")
        self.run_mode_var = StringVar(value="single")
        self._batch_config = None

        # Error-bar (std-dev of repeated readings per point) display controls
        self.show_errorbars_var = Tk.BooleanVar(value=True)
        self.errorbar_capsize_var = StringVar(value="4")
        self.errorbar_scale_var = StringVar(value="1.0")  # multiplier applied to error-bar yerr (whisker length)

        # Plot-Range Restriction: lets the user only see/record the I-V
        # curve between Vmin and Vmax, while the region below Vmin is
        # skipped through quickly (no per-point averaging, no plotting)
        # via _quick_skip_to_voltage(). Vmax then acts as the effective
        # stop point for the sweep (ramp-down triggers there instead of
        # waiting for the full End Voltage), independent of what End
        # Voltage is set to.
        self.restrict_plot_range_var = Tk.BooleanVar(value=False)
        self.plot_range_vmin_var = StringVar(value="")
        self.plot_range_vmax_var = StringVar(value="")
        self.plot_range_skip_delay_var = StringVar(value="0.1")
        self._restrict_plot_active = False
        self._restrict_plot_vmax = None
        self._restrict_plot_ascending = True

        self.calc_vbd_var = Tk.BooleanVar(value=True)
        self.show_geiger_var = Tk.BooleanVar(value=True)
        self.show_dcr_var = Tk.BooleanVar(value=False)
        self.show_rq_var = Tk.BooleanVar(value=False)
        self.rq_region_display_var = StringVar(value="both")  # "both", "region1", "region2" -- which Rq fit(s) to show
        self.rq_mode_var = StringVar(value="auto")            # "auto" or "manual"
        self.rq_r1_vmin_var = StringVar(value="")            # Manual region 1 V-min
        self.rq_r1_vmax_var = StringVar(value="")            # Manual region 1 V-max
        self.rq_r2_vmin_var = StringVar(value="")            # Manual region 2 V-min
        self.rq_r2_vmax_var = StringVar(value="")
        self.rq_unit_var = StringVar(value="auto")            # "auto", "ohm", "kohm", "mohm" -- Rq/R_total display unit
        self.rq_precision_var = StringVar(value="3")          # number of digits after the decimal point for Rq displays
        self.analysis_mode_var = StringVar(value="")   # "" (none selected), "forward", or "reverse"
        self.user_fit_params = {}

        # Tracks the last-clicked value for every Tk variable that backs a
        # "deselectable" Radiobutton group (see _make_deselectable below),
        # keyed by the Tk variable's internal name (var._name is unique per
        # StringVar/IntVar instance). Lets a second click on an
        # already-selected radio button clear the whole group back to "".
        self._radio_last_value = {}

        # --- Batch indicator state (canvas LED, set up in setup_gui) ---
        self._batch_led_canvas = None
        self._led_circle = None
        self._led_ring = None
        self._led_glow = None

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

        self.load_main_macro("modules/main.mac")
        self.load_batch_macro("modules/batch.mac")
        self.window.protocol("WM_DELETE_WINDOW", self.exits)

    def load_main_macro(self, filepath="main.mac"):
        if not os.path.exists(filepath): return
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'): continue
                    parts = re.split(r'[\s=:]+', line, 1)
                    if len(parts) < 2: continue
                    key, val = parts[0].strip().lower(), parts[1].strip()
                    
                    if key == 'power_supply_address': self.p_address.set(val)
                    elif key == 'module_name': self.module_name.set(val)
                    elif key == 'start_voltage': self.start_voltage.set(val)
                    elif key == 'end_voltage': self.end_voltage.set(val)
                    elif key == 'step_voltage': self.step_voltage.set(val)
                    elif key == 'down_step_voltage': self.down_step_voltage.set(val)
                    elif key == 'delay_time': self.delay_time.set(val)
                    elif key == 'current_limit': self.current_th.set(val)
                    elif key == 'meas_per_step': self.Nmeas.set(val)
                    elif key == 'polarity': self.user_answer.set(val)
                    elif key == 'analysis_mode':
                        self.analysis_mode_var.set(val.lower())
                        self.on_analysis_mode_selected()
                    elif key == 'plot_scale': self.scale_var.set(val.lower())
                    elif key == 'auto_y': self.auto_yscale_var.set(val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'ymin': self.ymin_var.set(val)
                    elif key == 'ymax': self.ymax_var.set(val)
                    elif key == 'show_errorbars': self.show_errorbars_var.set(val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'errorbar_capsize': self.errorbar_capsize_var.set(val)
                    elif key == 'errorbar_scale': self.errorbar_scale_var.set(val)
                    elif key == 'restrict_plot_range': self.restrict_plot_range_var.set(val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'plot_range_vmin': self.plot_range_vmin_var.set(val)
                    elif key == 'plot_range_vmax': self.plot_range_vmax_var.set(val)
                    elif key == 'plot_range_skip_delay': self.plot_range_skip_delay_var.set(val)
                    elif key == 'show_rq': self.show_rq_var.set(val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'calc_vbd': self.calc_vbd_var.set(val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'show_geiger': self.show_geiger_var.set(val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'show_dcr': self.show_dcr_var.set(val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'single_voltage': self.single_voltage.set(val)
                    elif key == 'rq_rbias': self.rq_rbias_var.set(val)
                    elif key == 'rq_ncells': self.rq_ncells_var.set(val)
                    elif key == 'rq_unit': self.rq_unit_var.set(val)
                    elif key == 'rq_precision': self.rq_precision_var.set(val)                    
                    elif key == 'rq_show': self.rq_region_display_var.set(val.lower()) 
                    elif key == 'arduino_enabled': 
                        is_on = val.lower() in ['1', 'true', 'yes', 'on']
                        self.var.set(1 if is_on else 0)
                        self.check_button_clicked(self.var) # Trigger GUI update
                    elif key == 'arduino_port':
                        self.arduino_ports.set(val)
                    elif key == 'run_mode':
                        self.run_mode_var.set(val.lower())
                    # ... rest of your existing logic ...                                       
                    # Inside the loop in load_main_macro
                    print(f"DEBUG: Found key='{key}', val='{val}'")
            print(f"[Macro] Loaded main settings from {filepath}")
        except Exception as e:
            print(f"[Macro Error] Could not load {filepath}: {e}")

    def load_batch_macro(self, filepath="batch.mac"):
        if not os.path.exists(filepath): return
        if not hasattr(self, '_batch_saved_state') or self._batch_saved_state is None:
            self._batch_saved_state = {}
        boards_dict = {}
        sipms_dict = {}
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'): continue
                    parts = re.split(r'[\s=:]+', line, 1)
                    if len(parts) < 2: continue
                    key, val = parts[0].strip().lower(), parts[1].strip()
                    
                    if key == 'ard_port': self._batch_saved_state['ard_port'] = val
                    elif key == 'baud': self._batch_saved_state['baud'] = val
                    elif key == 'inter_delay': self._batch_saved_state['inter_delay'] = val
                    elif key == 'inter_channel_delay': self._batch_saved_state['inter_channel_delay'] = val
                    elif key == 'autosave': self._batch_saved_state['autosave'] = (val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'overlay_fit_on_iv': self._batch_saved_state['overlay_fit_on_iv'] = (val.lower() in ['1', 'true', 'yes', 'on'])
                    elif key == 'grid_rows': self._batch_saved_state['grid_rows'] = val
                    elif key == 'grid_cols': self._batch_saved_state['grid_cols'] = val
                    elif key == 'batch_id': self._batch_saved_state['batch_id'] = val
                    elif key.startswith('board_'):
                        p = key.split('_')
                        if len(p) >= 3:
                            b_id, b_prop = p[1], "_".join(p[2:])
                            if b_id not in boards_dict: boards_dict[b_id] = {}
                            boards_dict[b_id][b_prop] = val
                    elif key.startswith('sipm_'):
                        p = key.split('_')
                        if len(p) >= 3:
                            s_id, s_prop = p[1], "_".join(p[2:])
                            if s_id not in sipms_dict: sipms_dict[s_id] = {}
                            sipms_dict[s_id][s_prop] = val
            if boards_dict:
                board_list = []
                for k in sorted(boards_dict.keys(), key=lambda x: int(x) if x.isdigit() else x):
                    d = boards_dict[k]
                    board_list.append({
                        "label": d.get("label", f"Board_{k}"),
                        "slave_id": int(d.get("slave_id", 8)),
                        "pin": int(d.get("pin", 2)),
                        "enabled": d.get("enabled", "true").lower() in ['1', 'true', 'yes', 'on'],
                        "participate": d.get("participate", "true").lower() in ['1', 'true', 'yes', 'on']
                    })
                self._batch_saved_state['boards'] = board_list
            if sipms_dict:
                sipm_list = []
                for k in sorted(sipms_dict.keys(), key=lambda x: int(x) if x.isdigit() else x):
                    d = sipms_dict[k]
                    sipm_list.append({
                        "label": d.get("label", f"SiPM_{k}"),
                        "board_label": d.get("board", ""),
                        "pin": int(d.get("pin", 2)),
                        "rbias": d.get("rbias", ""),
                        "ncells": d.get("ncells", ""),
                        "enabled": d.get("enabled", "true").lower() in ['1', 'true', 'yes', 'on'],
                        "participate": d.get("participate", "true").lower() in ['1', 'true', 'yes', 'on']
                    })
                self._batch_saved_state['sipms'] = sipm_list
            print(f"[Macro] Loaded batch settings from {filepath}")
        except Exception as e:
            print(f"[Macro Error] Could not load {filepath}: {e}")

    def run(self):
        self.window.mainloop()

    # ------------------------------------------
    # 5.2 GUI CONSTRUCTION
    # ------------------------------------------
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

    def setup_gui(self):
        """Builds the GUI layout ensuring it fits on screen without forced scrolling."""
        # ------------------------------------------------------------------
        # SECTION: STYLE / WINDOW INITIALISATION
        # PURPOSE: Apply ttk styles and set the root window background before
        #          any widget is created.
        # ------------------------------------------------------------------
        self._configure_styles()
        self.window.configure(bg=self.colors['bg_main'])

        # ------------------------------------------------------------------
        # SECTION: MAIN CONTAINER
        # PURPOSE: Root frame that holds the whole GUI and hosts the tab
        #          control (Notebook) below it.
        # ------------------------------------------------------------------
        # 1. Main Container
        self.master_frame = ttk.Frame(self.window, style='Main.TFrame')
        self.master_frame.pack(fill=Tk.BOTH, expand=True)

        # ------------------------------------------------------------------
        # SECTION: TAB CONTROL
        # PURPOSE: Creates the two top-level tabs: 'IV / HV Characterization'
        #          (tab1, live measurement) and 'Post Process' (tab3, offline
        #          CSV review).
        # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # SECTION: CONNECTION GROUP (sidebar)
        # PURPOSE: PSU VISA address entry + Connect / Load-macro buttons.
        # ------------------------------------------------------------------
        conn_group = ttk.LabelFrame(self.button_frame, text="CONNECTION", style='Group.TLabelframe', padding=5)
        conn_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        conn_inner = ttk.Frame(conn_group, style='Sidebar.TFrame')
        conn_inner.pack(fill=Tk.X)
        ttk.Label(conn_inner, text="Address:", style='Sidebar.TLabel').pack(side=Tk.LEFT)
        self.ps_address_screen = ttk.Entry(conn_inner, textvariable=self.p_address, width=15)
        self.ps_address_screen.pack(side=Tk.LEFT, fill=Tk.X, expand=True, padx=5)
        btn_action_inner = ttk.Frame(conn_group, style='Sidebar.TFrame')
        btn_action_inner.pack(fill=Tk.X, pady=(2, 0))
        ttk.Button(btn_action_inner, text="Connect", style='Action.TButton', command=self.search_or_set).pack(side=Tk.LEFT, fill=Tk.X, expand=True, padx=(0, 2))
        ttk.Button(btn_action_inner, text="Load mac", style='Action.TButton', command=lambda: self.load_main_macro("main.mac")).pack(side=Tk.LEFT, fill=Tk.X, expand=True, padx=(2, 0))

        # ------------------------------------------------------------------
        # SECTION: CONFIG GROUP (sidebar)
        # PURPOSE: Module name entry + Positive/Negative IV polarity
        #          selector + Forward/Reverse analysis-mode selector.
        # ------------------------------------------------------------------
        exp_group = ttk.LabelFrame(self.button_frame, text="CONFIG", style='Group.TLabelframe', padding=5)
        exp_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        ttk.Label(exp_group, text="Module Name:", style='Sidebar.TLabel').pack(anchor='w')
        ttk.Entry(exp_group, textvariable=self.module_name).pack(fill=Tk.X, pady=(0, 5))
        mode_frame = ttk.Frame(exp_group, style='Sidebar.TFrame')
        mode_frame.pack(fill=Tk.X)
        Radiobutton(mode_frame, text='Positive IV', variable=self.user_answer, value='HV', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self._make_deselectable(self.user_answer, 'HV', self.HVTEST)).pack(side=Tk.LEFT, expand=True)
        Radiobutton(mode_frame, text='Negative IV', variable=self.user_answer, value='IV', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self._make_deselectable(self.user_answer, 'IV', self.IVTEST)).pack(side=Tk.LEFT, expand=True)

        # --- Forward/Reverse Analysis Mode Selector ---
        # Starts with neither option selected; start_test_dispatch (via
        # RUN_IV_HV) blocks the run with a popup until one is chosen.
        
        mode_sel_frame = ttk.Frame(exp_group, style='Sidebar.TFrame')
        mode_sel_frame.pack(fill=Tk.X, pady=(5, 0))
        ttk.Label(mode_sel_frame, text="Mode:", style='Sidebar.TLabel').pack(side=Tk.LEFT)
        Radiobutton(mode_sel_frame, text='Forward', variable=self.analysis_mode_var, value='forward',
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'],
                    command=self._make_deselectable(self.analysis_mode_var, 'forward', self.on_analysis_mode_selected, on_deselect=self.toggle_analysis_modes)).pack(side=Tk.LEFT, expand=True)
        Radiobutton(mode_sel_frame, text='Reverse', variable=self.analysis_mode_var, value='reverse',
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'],
                    command=self._make_deselectable(self.analysis_mode_var, 'reverse', self.on_analysis_mode_selected, on_deselect=self.toggle_analysis_modes)).pack(side=Tk.LEFT, expand=True)

        # ------------------------------------------------------------------
        # SECTION: PARAMETERS (V) GROUP (sidebar)
        # PURPOSE: Voltage sweep Start/End/Ramp-Up/Ramp-Down/Delay entries,
        #          current-limit entry + Set button, and Meas/Step entry.
        # ------------------------------------------------------------------
        volt_group = ttk.LabelFrame(self.button_frame, text="PARAMETERS (V)", style='Group.TLabelframe', padding=5)
        volt_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        for i in range(3): volt_group.columnconfigure(i, weight=1)
        ttk.Label(volt_group, text="Start:", style='Sidebar.TLabel').grid(row=0, column=0, sticky='e')
        ttk.Entry(volt_group, textvariable=self.start_voltage, width=6).grid(row=0, column=1, sticky='ew', padx=2)
        ttk.Label(volt_group, text="End:", style='Sidebar.TLabel').grid(row=0, column=2, sticky='e')
        end_v_entry = ttk.Entry(volt_group, textvariable=self.end_voltage, width=6)
        end_v_entry.grid(row=0, column=3, sticky='ew', padx=2)
        end_v_entry.bind("<FocusOut>", self.enforce_forward_end_voltage_limit)
        ttk.Label(volt_group, text="Ramp Up:", style='Sidebar.TLabel').grid(row=1, column=0, sticky='e')
        ttk.Entry(volt_group, textvariable=self.step_voltage, width=6).grid(row=1, column=1, sticky='ew', padx=2)
        ttk.Label(volt_group, text="Ramp Dn:", style='Sidebar.TLabel').grid(row=1, column=2, sticky='e')
        ttk.Entry(volt_group, textvariable=self.down_step_voltage, width=6).grid(row=1, column=3, sticky='ew', padx=2)
        ttk.Label(volt_group, text="Delay(s):", style='Sidebar.TLabel').grid(row=2, column=0, sticky='e')
        ttk.Entry(volt_group, textvariable=self.delay_time, width=6).grid(row=2, column=1, sticky='ew', padx=2)
        ttk.Label(volt_group, text="I-Lim(\u00b5A):", style='Sidebar.TLabel').grid(row=2, column=2, sticky='e')
        ttk.Entry(volt_group, textvariable=self.current_th, width=6).grid(row=2, column=3, sticky='ew', padx=2)
        ttk.Button(volt_group, text="Set", style='Action.TButton', width=4,
                   command=self.apply_current_limit_now).grid(row=2, column=4, sticky='ew', padx=(2, 0))
        ttk.Label(volt_group, text="Meas/Step:", style='Sidebar.TLabel').grid(row=3,columnspan=2, column=0, sticky='e')
        e = ttk.Entry(volt_group, textvariable=self.Nmeas, width=6); e.bind("<FocusOut>", lambda _: self.validate_and_run()); e.grid(row=3, column=2,columnspan=2, sticky='ew', padx=2)
        
        # ------------------------------------------------------------------
        # SECTION: PLOT GROUP (sidebar)
        # PURPOSE: Log/Linear scale selector, Auto/Manual Y-range, error-bar
        #          cap/scale controls, and plot-range restriction controls.
        # ------------------------------------------------------------------
        plot_group = ttk.LabelFrame(self.button_frame, text="PLOT", style='Group.TLabelframe', padding=5)
        plot_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        scale_frame = ttk.Frame(plot_group, style='Sidebar.TFrame')
        scale_frame.pack(fill=Tk.X)
        Radiobutton(scale_frame, text='Log', variable=self.scale_var, value='log', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self._make_deselectable(self.scale_var, 'log', self.change_scale, on_deselect=self.change_scale)).pack(side=Tk.LEFT, expand=True)
        Radiobutton(scale_frame, text='Linear', variable=self.scale_var, value='linear', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self._make_deselectable(self.scale_var, 'linear', self.change_scale, on_deselect=self.change_scale)).pack(side=Tk.LEFT, expand=True)

        yrange_frame = ttk.Frame(plot_group, style='Sidebar.TFrame')
        yrange_frame.pack(fill=Tk.X, pady=(2, 0))
        Checkbutton(yrange_frame, text='Auto Y', variable=self.auto_yscale_var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self.change_scale).grid(row=0, column=0, sticky='w')
        ttk.Label(yrange_frame, text="Yn:", style='Sidebar.TLabel').grid(row=0, column=1, sticky='e', padx=(2,1))
        ttk.Entry(yrange_frame, textvariable=self.ymin_var, width=5).grid(row=0, column=2, sticky='ew', padx=1)
        ttk.Label(yrange_frame, text="Yx:", style='Sidebar.TLabel').grid(row=0, column=3, sticky='e', padx=(2,1))
        ttk.Entry(yrange_frame, textvariable=self.ymax_var, width=5).grid(row=0, column=4, sticky='ew', padx=1)
        ttk.Button(yrange_frame, text="Apply", style='Action.TButton', command=self.apply_y_range, width=5).grid(row=0, column=5, sticky='ew', padx=(2,0))

        errbar_frame = ttk.Frame(plot_group, style='Sidebar.TFrame')
        errbar_frame.pack(fill=Tk.X, pady=(2, 0))
        Checkbutton(errbar_frame, text='Err Bars', variable=self.show_errorbars_var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self.refresh_errorbars).grid(row=0, column=0, sticky='w')
        ttk.Label(errbar_frame, text="Cap:", style='Sidebar.TLabel').grid(row=0, column=1, sticky='e', padx=(2,1))
        errbar_capsize_entry = ttk.Entry(errbar_frame, textvariable=self.errorbar_capsize_var, width=4)
        errbar_capsize_entry.grid(row=0, column=2, sticky='ew', padx=1)
        errbar_capsize_entry.bind("<Return>", lambda _: self.refresh_errorbars())
        errbar_capsize_entry.bind("<FocusOut>", lambda _: self.refresh_errorbars())
        ttk.Label(errbar_frame, text="Scl:", style='Sidebar.TLabel').grid(row=0, column=3, sticky='e', padx=(2,1))
        errbar_scale_entry = ttk.Entry(errbar_frame, textvariable=self.errorbar_scale_var, width=4)
        errbar_scale_entry.grid(row=0, column=4, sticky='ew', padx=1)
        errbar_scale_entry.bind("<Return>", lambda _: self.on_errorbar_scale_changed())
        errbar_scale_entry.bind("<FocusOut>", lambda _: self.on_errorbar_scale_changed())

        # --- Plot Range Restriction: quick-skip below Vmin, stop at Vmax ---
        prange_frame = ttk.Frame(plot_group, style='Sidebar.TFrame')
        prange_frame.pack(fill=Tk.X, pady=(2, 0))
        Checkbutton(prange_frame, text='Plot Range', variable=self.restrict_plot_range_var,
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar']).grid(row=0, column=0, sticky='w')
        ttk.Label(prange_frame, text="Vn:", style='Sidebar.TLabel').grid(row=0, column=1, sticky='e', padx=(2,1))
        ttk.Entry(prange_frame, textvariable=self.plot_range_vmin_var, width=5).grid(row=0, column=2, sticky='ew', padx=1)
        ttk.Label(prange_frame, text="Vx:", style='Sidebar.TLabel').grid(row=0, column=3, sticky='e', padx=(2,1))
        ttk.Entry(prange_frame, textvariable=self.plot_range_vmax_var, width=5).grid(row=0, column=4, sticky='ew', padx=1)
        ttk.Label(prange_frame, text="Skip:", style='Sidebar.TLabel').grid(row=0, column=5, sticky='e', padx=(2,1))
        ttk.Entry(prange_frame, textvariable=self.plot_range_skip_delay_var, width=5).grid(row=0, column=6, sticky='ew', padx=1)

        # ------------------------------------------------------------------
        # SECTION: ANALYSIS GROUP (sidebar)
        # PURPOSE: Forward-bias (Quench Resistance) and reverse-bias
        #          (Breakdown Voltage, Geiger, DCR, Fit Params) analysis
        #          toggles and their configuration entry points.
        # ------------------------------------------------------------------
        analysis_group = ttk.LabelFrame(self.button_frame, text="ANALYSIS", style='Group.TLabelframe', padding=5)
        analysis_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)

        # --- Forward Characteristics subsection ---
        fwd_group = ttk.LabelFrame(analysis_group, text="Forward", style='Group.TLabelframe', padding=3)
        fwd_group.pack(fill=Tk.X, pady=(0, 2))

        # Quench Resistance: checkbox + "Config" button on one row
        rq_row = ttk.Frame(fwd_group, style='Sidebar.TFrame')
        rq_row.pack(fill=Tk.X)
        self.chk_rq = Checkbutton(rq_row, text="Quench R", variable=self.show_rq_var,
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'])
        self.chk_rq.pack(side=Tk.LEFT)
        self.btn_rq_config = ttk.Button(rq_row, text="Config", style='Action.TButton',
                                         width=6, command=self._open_rq_config_popup)
        self.btn_rq_config.pack(side=Tk.RIGHT, padx=(2, 0))

        # Mode indicator label (shows current mode compactly)
        self.lbl_rq_mode = ttk.Label(fwd_group, text="Mode: Auto", style='Sidebar.TLabel',
                                      font=('Segoe UI', 8, 'italic'))
        self.lbl_rq_mode.pack(anchor='w', padx=(18, 0))

        # ── Show Region selector (compact, single row) ───────────────────────
        rq_region_frame = ttk.Frame(fwd_group, style='Sidebar.TFrame')
        rq_region_frame.pack(fill=Tk.X, padx=(2, 0))
        ttk.Label(rq_region_frame, text="Show:", style='Sidebar.TLabel').pack(side=Tk.LEFT)
        self.rdo_rq_both = Radiobutton(rq_region_frame, text='Both', variable=self.rq_region_display_var, value='both',
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'],
                    command=self._make_deselectable(self.rq_region_display_var, 'both', self.refresh_rq_analysis_if_visible, on_deselect=self.refresh_rq_analysis_if_visible))
        self.rdo_rq_both.pack(side=Tk.LEFT, padx=(2, 2))
        self.rdo_rq_r1 = Radiobutton(rq_region_frame, text='R1', variable=self.rq_region_display_var, value='region1',
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'],
                    command=self._make_deselectable(self.rq_region_display_var, 'region1', self.refresh_rq_analysis_if_visible, on_deselect=self.refresh_rq_analysis_if_visible))
        self.rdo_rq_r1.pack(side=Tk.LEFT, padx=(0, 2))
        self.rdo_rq_r2 = Radiobutton(rq_region_frame, text='R2', variable=self.rq_region_display_var, value='region2',
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'],
                    command=self._make_deselectable(self.rq_region_display_var, 'region2', self.refresh_rq_analysis_if_visible, on_deselect=self.refresh_rq_analysis_if_visible))
        self.rdo_rq_r2.pack(side=Tk.LEFT)

        # --- Reverse Characteristics subsection ---
        rev_group = ttk.LabelFrame(analysis_group, text="Reverse", style='Group.TLabelframe', padding=3)
        rev_group.pack(fill=Tk.X, pady=(0, 3))
        rev_row = ttk.Frame(rev_group, style='Sidebar.TFrame')
        rev_row.pack(fill=Tk.X)
        self.chk_vbd = Checkbutton(rev_row, text="Bkdn V", variable=self.calc_vbd_var,
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'])
        self.chk_vbd.pack(side=Tk.LEFT, padx=(0, 6))

        self.chk_geiger = Checkbutton(rev_row, text="Geiger", variable=self.show_geiger_var,
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'])
        self.chk_geiger.pack(side=Tk.LEFT, padx=(0, 6))

        self.chk_dcr = Checkbutton(rev_row, text="DCR", variable=self.show_dcr_var,
                    bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar'], command=self.open_dcr_window)
        self.chk_dcr.pack(side=Tk.LEFT, padx=(0, 6))
        self.btn_fit_params = ttk.Button(rev_row, text="Fit Params", style='Action.TButton', width=8, command=self.open_param_window)
        self.btn_fit_params.pack(side=Tk.LEFT, padx=0)

        #self.btn_fit_params = ttk.Button(analysis_group, text="Fit Params", style='Action.TButton', command=self.open_param_window)
        #self.btn_fit_params.pack(fill=Tk.X, pady=(2,0))

        # ------------------------------------------------------------------
        # SECTION: MANUAL GROUP (sidebar)
        # PURPOSE: Manually set a single bias voltage or ramp it back to zero,
        #          outside of a full sweep.
        # ------------------------------------------------------------------
        manual_group = ttk.LabelFrame(self.button_frame, text="MANUAL", style='Group.TLabelframe', padding=5)
        manual_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        man_row = ttk.Frame(manual_group, style='Sidebar.TFrame')
        man_row.pack(fill=Tk.X)
        ttk.Entry(man_row, textvariable=self.single_voltage, width=8).pack(side=Tk.LEFT, fill=Tk.X, expand=True, padx=(0,5))
        ttk.Button(man_row, text="Set V", width=5, command=self.set_single_voltage).pack(side=Tk.LEFT, padx=1)
        ttk.Button(man_row, text="Zero", width=5, command=self.ramp_down_single_voltage).pack(side=Tk.LEFT, padx=1)

        # ------------------------------------------------------------------
        # SECTION: ENV GROUP (sidebar)
        # PURPOSE: Enable/disable the Arduino environmental (temp/humidity)
        #          sensor readout.
        # ------------------------------------------------------------------
        self.env_group = ttk.LabelFrame(self.button_frame, text="ENV", style='Group.TLabelframe', padding=5)
        self.env_group.pack(fill=Tk.X, padx=5, pady=2, expand=True)
        Checkbutton(self.env_group, text="Arduino", variable=self.var, bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=lambda: self.check_button_clicked(self.var)).pack(anchor='w')

        # ------------------------------------------------------------------
        # SECTION: ACTION FRAME (sidebar, bottom)
        # PURPOSE: Single/Batch run-mode selector and the main START TEST
        #          button that dispatches to the selected mode.
        # ------------------------------------------------------------------
        act_frame = ttk.Frame(self.button_frame, style='Sidebar.TFrame')
        act_frame.pack(fill=Tk.X, padx=5, pady=5, side=Tk.BOTTOM)

        # Run-mode selector: decides what the START TEST button below does.
        # Batch Run requires the batch config to have been set up & confirmed
        # via the BATCH MODE window first (self._batch_config).
        mode_frame = ttk.Frame(act_frame, style='Sidebar.TFrame')
        mode_frame.pack(fill=Tk.X, pady=(0, 3))
        ttk.Label(mode_frame, text="Run Mode:", style='Sidebar.TLabel').pack(side=Tk.LEFT, padx=(0, 5))
        Radiobutton(mode_frame, text='Single', variable=self.run_mode_var, value='single', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self._make_deselectable(self.run_mode_var, 'single')).pack(side=Tk.LEFT)
        Radiobutton(mode_frame, text='Batch', variable=self.run_mode_var, value='batch', bg=self.colors['bg_sidebar'], fg='white', selectcolor=self.colors['bg_sidebar'], activebackground=self.colors['bg_sidebar'], command=self._make_deselectable(self.run_mode_var, 'batch')).pack(side=Tk.LEFT)

        Button(act_frame, text="START TEST", bg=self.colors['success'], fg='white', font=('Segoe UI', 10, 'bold'), relief=Tk.FLAT, pady=2, command=self.start_test_dispatch).pack(fill=Tk.X, pady=1)
        ttk.Label(act_frame, text="Designed By: Dr. Tanay Dey", style='Sidebar.TLabel', font=('Segoe UI', 7)).pack(pady=0)

        # =========================================================================
        # TAB 1: MAIN AREA (RIGHT SIDE)
        # =========================================================================
        right_panel = Tk.Frame(self.tab1, bg='white')
        right_panel.pack(side=Tk.RIGHT, fill=Tk.BOTH, expand=True)

        # 1. Monitor Strip
        # ------------------------------------------------------------------
        # SECTION: MONITOR STRIP (main area, top)
        # PURPOSE: Container row holding the Real-time Readings display,
        #          Single Run controls, Batch Control, System Status LED,
        #          and the Exit button.
        # ------------------------------------------------------------------
        monitor_container = Tk.Frame(right_panel, bg='white')
        monitor_container.pack(side=Tk.TOP, fill=Tk.X, padx=10, pady=5)

        # Real-time Readings: no longer stretches to fill leftover space --
        # sized to its content so the Single Run / Batch Control groups can
        # sit beside it in the same row instead of being pushed to a second row.
        # ------------------------------------------------------------------
        # SUBSECTION: REAL-TIME READINGS DISPLAY
        # PURPOSE: Shows the live voltage/current readout text.
        # ------------------------------------------------------------------
        readout_frame = Tk.LabelFrame(monitor_container, text="Real-time Readings", bg='white', fg='#7f8c8d', font=("arial", 10))
        readout_frame.pack(side=Tk.LEFT, fill=Tk.Y, padx=(0, 10))
        self.labels1 = Label(readout_frame, textvariable=self.p_reading, bg='white', fg='#2980B9', font=("Noto Sans", 16, 'bold'), justify=Tk.LEFT, anchor="w")
        self.labels1.pack(fill=Tk.BOTH, padx=10, pady=5)

        # Single Run controls -- bordered to match the Batch Control group beside it.
        # ------------------------------------------------------------------
        # SUBSECTION: SINGLE RUN CONTROLS
        # PURPOSE: Pause / Stop / Simulate buttons for a single (non-batch)
        #          sweep.
        # ------------------------------------------------------------------
        ctrl_frame = Tk.LabelFrame(monitor_container, text="Single Run", bg='white', fg='#7f8c8d', font=("arial", 10))
        ctrl_frame.pack(side=Tk.LEFT, fill=Tk.Y, padx=(0, 10))
        ctrl_container = Tk.Frame(ctrl_frame, bg='white')
        ctrl_container.pack(padx=5, pady=4)
        self.pause = Button(ctrl_container, text='PAUSE', bg='#E0E0E0', relief=GROOVE, command=self.pause_plots)
        self.pause.pack(side=Tk.LEFT, padx=5)
        Button(ctrl_container, text='STOP', bg='#E74C3C', fg='white', relief=GROOVE, font=("arial", 10, "bold"), command=self.stop_run).pack(side=Tk.LEFT, padx=5)
        Button(ctrl_container, text='Simulate', bg='#E0E0E0', relief=GROOVE, command=self.simulation_run).pack(side=Tk.LEFT, padx=5)

        # Batch Control -- same row as Real-time Readings / Single Run now,
        # rather than its own row below. Scoped to the whole batch run
        # rather than a single channel's sweep (see batch_pause_flag /
        # batch_stop_flag in run_batch_sequence).
        batch_ctrl_frame = Tk.LabelFrame(monitor_container, text="Batch Control", bg='white', fg='#7f8c8d', font=("arial", 10))
        batch_ctrl_frame.pack(side=Tk.LEFT, fill=Tk.Y, padx=(0, 10))
        Button(batch_ctrl_frame, text='BATCH CONFIG', bg='#8E44AD', fg='white', font=("arial", 9, "bold"), relief=GROOVE, command=self.open_batch_config).pack(side=Tk.LEFT, padx=5, pady=4)
        self.batch_pause_btn = Button(batch_ctrl_frame, text='BATCH PAUSE', bg='#E0E0E0', relief=GROOVE, command=self.batch_pause_resume)
        self.batch_pause_btn.pack(side=Tk.LEFT, padx=5, pady=4)
        Button(batch_ctrl_frame, text='BATCH STOP', bg='#E74C3C', fg='white', relief=GROOVE, font=("arial", 10, "bold"), command=self.batch_stop).pack(side=Tk.LEFT, padx=5, pady=4)
        self._batch_status_var = StringVar(value="No batch running")

        # ═══════════════════════════════════════════════════════════════════════
        # SYSTEM STATUS INDICATOR  —  standalone, covers both Single & Batch runs
        # ═══════════════════════════════════════════════════════════════════════
        status_panel = Tk.LabelFrame(monitor_container, text="System Status",
                                     bg='white', fg='#7f8c8d', font=("arial", 10))
        status_panel.pack(side=Tk.LEFT, fill=Tk.Y, padx=(8, 4), pady=0)

        # Status text label above the LED
        self._led_status_text = Tk.StringVar(value="No Power Supply")
        self._led_status_label = Tk.Label(
            status_panel, textvariable=self._led_status_text,
            bg='white', fg='#C0392B',
            font=("Segoe UI", 9, "bold"), anchor='center')
        self._led_status_label.pack(pady=(6, 2), padx=10)

        # Canvas for the glowing LED circle
        _LED_SIZE = 64          # outer canvas dimension
        _LED_R1   = 30          # outer glow ring radius (from centre)
        _LED_R2   = 22          # main body radius
        _LED_R3   = 10          # bright inner highlight radius
        _CX = _LED_SIZE // 2
        _CY = _LED_SIZE // 2
        self._batch_led_canvas = Tk.Canvas(
            status_panel, width=_LED_SIZE, height=_LED_SIZE,
            bg='white', highlightthickness=0)
        self._batch_led_canvas.pack(pady=(0, 6))

        # Layer 1 — soft outer glow ring (large, semi-transparent via stipple)
        self._led_ring = self._batch_led_canvas.create_oval(
            _CX - _LED_R1, _CY - _LED_R1,
            _CX + _LED_R1, _CY + _LED_R1,
            fill='#E74C3C', outline='', stipple='gray25')

        # Layer 2 — main solid body
        self._led_circle = self._batch_led_canvas.create_oval(
            _CX - _LED_R2, _CY - _LED_R2,
            _CX + _LED_R2, _CY + _LED_R2,
            fill='#C0392B', outline='#7B241C', width=2)

        # Layer 3 — bright inner specular highlight (top-left)
        self._led_glow = self._batch_led_canvas.create_oval(
            _CX - _LED_R2 + 5, _CY - _LED_R2 + 5,
            _CX - _LED_R2 + 5 + _LED_R3 * 2,
            _CY - _LED_R2 + 5 + _LED_R3 * 2,
            fill='#EC7063', outline='')
        # ═══════════════════════════════════════════════════════════════════════

        # EXIT — placed beside the System Status panel so it's reachable
        # right where the operator is already watching the status LED.
        exit_panel = Tk.LabelFrame(monitor_container, text="", bg='white', fg='#7f8c8d', font=("arial", 10))
        exit_panel.pack(side=Tk.LEFT, fill=Tk.Y, padx=(4, 0), pady=0)
        Button(exit_panel, text="EXIT", bg=self.colors['warning'], fg='white', font=('Segoe UI', 10, 'bold'),
               relief=GROOVE, width=8, command=self.exits).pack(padx=8, pady=8, expand=True)

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
        self.plot1, = self.ax.plot([], [], 'o', color='#3498DB', markersize=4, label="Measured Data")
        self.plot2, = self.ax.plot([], [], 'x', color='#E74C3C', markersize=4, label=None)
        self.plot3, = self.ax.plot([], [], 'b', linestyle='None', label="Limit")
        self.ax2 = self.ax.twinx()
        self.plot4, = self.ax2.plot([], [], 'ro', linestyle='None', label=None)
        self.plot5, = self.ax2.plot([], [], 'b:', label="Temp")
        self.plot6, = self.ax2.plot([], [], 'g-.', label="Humidity")
        self.ax.set_yscale(self.scale_var.get() or 'linear')
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
        
        # Ensure analysis modes are toggled correctly on startup
        self.toggle_analysis_modes()

    def _make_deselectable(self, var, value, extra_command=None, on_deselect=None):
        """Returns a Radiobutton `command` callback that lets the user
        click an already-selected radio button again to clear the whole
        group (sets `var` back to "").

        Tkinter's Radiobutton has no built-in toggle-off behaviour -- by
        the time `command` fires, `var` already holds `value` whether this
        was a fresh selection or a re-click of the same option. So we keep
        a small per-variable memory of "what was selected before this
        click" in self._radio_last_value, keyed by the Tk variable's
        internal name. If the new value matches what was already selected,
        we know this was a re-click and clear the variable; otherwise we
        just record the new value and (optionally) run the radio button's
        normal command.

        `extra_command`, if given, is called only when the option is being
        newly selected (not on a deselecting click), so existing callbacks
        like HVTEST/IVTEST/on_analysis_mode_selected keep their original
        behaviour and simply don't fire when the user is clearing the
        group. `on_deselect`, if given, is called only when the click
        clears the group (e.g. to re-lock dependent controls).
        """
        var_key = str(var)

        def _callback():
            previous = self._radio_last_value.get(var_key)
            if previous == value:
                # Re-click on the already-selected option -> deselect.
                var.set('')
                self._radio_last_value[var_key] = ''
                if on_deselect is not None:
                    on_deselect()
            else:
                self._radio_last_value[var_key] = value
                if extra_command is not None:
                    extra_command()

        return _callback

    def toggle_analysis_modes(self):
        """Disables/Enables analysis checkbuttons based on the selected mode
        (Forward/Reverse). When neither is selected yet (startup default),
        every analysis control is disabled until the user picks one."""
        mode = self.analysis_mode_var.get()
        if mode == "forward":
            self.chk_rq.config(state=Tk.NORMAL)
            self.btn_rq_config.config(state=Tk.NORMAL)
            self.rdo_rq_both.config(state=Tk.NORMAL)
            self.rdo_rq_r1.config(state=Tk.NORMAL)
            self.rdo_rq_r2.config(state=Tk.NORMAL)
            self.chk_vbd.config(state=Tk.DISABLED)
            self.chk_geiger.config(state=Tk.DISABLED)
            self.chk_dcr.config(state=Tk.DISABLED)
            self.btn_fit_params.config(state=Tk.DISABLED)
        elif mode == "reverse":
            self.chk_rq.config(state=Tk.DISABLED)
            self.btn_rq_config.config(state=Tk.DISABLED)
            self.rdo_rq_both.config(state=Tk.DISABLED)
            self.rdo_rq_r1.config(state=Tk.DISABLED)
            self.rdo_rq_r2.config(state=Tk.DISABLED)
            self.chk_vbd.config(state=Tk.NORMAL)
            self.chk_geiger.config(state=Tk.NORMAL)
            self.chk_dcr.config(state=Tk.NORMAL)
            self.btn_fit_params.config(state=Tk.NORMAL)
        else:
            # Neither Forward nor Reverse selected yet -- lock both groups.
            self.chk_rq.config(state=Tk.DISABLED)
            self.btn_rq_config.config(state=Tk.DISABLED)
            self.rdo_rq_both.config(state=Tk.DISABLED)
            self.rdo_rq_r1.config(state=Tk.DISABLED)
            self.rdo_rq_r2.config(state=Tk.DISABLED)
            self.chk_vbd.config(state=Tk.DISABLED)
            self.chk_geiger.config(state=Tk.DISABLED)
            self.chk_dcr.config(state=Tk.DISABLED)
            self.btn_fit_params.config(state=Tk.DISABLED)

    def on_analysis_mode_selected(self):
        """Callback for the Forward/Reverse Mode radio buttons (now living
        in the CONFIG section). Updates which analysis controls are
        enabled, and then makes sure the Positive/Negative IV choice (and
        therefore the End Voltage) matches the newly selected mode:
          - Forward  -> Positive IV fills End = +FORWARD_MODE_MAX_V,
                        Negative IV fills End = -FORWARD_MODE_MAX_V.
          - Reverse  -> Positive IV fills End = +30, Negative IV fills
                        End = -30 (the original high-voltage behaviour).

        If Positive/Negative IV hasn't been chosen yet, the user is asked
        to pick one now (HVTEST/IVTEST below already prompt the +2/-2 or
        +30/-30 confirmation once that radio is clicked).
        """
        self.toggle_analysis_modes()

        mode = self.analysis_mode_var.get()
        if mode == "":
            return

        current_choice = self.user_answer.get()
        if current_choice == '':
            # Nothing picked yet -- ask the user to choose Positive/Negative IV now.
            msg.showinfo(
                "Select IV Polarity",
                f"{'Forward' if mode == 'forward' else 'Reverse'} mode selected.\n"
                "Please choose Positive IV or Negative IV."
            )
            return

        # Positive/Negative IV was already selected -- re-apply it under
        # the newly chosen Forward/Reverse mode so End Voltage tracks
        # correctly (e.g. switching Reverse -> Forward pulls 30V back to 2V).
        if current_choice == 'HV':
            self.HVTEST()
        elif current_choice == 'IV':
            self.IVTEST()

    def check_forward_voltage_limit(self, start_voltage_num, end_voltage_num):
        """Returns True if the run may proceed, False if it was blocked.
        Only enforced when Forward analysis mode is selected: a Forward
        sweep is meant to stay within FORWARD_MODE_MAX_V of 0 V, since it's
        characterizing the diode's forward turn-on / quenching resistance,
        not a high-voltage reverse breakdown sweep.

        If End Voltage exceeds the limit, it is auto-clamped back to
        +-FORWARD_MODE_MAX_V (sign preserved) instead of just blocking the
        run, and the user is warned that this happened.
        """
        if self.analysis_mode_var.get() != "forward":
            return True
        max_abs_v = max(abs(start_voltage_num), abs(end_voltage_num))
        if max_abs_v > self.FORWARD_MODE_MAX_V:
            clamped = self.FORWARD_MODE_MAX_V if end_voltage_num >= 0 else -self.FORWARD_MODE_MAX_V
            self.end_voltage.set(str(clamped))
            msg.showwarning(
                "Forward Mode Voltage Limit",
                f"Forward mode is selected, but the requested voltage "
                f"({max_abs_v:.2f} V) exceeds the \u00b1{self.FORWARD_MODE_MAX_V:.1f} V "
                f"limit for forward-bias sweeps.\n\n"
                f"End Voltage has been reset to {clamped:.1f} V."
            )
            # Value is now corrected in-place -- let the caller re-read
            # end_voltage and proceed rather than forcing a second click.
            return True
        return True

    def enforce_forward_end_voltage_limit(self, event=None):
        """<FocusOut> handler on the End Voltage entry: if Forward mode is
        active and the user types something outside +-FORWARD_MODE_MAX_V,
        clamp it back (sign preserved) and warn, instead of waiting until
        START TEST is pressed."""
        if self.analysis_mode_var.get() != "forward":
            return
        flagE, ev = self.is_number(self.end_voltage.get())
        if not flagE:
            return
        if abs(ev) > self.FORWARD_MODE_MAX_V:
            clamped = self.FORWARD_MODE_MAX_V if ev >= 0 else -self.FORWARD_MODE_MAX_V
            self.end_voltage.set(str(clamped))
            msg.showwarning(
                "Forward Mode Voltage Limit",
                f"Forward mode sweeps must stay within \u00b1{self.FORWARD_MODE_MAX_V:.1f} V.\n"
                f"End Voltage has been reset to {clamped:.1f} V."
            )

    # ----------------------------------------------------------------------
    # SECTION: FIT-PARAMETERS POPUP
    # PURPOSE: Small Toplevel window letting the user override the initial
    #          guesses (V_bd, V_cr, p, A, leak_a, leak_b) fed into the
    #          breakdown-curve fit optimizer (optimize_fit).
    # ----------------------------------------------------------------------
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

    # ----------------------------------------------------------------------
    # SECTION: DCR-PARAMETERS POPUP
    # PURPOSE: Toplevel window collecting the microcell capacitance
    #          (C_ucell) used to convert fitted amplitude into a Dark
    #          Count Rate (DCR) estimate.
    # ----------------------------------------------------------------------
    def open_dcr_window(self):
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

    def show_placeholder(self):
        if self.image_label_vi:
            self.image_label_vi.lift()

    def hide_placeholder(self):
        if self.plot_frame:
            self.plot_frame.lift()

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

    def exits(self, event=None):
        self.window.quit()

    def _set_led(self, body_color, outline_color, glow_color, ring_color,
                 status_text, text_color):
        """Update the System Status LED and text label."""
        if self._batch_led_canvas is None:
            return
        self._batch_led_canvas.itemconfig(
            self._led_ring,   fill=ring_color)
        self._batch_led_canvas.itemconfig(
            self._led_circle, fill=body_color, outline=outline_color)
        self._batch_led_canvas.itemconfig(
            self._led_glow,   fill=glow_color)
        if hasattr(self, '_led_status_text'):
            self._led_status_text.set(status_text)
        if hasattr(self, '_led_status_label'):
            self._led_status_label.config(fg=text_color)

    def show_green_light(self):
        """I-V sweep is actively running → bright green."""
        self._set_led(
            body_color='#1D8348', outline_color='#0B5345',
            glow_color='#58D68D', ring_color='#A9DFBF',
            status_text='Running', text_color='#1D8348')

    def show_yellow_light(self):
        """Sweep paused / on hold / ramp-down → amber."""
        self._set_led(
            body_color='#D4AC0D', outline_color='#9A7D0A',
            glow_color='#F9E79F', ring_color='#FCF3CF',
            status_text='Paused "\n"or Hold', text_color='#9A7D0A')

    def show_red_light(self):
        """No power supply connected → red."""
        self._set_led(
            body_color='#C0392B', outline_color='#7B241C',
            glow_color='#EC7063', ring_color='#FADBD8',
            status_text='No Power Supply', text_color='#C0392B')

    def show_idle_light(self):
        """Instrument connected but sweep not started → blue/cyan idle."""
        self._set_led(
            body_color='#1A5276', outline_color='#154360',
            glow_color='#5DADE2', ring_color='#D6EAF8',
            status_text='Idle / Ready', text_color='#1A5276')

    def show_complete_light(self):
        """Batch or single run finished successfully → teal."""
        self._set_led(
            body_color='#117A65', outline_color='#0E6655',
            glow_color='#76D7C4', ring_color='#D1F2EB',
            status_text='Complete', text_color='#117A65')

    # ------------------------------------------
    # 5.3 HARDWARE INTERFACE (HAL)
    # ------------------------------------------
    def search(self):
        self.rm = visa.ResourceManager()
        
        # FIX: Do not wipe and reset the canvas if we are recovering from a disconnect mid-batch
        if not self.awaiting_reconnect:
            self.plot_VI_graph(-1, 1)
            
        try:
            self.location = self.rm.list_resources()
            print(self.location)
            value, address = self.find_powersupply(self.location)
            if (value == 1):
                self.instrument = self.rm.open_resource(address)
                print('Power supply Detected at address:: ', address)
                self.p_address.set(address)
                self.window.after(0, self.show_idle_light)
                self.search_flag = 1
            else:
                msg.showwarning("warning", "Power supply is not detected.")
                self.p_address.set('')
                self.window.after(0, self.show_red_light)
        except visa.Error as e:
            self.window.after(0, self.show_red_light)
            self.p_address.set('')
            msg.showerror("Error", f"PyVISA error: {e}")

    def measure_voltage(self):
        self.instrument.write("*CLS")
        self.instrument.write(":OUTP ON")
        self.instrument.write("*WAI")
        data = self.instrument.query("READ?")
        vals = data.split(',')
        return float(vals[0])

    def measure_current(self):
        self.instrument.write(":SENS:FUNC 'CURR'")
        self.instrument.write(":CURR:RANG:AUTO ON")
        data = self.instrument.query("READ?")
        vals = data.split(',')
        return float(vals[1])

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

    def set_current_threshold(self, threshold):
        # NOTE: ":SOUR:VOLT:ILIM" is only valid on Keithley 2450/2470-series
        # SMUs. On a 2410 (and the 2400-series generally) that header does
        # not exist at all and the instrument returns SCPI error -113
        # "Undefined header" regardless of source function/state. The
        # 2400-series sets compliance via the *sense* subsystem's current
        # protection level while sourcing voltage.
        self.instrument.write(":SENS:CURR:PROT %f" % threshold)

    def init_arduino(self):
        ports = serial.tools.list_ports.comports()
        self.all_ports = []
        words_to_find = ["ACM", "VID", "PID", "SER", "LOCATION","Arduino", "CH340", "CP210", "FT232", "ttyACM", "USB Serial"]
        #identifiers = ["Arduino", "CH340", "CP210", "FT232", "ttyACM", "USB Serial"]
        address_powersupply=''#self.address_powersupply#self.address_powersupply#self.find_keithley()
        print(self.address_powersupply,'init arduino')
        if (self.address_powersupply!=''):
          address_powersupply= self.extract_port(self.address_powersupply) 
        else:
             # Power supply is missing, trigger the popup
            msg.showerror(
                title="Hardware Error", 
                message="Please connect the power supply first!"
            )
            return  
        #clean_port = visa_string.replace("ASRL", "").replace("::INSTR", "")
        for port in ports:
            if self.search_all_words(port.device, words_to_find) or self.search_all_words(port.description, words_to_find):
                if address_powersupply != port.device:
                    print(address_powersupply,port.device,'arduino check')
                    self.all_ports.append(port.device)
                   
        if len(self.all_ports) < 1:
            msg.showwarning("warning", "Arduino is not found")
            self.var.set(0)
            self.stop_run()
        elif self.var.get() == 1:
            self.arduino_ports.set(self.all_ports[0])
            if self.ard_flag == 0:
                self.ser = serial.Serial(self.arduino_ports.get(), self.baud_rate)

    def _batch_send_arduino(self, ser_obj, command):
        try:
            ser_obj.write((command.strip() + "\n").encode("utf-8"))
            time.sleep(0.15)
            reply = ""
            while ser_obj.in_waiting > 0:
                reply += ser_obj.readline().decode("utf-8",
                                                    errors="replace").strip()
            return reply
        except Exception as e:
            return f"ERROR: {e}"

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
        count=0
        
        for loc in location:
             
             if 'USB' in loc:
                 count=count+1
        print("\n\ncount=",count)         
        if count >1:
          selected = self.manual_port_selection(location)
          if selected:
            flag = 1
            address_powersupply = selected
          self.address_powersupply=address_powersupply
          return flag, address_powersupply
           
                
        for loc in location:
            if 'USB' in loc:
                flag = 1
                address_powersupply = loc
                self.address_powersupply=address_powersupply
                return flag, address_powersupply
                
        #flag,address_powersupply=self.find_keithley()
        # Manual selection if not found
        if (flag==0):
            selected = self.manual_port_selection(location)

            if selected:
                flag = 1
                address_powersupply = selected
            self.address_powersupply=address_powersupply
            return flag, address_powersupply
        else:
           self.address_powersupply=address_powersupply
           return flag,address_powersupply

    def find_powersupply1(self, location):
        flag = 0
        index = location.find('USB')
        if (index != -1): flag = 1
        return flag, location

    def search_or_set(self):
        # FIX: Do not wipe and reset the canvas if we are recovering from a disconnect mid-batch
        if not self.awaiting_reconnect:
            self.plot_VI_graph(0, 0)
            
        self.sim_flag = 0
        if (self.is_blank_string(self.p_address.get()) == False): self.set_address()
        else: self.search()            

    def is_blank_string(self, s):
        return not s or s.isspace()

    def find_keithley(self):
        rm = visa.ResourceManager()
        resources = rm.list_resources()
    
    # 1. Filter specifically for USB serial adapters to avoid phantom motherboard ports
        serial_ports = [res for res in resources if 'ttyUSB' in res or 'ttyACM' in res]
    
    # Fallback just in case standard ASRL mapping is used
        if not serial_ports:
            serial_ports = [res for res in resources if 'ASRL' in res]
        
        for port in serial_ports:
            try:
                print(f"Testing port: {port}...")
            # Open with a short timeout so the script doesn't hang
                inst = rm.open_resource(port, baud_rate=9600, timeout=1000)
                inst.write_termination = '\n'
                inst.read_termination = '\n'
            
            # Send standard identification query
                idn = inst.query('*IDN?')
                print(idn)
            # Check if it's a Keithley
                if 'KEITHLEY' in idn.upper():
                    print(f"SUCCESS! Found Keithley on {port}: {idn}")
                    return 1,port # Return the active connection
                else:
                    inst.close()
                
            except Exception as e:
            # 2. Catch ALL exceptions (like the termios Input/output error)
            # This ensures the loop gracefully skips bad ports instead of crashing.
                print(f"  -> Skipped {port} due to: {type(e).__name__}")
                pass
            
        print("Could not find a responding Keithley on any serial port.")
        return 0 , None

    def set_address(self):
        self.rm = visa.ResourceManager()
        try:
            if (self.is_blank_string(self.p_address.get()) == False):
                flag, address = self.find_powersupply1(self.p_address.get())
                if flag == 1:
                    self.instrument = self.rm.open_resource(address)
                    msg.showinfo("Information", "Power supply address set::\n" + self.p_address.get())
                    self.p_address.set(address)
                    self.window.after(0, self.show_idle_light)
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

    def get_temp_dir(self):
        temp_dir = tempfile.gettempdir()
        return os.path.join(temp_dir, '_MEI<some_random_string>')

    def get_temp_dir(self):
        temp_dir = tempfile.gettempdir()
        return os.path.join(temp_dir, '_MEI<some_random_string>')

    def check_output_state(self):
        output_state = self.instrument.query(":OUTPUT:STATE?")
        if output_state.strip() == "1": return 1
        else: return 0

    def apply_current_limit_now(self):
        """Push the 'Curr Lim (uA)' field to the instrument's compliance
        immediately, independent of starting a sweep. Lets the user set/
        change compliance on the fly (e.g. right after connecting, or
        mid-session via the 'Set V' controls) without having to start an
        I-V run first, since that was previously the only place
        set_current_threshold() got called."""
        flag, current_th_num = self.is_number(self.current_th.get())
        if not flag:
            msg.showwarning('warning', 'Please enter a valid number (in µA) for Curr Lim.')
            return
        if getattr(self, 'instrument', None) is None:
            msg.showwarning('warning', 'Instrument is not connected.')
            return
        try:
            threshold_A = current_th_num * 1e-6
            self.set_current_threshold(threshold_A)
            self.curr_th = threshold_A
            print(f"[Curr Lim] Compliance set to {current_th_num} µA ({threshold_A:.6g} A)")
        except visa.VisaIOError as error:
            msg.showerror("Error", f"Failed to set current compliance: {error}")

    def set_output_off(self):
        self.instrument.write("OUTP OFF")

    def extract_port(self,visa_string):
    # This pattern looks for EITHER a Linux port (/dev/...) OR a Windows port (COM...)
        pattern = r"(/dev/[a-zA-Z0-9_]+|COM\d+)"
    
        match = re.search(pattern, visa_string)
    
        if match:
            return match.group(0) # Returns just the found port
        else:
            return None # Returns None if no port pattern was found

    def search_all_words(self, my_string, words):
        for word in words:
            if word in my_string: return True
        return False

    def batch_init_arduino(self):
        ports = serial.tools.list_ports.comports()
        self.all_ports = []
        words_to_find = ["ACM", "VID", "PID", "SER", "LOCATION","Arduino", "CH340", "CP210", "FT232", "ttyACM", "USB Serial"]
        #identifiers = ["Arduino", "CH340", "CP210", "FT232", "ttyACM", "USB Serial"]
        address_powersupply=''#self.address_powersupply#self.find_keithley()
        print(self.address_powersupply,'batch init arduino')        
        if (self.address_powersupply!=''):
          address_powersupply= self.extract_port(self.address_powersupply) 
        else:
             # Power supply is missing, trigger the popup
            msg.showerror(
                title="Hardware Error", 
                message="Please connect the power supply first!"
            )
            return  
        #clean_port = visa_string.replace("ASRL", "").replace("::INSTR", "")
        for port in ports:
            if self.search_all_words(port.device, words_to_find) or self.search_all_words(port.description, words_to_find):
                if address_powersupply != port.device:
                    print(address_powersupply,port.device,'arduino check')
                    self.all_ports.append(port.device)
                   
        if len(self.all_ports) < 1:
            msg.showwarning("warning", "Arduino is not found")
            return None
        else:
            return  self.all_ports   

    def arduino_port_on_select(self, event):
        selected_port = self.arduino_port_list.get()
        self.arduino_ports.set(selected_port)
        print("Selected:", selected_port)

        # Actually switch the live serial connection to the newly selected port.
        # Previously this only updated the StringVar, so self.ser stayed bound
        # to whatever port init_arduino() opened first (all_ports[0]).
        try:
            if getattr(self, 'ser', None) and self.ser.is_open:
                self.ser.close()
        except Exception as e:
            print("Error closing previous serial port:", e)

        try:
            self.ser = serial.Serial(selected_port, self.baud_rate)
            self.ard_flag = 0
        except Exception as e:
            msg.showerror(
                title="Hardware Error",
                message=f"Could not connect to {selected_port}:\n{e}"
            )
            self.ser = None

    def run_arduino(self):
        try:
            self.ser.write(b"all\n")
            l1 = self.ser.readline().decode('utf-8', errors='ignore').strip()
            temp = '-998'
            humid = '-998'
            numbers = re.findall(r'\d+\.\d+', str(l1))
            print('temp/humid',numbers)
            if len(numbers) >= 2: temp, humid = numbers[0], numbers[1]
            return str(temp), str(humid)
        except Exception:
            return -999, -999

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

    def chk_polarity(self, voltage, pol_voltage):
        if voltage > pol_voltage: return 1
        else: return 0

    def validate_and_run(self):
        try:
            iterations = int(self.Nmeas.get())
            print(iterations)
            if iterations <= 0:
                self.Nmeas.set("5")
                raise ValueError("Number must be greater than zero, reseting it to 5")
            print(f"Starting with {iterations} measurements.")

        except ValueError:
            self.Nmeas.set("5")
            msg.showwarning("Input Error", "Please enter a valid whole number for 'No. Meas Per Step'.,  reseting it to 5")

    def measure_all(self):
       total_vol = 0.0
       total_curr = 0.0
       iterations = int(self.Nmeas.get())
       curr_readings = []

       self.instrument.write("*CLS")
       self.instrument.write(":OUTP ON")
       self.instrument.write("*WAI")

       try:
           for i in range(iterations):
               data = self.instrument.query("READ?")
               vals = data.split(',')
               total_vol += float(vals[0])
               curr_val = float(vals[1])
               total_curr += curr_val
               curr_readings.append(curr_val)

       except Exception as e:
           print(f"Error during measurement: {e}")

       curr_std = float(np.std(curr_readings)) if len(curr_readings) > 1 else 0.0
       return total_vol / iterations, total_curr / iterations, curr_std

    def set_plot_on_or_off(self, val=1):
        self.plt_flag = val

    def is_number(self, num):
        try: return True, float(num)
        except ValueError: return False

    # ------------------------------------------
    # 5.4 MEASUREMENT & STATE MACHINE LOOP
    # ------------------------------------------
    def start_process(self, event=None):
        self.sim_flag = 0
        self.run_flag = self.RUN_IV_HV()
        if (not self.run_flag): return
        self.run_index = 0
        self.warn_flag = 0
        try:
            if hasattr(self, '_warn_annot') and self._warn_annot is not None:
                self._warn_annot.remove()
                self._warn_annot = None
        except Exception:
            pass
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
        self.curr_std_arr = []

        flag1, current_th_num = self.is_number(self.current_th.get())
        flag2, start_voltage_num = self.is_number(self.start_voltage.get())
        flag3, end_voltage_num = self.is_number(self.end_voltage.get())
        flag4, step_voltage_num = self.is_number(self.step_voltage.get())
        flag6, down_step_vol_num = self.is_number(self.down_step_voltage.get())
        flag5, delay_time_num = self.is_number(self.delay_time.get())

        self.start_vol = start_voltage_num
        self.end_vol = end_voltage_num
        self.step_vol = step_voltage_num
        if not self.batch_mode_active:
            self.down_step_vol = down_step_vol_num
        self.time_delay = delay_time_num

        if not self.batch_mode_active:
            self.plot_VI_graph(-1, 1)

        #if (self.user_answer.get() == 'HV'): self.ax.set_ylim(0.001, current_th_num * 1e3 + 10)

        self.polarinit = self.chk_polarity(self.end_vol, self.start_vol)
        current_th_num = current_th_num * 0.000001
        self.curr_th = current_th_num
        self.set_current_threshold(current_th_num)
        self.show_green_light()
        if self.var.get()==1:
           self.ax2.set_ylim(0, 80)

        # --- Plot Range Restriction (Vmin -> Vmax) ---
        # If enabled, quickly skip (no recording/plotting) from Start
        # Voltage up to just below Vmin, then let the normal per-step
        # sweep take over from there. Vmax becomes the effective stop
        # point for recording -- auto_run_process() checks
        # self._restrict_plot_vmax the same way it checks End Voltage.
        self._restrict_plot_active = False
        self._restrict_plot_vmax = None
        self._restrict_plot_ascending = self.end_vol >= self.start_vol
        if self.restrict_plot_range_var.get():
            flagA, vmin_num = self.is_number(self.plot_range_vmin_var.get())
            flagB, vmax_num = self.is_number(self.plot_range_vmax_var.get())
            if flagA and flagB and abs(vmax_num - vmin_num) > 1e-9:
                lo, hi = min(vmin_num, vmax_num), max(vmin_num, vmax_num)
                ascending = self._restrict_plot_ascending
                skip_target = (lo - self.step_vol) if ascending else (hi + self.step_vol)
                # Never skip past the configured sweep bounds.
                if ascending:
                    skip_target = max(self.start_vol, min(skip_target, self.end_vol))
                else:
                    skip_target = min(self.start_vol, max(skip_target, self.end_vol))
                self._restrict_plot_active = True
                self._restrict_plot_vmax = hi if ascending else lo
                flagC, skip_delay_num = self.is_number(self.plot_range_skip_delay_var.get())
                if not flagC or skip_delay_num < 0:
                    skip_delay_num = 0.1
                if hasattr(self, '_batch_status_var') and self.batch_mode_active:
                    self._batch_status_var.set("Quick-skipping to Plot Range start\u2026")
                if abs(skip_target - self.start_vol) > 1e-6:
                    self._quick_skip_to_voltage(skip_target, self.step_vol, quick_delay=skip_delay_num)
                    self.start_vol = skip_target
            else:
                msg.showwarning('warning',
                                 'Plot Range is enabled but Vmin/Vmax are invalid or equal.\n'
                                 'Please enter two different numbers, or uncheck Plot Range.')

        self.auto_run_process()

    def auto_run_process(self):
        # Safety guard: a previously-queued Tk `after` callback for this
        # function can still be sitting in the event queue and get
        # dispatched by window.update() while a ramp-down (or batch
        # teardown) is already in progress elsewhere. If that happens,
        # bail out immediately instead of issuing another up-step
        # setVoltage call on top of the active ramp-down.
        if getattr(self, 'rmp_dwn_flag', 0) == 1 or self.run_flag == 0:
            return

        temp, humid = 0.0, 0.0
        temp1, humid1 = '', ''
        if self.var.get() == 1:
            temp1, humid1 = self.run_arduino()
            if temp1 == '-999': self.stop_run()
            elif temp1 == '-998':
                retries=0 
                while temp1 == '-998' and retries < 20:
                    temp1, humid1 = self.run_arduino()
                    time.sleep(0.3)
                    #retries += 1
            temp = float(temp1)
            humid = float(humid1)

        try:
            if abs(self.end_vol - self.start_vol) >= 1e-3:       
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

                ################################################# Warning msg  block #######################################
                
                volt_tol = max(3 * self.step_vol, 0.05)           # V
                curr_tol_nA = max(self.curr_th * 1e9 * 0.02, 1.0)  # nA
                if (diff_V >= volt_tol or diff_I <= curr_tol_nA) and self.warn_flag == 0:
                    warning_message = (
                        'WARNING: Current limit reached\n'
                        f'Last measured current: {current_tmp_store:.1f} nA\n'
                        f'Voltage diff: {diff_V:.3f} V   Current diff: {diff_I:.3f} nA'
                    )
                    print(warning_message)
                    self.warn_flag = 1
                    self.run_flag = 0  # clear so batch loop exits on current limit too

                    self.xp.append(voltage_tmp)
                    self.xp_ap.append(self.start_vol)
                    self.yp.append(current_tmp_store)
                    self.curr_std_arr.append(0.0)  # keep arrays same length for error bars
                    self.temp_arr.append(float(temp))
                    self.humid_arr.append(float(humid))
                    self.time_arr.append(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

                    if self.scale_var.get() == 'log':
                        y_plot = [abs(v) for v in self.yp]
                    else:
                        y_plot = self.yp
                    self.plot1.set_data(self.xp, y_plot)
                    self.plot2.set_data(self.xp_ap, y_plot)
                    self.plot3.set_label(warning_message)
                    self.plot3.set_data([voltage_tmp], [y_plot[-1]])
                    h1, l1 = self.ax.get_legend_handles_labels()
                    h2, l2 = self.ax2.get_legend_handles_labels()
                    if h1 or h2:
                        self.ax.legend(h1 + h2, l1 + l2, fontsize=6, framealpha=0.9)
                    self._apply_yscale_and_limits()
                    self.refresh_errorbars()
                    self.ax.relim()
                    self.ax.autoscale_view()
                    self.figure_canvas.draw()

                    if self.batch_mode_active == False:
                        self.ramp_down_complete = False
                        self.ramp_down_zero(self.down_step_vol, self.time_delay)
                        # NOTE: save_results() is no longer called here -- it
                        # now runs inside _wait_ramp_then_analysis, *after*
                        # the forward/reverse fit has been computed, so the
                        # saved Analysis PNG actually contains the fit
                        # (previously it saved before the fit existed).
                        self.window.after(200, self._wait_ramp_then_analysis)
                        self.window.after(0, self.show_yellow_light)
                    else:
                        # Ramp-down/save for THIS channel is already handled
                        # by run_batch_sequence()'s wait-loop as soon as it
                        # sees run_flag == 0 (same as every other end-of-
                        # sweep path in batch mode) -- don't duplicate it.
                        cont = msg.askquestion(
                            "Current Limit Reached",
                            warning_message + "\n\nDo you want to continue with the rest of the batch?"
                        ).lower()
                        if cont in ('no', 'n'):
                            self.batch_stop_flag = 1
                    return
  
                ##################################################################################################################
                    
                if abs(voltage_tmp)>0.0001 and abs(current_tmp)>0:
                 voll_avg,curr_avg,curr_std=self.measure_all()
                 self.xp.append(voll_avg) 
                 self.xp_ap.append(self.start_vol)
                 self.yp.append(curr_avg * 1e9) 
                 self.curr_std_arr.append(curr_std * 1e9)  
                 self.temp_arr.append(float(temp))
                 self.humid_arr.append(float(humid))
                 self.time_arr.append(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
                 if self.scale_var.get() == 'log':
                     y_plot = [abs(v) for v in self.yp]
                 else:
                     y_plot = self.yp
                 self.plot1.set_data(self.xp, y_plot)
                 self.plot2.set_data(self.xp_ap, y_plot)
                 self.plot5.set_data(self.xp, self.temp_arr)
                 self.plot6.set_data(self.xp, self.humid_arr)
                 self._apply_yscale_and_limits()
                 self.refresh_errorbars()
                 self.figure_canvas.draw()

                if self.check_output_state() == 1:
                    if self.pause_plot == 0 and self.stop_flag == 0:
                        self.window.after(100, self.auto_run_process)
                    elif self.stop_flag == 1:
                        self.ramp_down_complete = False
                        if self.batch_mode_active==False:
                            try:
                                self.ramp_down_zero(self.down_step_vol, self.time_delay)
                            except Exception as e:
                                print(f"[Stop] ramp_down_zero failed: {e}; forcing output OFF directly.")
                                try:
                                    self.instrument.write("OUTP OFF")
                                except Exception as e2:
                                    print(f"[Stop] forced OUTP OFF also failed: {e2}")
                                self.rmp_dwn_flag = 0
                                self.ramp_down_complete = True
                                
                            # NOTE: save_results() now runs inside
                            # _wait_ramp_then_analysis, after the fit is
                            # computed (see note above).
                            self.window.after(200, self._wait_ramp_then_analysis)
                            self.window.after(0, self.show_yellow_light)
                        return
                    else:
                        self.window.after(0, self.show_yellow_light)
                        return

                restrict_done = False
                if self._restrict_plot_active and self._restrict_plot_vmax is not None:
                    if self._restrict_plot_ascending:
                        restrict_done = voltage_tmp >= self._restrict_plot_vmax - 1e-3
                    else:
                        restrict_done = voltage_tmp <= self._restrict_plot_vmax + 1e-3

                if (abs(voltage_tmp - self.end_vol) < 1e-3 or abs(self.start_vol - self.end_vol) < 1e-3
                        or restrict_done):
                    self.run_flag = 0  
                    self.ramp_down_complete = False
                    if self.batch_mode_active==False:
                        self.ramp_down_zero(self.down_step_vol, self.time_delay)
                        
                        # NOTE: save_results() now runs inside
                        # _wait_ramp_then_analysis, after the fit is
                        # computed (see note above).
                        self.window.after(200, self._wait_ramp_then_analysis)
                        self.window.after(0, self.show_yellow_light)
                    return

                if True:
                    if self.start_vol <= self.end_vol:
                        self.start_vol = min(self.start_vol + self.step_vol, self.end_vol)
                    else:
                        self.start_vol = max(self.start_vol - self.step_vol, self.end_vol)

            else:
                self.run_flag = 0  
                self.ramp_down_complete = False
                if self.batch_mode_active==False:
                    self.ramp_down_zero(self.down_step_vol, self.time_delay)
                    
                    # NOTE: save_results() now runs inside
                    # _wait_ramp_then_analysis, after the fit is
                    # computed (see note above).
                    self.window.after(200, self._wait_ramp_then_analysis)
                    self.window.after(1000, self.show_yellow_light)
                return

        except Exception as e:
            if self._looks_like_instrument_disconnect(e):
                self.handle_disconnect_during_sweep(e)
            else:
                msg.showerror("Error", f"{e}")                

    def ramp_up(self, voltage, vol_step=.50, sec_t=0.01):
        # Starting a new ramp: bump the generation token so that any
        # earlier ramp_up_run() chain still pending via window.after()
        # (e.g. an initial ramp-to-start_vol that a Batch Stop interrupts)
        # recognizes itself as stale on its next tick and quietly exits
        # instead of continuing to drive the instrument in parallel with
        # this new ramp. See _ramp_gen comment in __init__ for the full
        # explanation of the bug this prevents.
        self._ramp_gen += 1
        my_gen = self._ramp_gen
        self.instrument.write("*WAI")
        voltage_r1 = self.measure_voltage()
        current_r1 = self.measure_current() * 1000000000.0
        if self.var.get() == 1:
            temp, humid = self.run_arduino()
            self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(self.measure_voltage()) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r1) + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
            self.labels1.config(text=self.p_reading.get())
        else:
            self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(self.measure_voltage()) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r1))
            self.labels1.config(text=self.p_reading.get())

        indx = 0
        polar1 = self.chk_polarity(voltage_r1, voltage)
        if voltage_r1 > voltage: vol_step = -1.0 * vol_step
        if abs(voltage) < 0.5 and abs(voltage_r1) < 1:
            self.setVoltage(voltage)
            if (self.rmp_dwn_flag == 1):
                self.instrument.write("OUTP OFF")
                self.rmp_dwn_flag = 0
                # BUG FIX: this branch (target ~0V and already <1V) is
                # the ramp's *completion* just as much as the three
                # completion points inside ramp_up_run -- but unlike
                # those, it never set ramp_down_complete. Since
                # ramp_down_zero() always targets 0V, abs(voltage)<0.5
                # is always true, so this branch fires every time a
                # ramp-down starts from below 1V (e.g. a warning/STOP
                # that fires early in a sweep). Without this flag,
                # _wait_ramp_then_analysis (single mode) and the batch
                # loop's ramp-down wait both sit polling until their
                # 3600s timeout instead of proceeding immediately.
                self.ramp_down_complete = True
        else:
            # NOTE: previously this rounded vol_step to 1 decimal place
            # (round(vol_step, 1)). For any step under 0.05 V (e.g. the
            # common 0.02 V "Ramp Up" step), that rounds to exactly 0.0,
            # so ramp_up_run's "voltage_r1 = voltage_r1 + vol_step" never
            # advances -- the instrument gets re-commanded to the same
            # frozen voltage on every 1s tick, forever, while the outer
            # sweep loop keeps advancing self.start_vol regardless. Pass
            # vol_step through unrounded so the ramp actually progresses.
            self.ramp_up_run(voltage_r1, voltage, vol_step, polar1, sec_t, my_gen)

    def ramp_down_zero(self, v_step=1.0, delay_t=0.01):
        self.instrument.write("*WAI")
        voltage_r1 = self.measure_voltage()
        curr_r1 = self.measure_current()
        print('Threshod crossed Curr:: ', curr_r1*1e9, ' VOlTAGE:: ', voltage_r1)
        self.rmp_dwn_flag = 1
        end_volt = 0
        diff = abs(end_volt - voltage_r1)
        if (diff <= v_step): 
            v_step = v_step / 2.0
            print("Your Vstep ",v_step)
        self.ramp_up(end_volt, v_step, delay_t)

    def pause_plots(self, event=None):
        if (self.stop_flag == 0):
            if self.pause_plot == 0:
                self.pause_plot = 1
                self.pause.config(text='RESUME')
            else:
                # If we're paused because the power supply disconnected
                # mid-sweep, RESUME needs to search for it again and
                # reapply compliance before continuing -- just calling
                # auto_run_process() would be a no-op here since run_flag
                # was cleared when the disconnect was detected.
                if self.awaiting_reconnect:
                    self.attempt_reconnect_and_resume()
                    return
                self.pause_plot = 0
                self.pause.config(text='PAUSE')
                if self.sim_flag == 1: self.simulation()
                elif self.run_flag == 1: self.auto_run_process()
        else:
            msg.showwarning("warning", 'Run is stopped.Can\'t resume. Please start again.')

    def stop_run(self, event=None):
        should_save = (self.run_flag == 1) or (self.sim_flag == 1) or self.awaiting_reconnect
        self.awaiting_reconnect = False
        self.disconnect_resume_v = None
        self.pause_plot = 1
        self.stop_flag = 1
        self.pause.config(text='PAUSE', bg='#E0E0E0')
        if self.run_flag == 1:
            try:
                flag4, step_volt = self.is_number(self.step_voltage.get())
                flag5, delay_t = self.is_number(self.delay_time.get())
                if self.instrument: self.instrument.write("*WAI")
                self.ramp_down_complete = False
                self.ramp_down_zero(step_volt, delay_t)
            except Exception as e: print(f"Ramp Down Error: {e}")
        self.run_flag = 0
        self.sim_flag = 0
        
        if should_save:
            # FIX: Only run the single-mode finishing sequence if NOT in a batch
            if not self.batch_mode_active:
                # NOTE: save_results() is no longer scheduled separately
                # here. It now runs inside _wait_ramp_then_analysis, once
                # the ramp-down has actually finished and the forward/
                # reverse fit has been computed and drawn -- previously
                # it fired 100ms *before* _wait_ramp_then_analysis even
                # started polling for ramp completion, so the saved
                # Analysis PNG never contained the fit for a single
                # (non-batch) forward or reverse I-V run.
                self.window.after(600, self._wait_ramp_then_analysis)            

    def ramp_up_run(self, voltage_r1, voltage, vol_step, polar1, sec_t, my_gen=None):
        # Stale-chain guard: if a newer ramp has started since this
        # particular window.after() callback was scheduled, my_gen no
        # longer matches self._ramp_gen. Bail out silently -- do NOT touch
        # the instrument, rmp_dwn_flag, or ramp_down_complete, since those
        # belong to whichever ramp is current now. (my_gen defaults to
        # None for any external/legacy caller that doesn't pass it, in
        # which case the guard is skipped -- but every internal call site
        # in this file now supplies it.)
        if my_gen is not None and my_gen != self._ramp_gen:
            return
        if abs(voltage_r1 - voltage) > 1e-2:
            voltage_r1 = voltage_r1 + vol_step
            polar2 = self.chk_polarity(voltage_r1, voltage)
            if polar1 != polar2:
                voltage_r2, current_r2 = self.setVoltage(voltage)
                current_r2 = current_r2 * 1000000000.0
                if self.var.get() == 1:
                    temp, humid = self.run_arduino()
                    self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2) + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                    self.labels1.config(text=self.p_reading.get())
                else:
                    self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2))
                    self.labels1.config(text=self.p_reading.get())
                if (self.rmp_dwn_flag == 1):
                    self.instrument.write("OUTP OFF")
                    self.rmp_dwn_flag = 0
                    self.ramp_down_complete = True   #  ADD THIS
                return
            voltage_r2, current_r2 = self.setVoltage(voltage_r1)
            current_r2 = current_r2 * 1000000000.0
            time.sleep(sec_t)
            if self.var.get() == 1:
                temp, humid = self.run_arduino()
                self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2) + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                self.labels1.config(text=self.p_reading.get())
            else:
                self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2))
                self.labels1.config(text=self.p_reading.get())
            self.instrument.write("*WAI")
            if (abs(voltage_r2 - voltage) <= 1e-2):
                voltage_r2, current_r2 = self.setVoltage(voltage)
                current_r2 = current_r2 * 1000000000.0
                if self.var.get() == 1:
                    temp, humid = self.run_arduino()
                    self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2) + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                    self.labels1.config(text=self.p_reading.get())
                else:
                    self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2))
                    self.labels1.config(text=self.p_reading.get())
                if (self.rmp_dwn_flag == 1):
                    self.instrument.write("OUTP OFF")
                    self.ramp_down_complete = True
                    self.rmp_dwn_flag = 0
                    return
                return
            if (self.rmp_dwn_flag == 1): print('Ramping Down\nCurr:: ', current_r2, ' VOlTAGE:: ', voltage_r2)
            self.window.after(1000, lambda: self.ramp_up_run(voltage_r1, voltage, vol_step, polar1, sec_t, my_gen)) # Tanay Dey
        else:
            voltage_r2, current_r2 = self.setVoltage(voltage)
            current_r2 = current_r2 * 1000000000.0
            if self.var.get() == 1:
                temp, humid = self.run_arduino()
                self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2) + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
                self.labels1.config(text=self.p_reading.get())
            else:
                self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage_r2) + '\n' + 'CURRENT::' + self._fmt_current_nA(current_r2))
                self.labels1.config(text=self.p_reading.get())
            if (self.rmp_dwn_flag == 1):
                self.instrument.write("OUTP OFF")
                self.rmp_dwn_flag = 0
                self.ramp_down_complete = True
            return

    def _quick_skip_to_voltage(self, target_voltage, step_vol, quick_delay=0.1):
        """Ramp from the instrument's current voltage to `target_voltage`
        in `step_vol`-sized increments, used by the 'Plot Range
        (Vmin\u2192Vmax)' feature to move through the region below Vmin
        before recording/plotting starts.

        Unlike the normal per-step sweep in auto_run_process(), this does
        NOT call measure_all(), does NOT append to xp/yp/curr_std_arr, and
        does NOT touch the plot -- it only commands the source voltage and
        waits `quick_delay` seconds between steps. `quick_delay` comes from
        the dedicated "Skip Delay" field beside Vmax in the GUI (independent
        of the main sweep's "Delay (Sec)"), letting the skip phase move
        faster (or slower) than the recorded sweep. Blocking (pumps the Tk
        event loop so the GUI stays responsive), and aborts early on
        STOP / Batch Stop.
        """
        if getattr(self, 'instrument', None) is None:
            return
        try:
            quick_delay = float(quick_delay)
        except (TypeError, ValueError):
            quick_delay = 0.1
        if quick_delay < 0:
            quick_delay = 0.1

        try:
            voltage_now = self.measure_voltage()
        except Exception as e:
            print(f"[QuickSkip] Could not read starting voltage: {e}")
            return

        step = abs(step_vol) if step_vol else 0.5
        if step <= 0:
            step = 0.5
        direction = 1.0 if target_voltage >= voltage_now else -1.0
        v = voltage_now

        while (direction > 0 and v < target_voltage - 1e-6) or \
              (direction < 0 and v > target_voltage + 1e-6):
            if self.stop_flag == 1 or getattr(self, 'batch_stop_flag', 0):
                print("[QuickSkip] Aborted by Stop/Batch Stop.")
                return
            v = v + direction * step
            if (direction > 0 and v > target_voltage) or \
               (direction < 0 and v < target_voltage):
                v = target_voltage
            try:
                self.setVoltage(v)
            except Exception as e:
                print(f"[QuickSkip] setVoltage failed: {e}")
                return
            time.sleep(quick_delay)
            try:
                self.window.update()
            except Exception:
                pass

        try:
            self.setVoltage(target_voltage)
        except Exception as e:
            print(f"[QuickSkip] Final setVoltage failed: {e}")

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

    def start_test_dispatch(self, event=None):
        """Entry point for the main-GUI START TEST button. Branches to a
        single-channel run or a full batch run depending on the Run Mode
        radio buttons, replacing the old separate "START BATCH RUN" button
        that used to live inside the Batch Mode config window."""
        if self.run_mode_var.get() == 'batch':
            if self.batch_mode_active:
                msg.showwarning("Batch Mode", "A batch run is already in progress.")
                return
            if not self._batch_config:
                msg.showwarning("Batch Mode",
                                 "No batch configuration found.\n"
                                 "Open BATCH MODE, set up your boards/SiPMs, "
                                 "and click 'Confirm Batch Setup' first.")
                return
            self.run_batch_sequence()
        else:
            self.start_process()

    def sensel_current(self,indx):
        return self.current_array_sim[indx]*1e-9

    def simulation_run(self, event=None):
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
        self.curr_std_arr = []
        self.run_index = 0
        self.simulation()

    def simulation(self):
        self.window.after(0, self.show_green_light)
        temp = '25.4'
        humid = '60.62'

        voltage = self.voltage_array_sim[self.run_index]
        self.xp.append(voltage)
        self.temp_arr.append(float(temp))
        self.humid_arr.append(float(humid))
        self.time_arr.append(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

        # --- Dummy Gaussian measurement noise (simulation mode only) ---
        # The simulated I-V curve (self.current_array_sim) is a fixed,
        # noise-free lookup table, so on its own it can never feed the
        # error-bar code on the live/breakdown/Rq plots (curr_std_arr would
        # stay empty). To exercise those error bars in Simulation mode, a
        # fake per-point std-dev is generated here -- 4% of the current's
        # own magnitude plus a small absolute floor (so near-zero points
        # still get a visible whisker) -- and the plotted current itself is
        # jittered by a Gaussian sample with that std, so simulated points
        # scatter the way real noisy readings would.
        cur_clean = self.sensel_current(self.run_index) * 1e9  # nA, noise-free
        sim_std_nA = max(abs(cur_clean) * 0.04, 0.05)
        cur = float(np.random.normal(loc=cur_clean, scale=sim_std_nA))
        self.curr_std_arr.append(sim_std_nA)

        self.yp.append(cur)
        self.plot1.set_data(self.xp, self.yp)
        self.plot5.set_data(self.xp, self.temp_arr)
        self.plot6.set_data(self.xp, self.humid_arr)
        self.refresh_errorbars()

        if self.var.get() == 1:
            self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage) + '\n' + 'CURRENT::' + self._fmt_current_nA(cur) + "\n" + 'Temp:: ' + temp + ' \u00B0C  Humid:: ' + humid + ' %')
            self.labels1.config(text=self.p_reading.get())
        else:
            self.p_reading.set('VOLTAGE:: ' + self._fmt_voltage(voltage) + '\n' + 'CURRENT::' + self._fmt_current_nA(cur))
            self.labels1.config(text=self.p_reading.get())

        self.ax.relim()
        self.ax.autoscale_view()
        self.figure_canvas.draw()
        time.sleep(self.time_delay)
        self.run_index = self.run_index + 1

        if self.run_index < len(self.voltage_array_sim) and self.warn_flag == 0:
            if self.pause_plot == 0 and self.stop_flag == 0: self.window.after(100, self.simulation)
            else: self.window.after(0, self.show_yellow_light)
        else:
            xl=[]
            yl=[]
            xl.append(cur)
            yl.append(voltage)
            warning_message = 'WARNING: The current limit has been reached \n The voltage ramp-up is stopped \n Value of last measured current is '+str(round(cur,1)) +' nA'
            self.sim_flag = 0
            self.warn_flag = 1
            '''self.plot3.set_label(warning_message)
            self.plot3.set_data(xl, yl)

            self.window.after(0, self.show_red_light)
            if self.calc_vbd_var.get(): self.window.after(100, self.run_breakdown_analysis)
            self.save_results()'''
            self.plot3.set_label(warning_message)
            self.plot3.set_data(xl, yl)

            self.window.after(0, self.show_red_light)
            
            # --- FIX: Mode-Aware Analysis Trigger for Simulation ---
            # Run the fit synchronously (not via window.after) so it has
            # actually completed and been drawn onto self.fig_analysis
            # *before* save_results() saves it below. Previously the fit
            # was scheduled 100ms in the future while save_results() ran
            # immediately, so the saved Analysis PNG for a single
            # (non-batch) forward or reverse simulated run never
            # contained the fit.
            mode = self.analysis_mode_var.get()
            if mode == "reverse" and self.calc_vbd_var.get():
                try:
                    self.run_breakdown_analysis()
                except Exception as e:
                    print(f"[Analysis] Breakdown analysis error (ignored): {e}")
            elif mode == "forward" and self.show_rq_var.get():
                try:
                    self.run_quench_resistance_analysis()
                except Exception as e:
                    print(f"[Analysis] Rq analysis error (ignored): {e}")
            # -------------------------------------------------------
            
            self.save_results()            

    # ── Batch-level Pause/Resume/Stop (main GUI) ─────────────────────────
    # These are independent of the per-channel pause_plot/stop_flag, which
    # get reset to 0 at the start of every SiPM in run_batch_sequence (see
    # the per-channel reset block). batch_pause_flag / batch_stop_flag
    # persist across that reset, so they can actually control the whole
    # batch rather than just whichever channel happens to be running.
    def batch_pause_resume(self, event=None):
        if not self.batch_mode_active:
            msg.showinfo("Batch", "No batch is currently running.")
            return
        if self.batch_stop_flag:
            msg.showwarning("Batch", "Batch is stopping. Can't resume.")
            return
        # If the batch is paused because the power supply disconnected
        # mid-sweep, BATCH RESUME needs to search for it again and
        # reapply compliance before continuing, same as the single-mode
        # RESUME button -- toggling batch_pause_flag alone would leave
        # run_flag at 0 and the channel would just sit there.
        if self.awaiting_reconnect:
            self.attempt_reconnect_and_resume()
            return

        if not self.batch_pause_flag:
            self.batch_pause_flag = 1
            # Also pause whichever channel sweep is live right now, so the
            # pause is felt immediately rather than only at the next
            # inter-channel/inter-board checkpoint.
            self.pause_plot = 1
            if hasattr(self, 'batch_pause_btn'):
                self.batch_pause_btn.config(text='BATCH RESUME', bg='#F1C40F')
            if hasattr(self, '_batch_status_var') and self._batch_status_var:
                self._batch_status_var.set("Batch PAUSED — will resume after current step")
        else:
            self.batch_pause_flag = 0
            self.pause_plot = 0
            if hasattr(self, 'batch_pause_btn'):
                self.batch_pause_btn.config(text='BATCH PAUSE', bg='#E0E0E0')
            # Resume whatever live sweep was paused, same as the per-run
            # Resume button does.
            if self.sim_flag == 1:
                self.simulation()
            elif self.run_flag == 1:
                self.auto_run_process()

    def batch_stop(self, event=None):
        if not self.batch_mode_active:
            msg.showinfo("Batch", "No batch is currently running.")
            return
        confirm = msg.askyesno("Stop Batch",
                                "Stop the entire batch run?\n"
                                "The current SiPM will ramp down safely "
                                "before the batch exits.")
        if not confirm:
            return
        self.batch_stop_flag = 1
        # Release any pause so the loop can reach a checkpoint and unwind.
        self.batch_pause_flag = 0
        self.pause_plot = 0
        self.stop_flag = 1
        self.awaiting_reconnect = False
        self.disconnect_resume_v = None
        if hasattr(self, 'batch_pause_btn'):
            self.batch_pause_btn.config(text='BATCH PAUSE', bg='#E0E0E0')
        if hasattr(self, '_batch_status_var') and self._batch_status_var:
            self._batch_status_var.set("Batch stopping… ramping down current channel")

    def _batch_checkpoint(self):
        """Call inside every wait loop in run_batch_sequence. Blocks (while
        still pumping the Tk event loop) if batch_pause_flag is set, and
        raises BatchStopRequested if batch_stop_flag is set so the caller
        can unwind via the existing try/finally cleanup."""
        if self.batch_stop_flag:
            raise BatchStopRequested()
        while self.batch_pause_flag and not self.batch_stop_flag:
            self.window.update()
            time.sleep(0.1)
        if self.batch_stop_flag:
            raise BatchStopRequested()

    def _looks_like_instrument_disconnect(self, error):
        """Heuristic used to tell 'the power supply was physically
        unplugged / powered off mid-sweep' apart from an ordinary bug
        (a plotting error, a bad value, etc.) so only genuine
        communication failures trigger the reconnect-and-resume prompt.

        pyvisa raises VisaIOError for I/O-level failures -- timeout,
        resource no longer available, device not responding -- which is
        exactly what happens when a USB/GPIB SMU disconnects. Some
        backends surface a vanished USB-serial device node as a plain
        OSError instead, so that's treated the same way."""
        return isinstance(error, (visa.VisaIOError, OSError))

    def handle_disconnect_during_sweep(self, error):
        """Called from auto_run_process's except block when the
        power supply appears to have disconnected mid-I-V (single-mode
        or batch-mode -- both run through auto_run_process). Freezes the
        run and turns the PAUSE button into a RESUME button so the user
        can physically reconnect the power supply in their own time and
        then press RESUME to pick the sweep back up from exactly the
        voltage it was interrupted at, instead of restarting it.

        Previously this immediately popped an "search again now?" dialog.
        That fires before the user has necessarily plugged the supply
        back in, so answering "No" (because it isn't reconnected yet)
        permanently stopped the run, and there was no way to resume once
        it *was* reconnected -- the RESUME button just silently did
        nothing because run_flag had already been cleared. Now it just
        freezes and waits; all the actual reconnect work happens in
        attempt_reconnect_and_resume(), triggered by the RESUME /
        BATCH RESUME button.

        This relies on self.start_vol/self.end_vol/self.step_vol/
        self.down_step_vol/self.time_delay/self.curr_th being left
        untouched here; none of them are reset, so simply calling
        auto_run_process() again after reconnecting picks the sweep
        back up mid-stream. In batch mode, run_batch_sequence's
        per-channel wait loop checks self.awaiting_reconnect and holds
        in place (rather than tearing the channel down for ramp-down)
        until this resolves."""
        print(f"[Disconnect] Instrument communication error during I-V: {error}")

        # Freeze immediately so no further auto_run_process/ramp calls
        # try to talk to the now-unreachable instrument.
        self.run_flag = 0
        self.pause_plot = 1
        self.awaiting_reconnect = True
        self.disconnect_resume_v = self.start_vol
        self.pause.config(text='RESUME', bg='#F1C40F')
        if hasattr(self, 'batch_pause_btn') and self.batch_mode_active:
            self.batch_pause_btn.config(text='BATCH RESUME', bg='#F1C40F')
            if hasattr(self, '_batch_status_var') and self._batch_status_var:
                self._batch_status_var.set(
                    f"Batch PAUSED — power supply disconnected at "
                    f"{self.disconnect_resume_v:.3f} V")
        self.search_flag = 0
        self.instrument = None
        self.window.after(0, self.show_red_light)

        msg.showwarning(
            "Power Supply Disconnected",
            "Communication with the power supply was lost:\n"
            f"{error}\n\n"
            f"The sweep is paused at {self.disconnect_resume_v:.3f} V -- "
            "it has NOT been stopped.\n\n"
            "Reconnect the power supply, then press RESUME "
            f"{'(or BATCH RESUME) ' if self.batch_mode_active else ''}"
            "to search for it again and continue the sweep from where it "
            "left off. Press STOP instead if you want to end the run and "
            "keep the data collected so far."
        )

    def attempt_reconnect_and_resume(self):
        """Triggered by the RESUME (single-mode) or BATCH RESUME
        (batch-mode) button while self.awaiting_reconnect is set. Tries
        to (re)find the power supply and, if found, reapplies the
        current-limit compliance (a freshly opened VISA session doesn't
        retain the old session's setting) and resumes the sweep from
        self.disconnect_resume_v -- the voltage handle_disconnect_
        during_sweep froze at -- by simply calling auto_run_process()
        again, since start_vol/end_vol/step_vol were never reset.

        Returns True on a successful reconnect+resume, False otherwise
        (leaving the button on RESUME so the user can retry once the
        cable is actually back in)."""
        resume_v = self.disconnect_resume_v if self.disconnect_resume_v is not None else self.start_vol

        # Reuses the normal auto-detect connect flow.
        self.search()
        if self.search_flag != 1 or self.instrument is None:
            msg.showerror("Reconnect Failed",
                           "Could not find the power supply. Check the "
                           "cable/connection, then press RESUME "
                           f"{'/ BATCH RESUME ' if self.batch_mode_active else ''}"
                           "again -- the sweep is still paused, not lost.")
            self.window.after(0, self.show_red_light)
            return False

        try:
            self.set_current_threshold(self.curr_th)
        except Exception as e:
            print(f"[Disconnect] Could not reapply current limit after reconnect: {e}")

        msg.showinfo("Reconnected",
                      f"Power supply reconnected. Resuming I-V sweep from "
                      f"{resume_v:.3f} V.")

        self.awaiting_reconnect = False
        self.disconnect_resume_v = None
        self.run_flag = 1
        self.pause_plot = 0
        self.stop_flag = 0
        self.pause.config(text='PAUSE', bg='#E0E0E0')
        if hasattr(self, 'batch_pause_btn') and self.batch_mode_active:
            self.batch_pause_btn.config(text='BATCH PAUSE', bg='#E0E0E0')
            if hasattr(self, '_batch_status_var') and self._batch_status_var:
                self._batch_status_var.set(
                    f"Batch RESUMED — continuing from {resume_v:.3f} V")
        self.window.after(0, self.show_green_light)
        self.auto_run_process()
        return True                

    def HVTEST(self):
        is_forward = (self.analysis_mode_var.get() == "forward")
        if is_forward:
            user_response = msg.askquestion(
                "Positive IV TEST (Forward)",
                "Positive IV test selected, Forward mode is active \n"
                f"End voltage will be set to +{self.FORWARD_MODE_MAX_V:.0f} V.\n Do you want to continue?"
            ).lower()
        else:
            user_response = msg.askquestion("Positive IV TEST", "Positive IV test selected \n Do you want to continue?").lower()
        if user_response in ('no', 'n'):
            self.user_answer.set('')
            self._radio_last_value[str(self.user_answer)] = ''
            return
        else:
            self.start_voltage.set('0')
            if is_forward:
                self.end_voltage.set(str(self.FORWARD_MODE_MAX_V))
                self.step_voltage.set('0.05')
                self.down_step_voltage.set('0.5')
            else:
                self.end_voltage.set('30')
                self.step_voltage.set('0.5')
                self.down_step_voltage.set('5')
            self.current_th.set('10')
            self.delay_time.set('0.5')

    def IVTEST(self):
        is_forward = (self.analysis_mode_var.get() == "forward")
        if is_forward:
            user_response = msg.askquestion(
                "Negative IV TEST (Forward)",
                "Negative IV test selected, Forward mode is active \n"
                f"End voltage will be set to -{self.FORWARD_MODE_MAX_V:.0f} V.\n Do you want to continue?"
            ).lower()
        else:
            user_response = msg.askquestion("Negative IV TEST", "Negative IV test selected \n Do you want to continue?").lower()
        if user_response in ('no', 'n'):
            self.user_answer.set('')
            self._radio_last_value[str(self.user_answer)] = ''
            return
        else:
            self.start_voltage.set('0')
            if is_forward:
                self.end_voltage.set(str(-self.FORWARD_MODE_MAX_V))
                self.step_voltage.set('0.05')
                self.down_step_voltage.set('0.5')
            else:
                self.end_voltage.set('-30')
                self.step_voltage.set('0.5')
                self.down_step_voltage.set('5')
            self.current_th.set('10')
            self.delay_time.set('0.5')
            return

    def RUN_IV_HV(self):
        self.run_time_flag = 0
        if (self.is_blank_string(self.p_address.get()) == True or self.search_flag == 0):
            msg.showwarning('warning', 'Power supply is not detected \n SEARCH OR SET SOURCE ADDRESS')
            return 0
        if (self.user_answer.get() == ''):
            msg.showwarning('warning', 'Please choose any option from \nTEST TYPE HV/IV')
            return 0
        if (self.analysis_mode_var.get() == ''):
            msg.showwarning('warning', 'Please select a Mode: Forward or Reverse before starting the test.')
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

            if not self.check_forward_voltage_limit(start_voltage_num, end_voltage_num):
                return 0

            self.run_time_flag = 1
            return 1
        except Exception as e:
            msg.showerror("Error", f"{e}")
            return 0

    def _wait_ramp_then_analysis(self, timeout=3600):
        if not hasattr(self, '_ramp_wait_start'):
            self._ramp_wait_start = time.time()
            
        if self.ramp_down_complete:
            del self._ramp_wait_start
            
            # --- Automatic Analysis Trigger based on Mode ---
            if not self.batch_mode_active:
                mode = self.analysis_mode_var.get()
                if mode == "reverse" and self.calc_vbd_var.get():
                    try:
                        self.run_breakdown_analysis()
                    except Exception as e:
                        print(f"[Analysis] Breakdown analysis error (ignored): {e}")
                elif mode == "forward" and self.show_rq_var.get():
                    try:
                        self.run_quench_resistance_analysis()
                    except Exception as e:
                        print(f"[Analysis] Rq analysis error (ignored): {e}")

            # --- FIX: save AFTER the fit above has run, not before ---
            # save_results() writes self.fig_analysis out to PNG. Every
            # caller that schedules/polls this function (auto_run_process,
            # stop_run) used to call save_results() itself *before*
            # scheduling this function, so the forward (Rq) or reverse
            # (Vbd) fit above hadn't been computed yet -- self.fig_analysis
            # was still empty or held a stale fit from a previous run, and
            # the saved "*_Analysis_Graph.png" for a single (non-batch)
            # run never showed the current fit. Doing the save here, once
            # the fit has actually been drawn onto self.fig_analysis, is
            # what batch mode already does correctly via _batch_autosave
            # (called after its own analysis step).
            if not self.batch_mode_active:
                self.save_results()

            # ── Single (non-batch) I-V run finished: tell the user, then
            #    drop the status from "Complete" to "Paused / Hold" once
            #    they press OK. (Batch mode shows its own pop-up at the
            #    end of run_batch_sequence, so it's excluded here.)
            if not self.batch_mode_active:
                self.show_complete_light()
                msg.showinfo("I-V Run Complete", "Single I-V run complete!")
                self.show_yellow_light()
                
        elif time.time() - self._ramp_wait_start > timeout:
            del self._ramp_wait_start
            print("[Analysis] Ramp-down wait timed out; skipping analysis.")
        else:
            self.window.after(200, self._wait_ramp_then_analysis)

    def _fmt_voltage(self, v):
        """Auto-scale a voltage to V / mV / \u00b5V for the real-time display."""
        try:
            v = float(v)
            av = abs(v)
            if av < 1e-6:
                return "0.000 mV"
            if av < 1e-3:
                return f"{v * 1e6:.3f} \u00b5V"
            if av < 1:
                return f"{v * 1e3:.3f} mV"
            return f"{v:.3f} V"
        except (ValueError, TypeError):
            return "0.000 mV"

    def _fmt_current_nA(self, i_nA):
        """Auto-scale a current (given in nA, matching this app's existing
        convention) to pA / nA / \u00b5A / mA for the real-time reading
        display, so small currents don't round to '0.000 nA' and vanish."""
        ai = abs(i_nA)
        if ai == 0:
            return "0.000 nA"
        if ai < 1:
            return f"{i_nA * 1000.0:.3f} pA"
        if ai >= 1e6:
            return f"{i_nA / 1e6:.3f} mA"
        if ai >= 1e3:
            return f"{i_nA / 1e3:.3f} \u00b5A"
        return f"{i_nA:.3f} nA"

    def save_results(self):
        try:
            std_tmp = list(self.curr_std_arr) if hasattr(self, 'curr_std_arr') and len(self.curr_std_arr) == len(self.xp) else [0.0] * len(self.xp)
            alldata_tmp = pd.DataFrame({
                "VOLTS":           self.xp,
                "CURRNT_NAMP":     self.yp,
                "CURR_STD_NAMP":   std_tmp,
                "TEMP_DEGC":       self.temp_arr,
                "RH_PRCNT":        self.humid_arr,
                "TIME":            self.time_arr,
            })
            alldata_tmp.to_csv('temp.csv', index=False)
        except Exception:
            pass

        if self.batch_mode_active:
            return

        user_response = msg.askquestion("Save results",
                                        "Do You Want to save results?").lower()
        if user_response in ('no', 'n'):
            return

        outfile = self.module_name.get()
        outfile = re.sub(r'\s+', '', outfile)
        outfile = "".join(outfile.split())
        outfile = outfile.replace(":", "")
        self.current_datetimes.set(datetime.now().strftime("%d-%m-%Y-%H-%M"))
        directory = './Results/' + str(self.current_datetimes.get()) + '_' + outfile

        if os.path.exists(directory):
            user_response = msg.askquestion("Path Clashes",
                                            "Same Module \nDo You Want to Continue?").lower()
            if user_response in ('yes', 'y'):
                while os.path.exists(directory):
                    directory = directory + '_clone'
                os.makedirs(directory)
            else:
                return
        else:
            os.makedirs(directory)

        outfile = directory + '/' + str(self.current_datetimes.get()) + '_' + outfile + '_Result'
        log_file = outfile + '_Log.csv'

        std_log = list(self.curr_std_arr) if hasattr(self, 'curr_std_arr') and len(self.curr_std_arr) == len(self.xp) else [0.0] * len(self.xp)
        alldata = pd.DataFrame({
            "VOLTS":           self.xp,
            "CURRNT_NAMP":     self.yp,
            "CURR_STD_NAMP":   std_log,
            "TEMP_DEGC":       self.temp_arr,
            "RH_PRCNT":        self.humid_arr,
            "TIME":            self.time_arr,
        })
        alldata.to_csv(log_file, index=False)

        meas_plot = outfile + '_IV_Graph.png'
        self.figure.savefig(meas_plot)

        try:
            analysis_plot = outfile + '_Analysis_Graph.png'
            self.fig_analysis.savefig(analysis_plot)
        except Exception:
            pass

    # ------------------------------------------
    # 5.5 PLOTTING & ANALYSIS OVERLAYS
    # ------------------------------------------
    def plot_VI_graph(self, voltage_start, voltage_end):
        if hasattr(self, 'keithley_img_frame'): self.keithley_img_frame.pack_forget()
        if hasattr(self, 'analysis_artists'): self.analysis_artists = []

        self.figure.clf()
        self.errbar_container = None  # old artist was destroyed along with the cleared figure
        self.figure.subplots_adjust(left=0.12, right=0.88, top=0.90, bottom=0.25)

        self.ax = self.figure.add_subplot(111)
        self.ax2 = self.ax.twinx()

        self.plot1, = self.ax.plot([], [], 'o-', color='#3498DB', markersize=4, label="Measured I-V Data")
        self.plot2, = self.ax.plot([], [], 'x', color='#E74C3C', markersize=4, label="Validation of applied vs measured voltage")
        self.plot3, = self.ax.plot([], [], 'b', linestyle='None', label=None)
        self.plot4, = self.ax2.plot([], [], 'ro', linestyle='None', label=None)
        self.plot5, = self.ax2.plot([], [], 'bd', label="Temp")
        self.plot6, = self.ax2.plot([], [], 'ms', label="Humidity")

        if (self.sim_flag == 0):
            if self.var.get() == 1:
                self.ax2.set_visible(True)
                self.ax2.yaxis.set_label_position("right")
                self.ax2.yaxis.tick_right()
                #self.ax2.set_ylabel('Temp (\u00B0C) / Humidity (%)', color='m', fontsize=11, fontweight='bold', rotation=270, labelpad=20)
                self.ax2.text(1.08, 0.52, 'Temp (\u00B0C)', color='b', 
              fontsize=7, fontweight='bold', rotation=270, 
              transform=self.ax2.transAxes, ha='left', va='bottom')
                self.ax2.text(1.08, 0.48, ' / Humidity (%)', color='m', 
              fontsize=7, fontweight='bold', rotation=270, 
              transform=self.ax2.transAxes, ha='left', va='top')

            else:
                self.ax2.set_visible(False)
            self.ax.set_ylabel("Current\nin nA", fontsize=14, fontweight='bold', color='red', labelpad=6)
        else:
            if self.var.get() == 0:
                self.ax2.set_visible(False)
            else:
                self.ax2.set_visible(True)
                self.ax2.yaxis.set_label_position("right")
                self.ax2.yaxis.tick_right()
                #self.ax2.set_ylabel('Temp (\u00B0C) / Humidity (%)', color='m', fontsize=11, fontweight='bold', rotation=270, labelpad=20)
                self.ax2.text(1.08, 0.52, 'Temp (\u00B0C)', color='b', 
              fontsize=7, fontweight='bold', rotation=270, 
              transform=self.ax2.transAxes, ha='left', va='bottom')
                self.ax2.text(1.08, 0.48, ' / Humidity (%)', color='m', 
              fontsize=7, fontweight='bold', rotation=270, 
              transform=self.ax2.transAxes, ha='left', va='top')
            self.ax.set_ylabel("Current\nin nA", fontsize=14, fontweight='bold', color='#3498DB', labelpad=6)

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

        h1, l1 = self.ax.get_legend_handles_labels()
        h2, l2 = self.ax2.get_legend_handles_labels()
        if h1 or h2:
            self.ax.legend(h1+h2, l1+l2, bbox_to_anchor=(0.5, -0.15), loc='upper center', ncol=3, borderaxespad=0., fontsize=10, framealpha=0.9)

        self.figure_canvas = FigureCanvasTkAgg(self.figure, master=self.tab_measure)
        self.figure_canvas.get_tk_widget().pack(anchor="center", fill=Tk.BOTH, expand=True)
        self.figure_canvas.draw()

        self.plot_notebook.select(self.tab_measure)

    def run_breakdown_analysis(self, target_ax_iv=None, target_ax_prob=None, title=None):
        if not self.xp or len(self.xp) < 5: return
        volts = np.array(self.xp)
        if self.sim_flag == 1:
        	currents_nA = np.array(self.yp)*1000000000.
        else:
                currents_nA = np.array(self.yp)

        # Per-point current std-dev (nA), index-aligned with self.xp/self.yp,
        # used to weight the fit so noisier points count for less.
        current_std = None
        if hasattr(self, 'curr_std_arr') and self.curr_std_arr and len(self.curr_std_arr) == len(volts):
            current_std = np.array(self.curr_std_arr, dtype=float)

        v_bd_deriv = find_vbd_derivative(volts, currents_nA)
        popt, success, perr = optimize_fit(volts, currents_nA, v_bd_deriv, user_params=self.user_fit_params, current_std=current_std)
        self.plot_analysis_results(volts, currents_nA, v_bd_deriv, popt, success, target_ax_iv, target_ax_prob, title, perr=perr, current_std=current_std)        

    def run_quench_resistance_analysis(self, title=None):
        """
        Extract the SiPM quenching resistance Rq from a forward-bias I-V curve.
        Uses a ROOT-style Chi-square minimization fit (I vs V) via scipy.optimize.curve_fit.
        Rq is calculated as 1/slope, and its error is propagated as sigma_m / m^2.
        """
        if not self.xp or len(self.xp) < 6:
            msg.showwarning("Rq Analysis",
                            "Need at least 6 data points. Run a forward-bias I-V sweep first.")
            return

        volts = np.array(self.xp, dtype=float)
        currents_nA = np.array(self.yp, dtype=float)
        if self.sim_flag == 1:
            currents_nA = currents_nA * 1e9

        if hasattr(self, 'curr_std_arr') and self.curr_std_arr and len(self.curr_std_arr) == len(volts):
            curr_std_nA = np.array(self.curr_std_arr, dtype=float)
        else:
            curr_std_nA = np.zeros_like(volts)

        pos_mask = currents_nA > 0
        if np.sum(pos_mask) < 6:
            msg.showwarning("Rq Analysis",
                            "Fewer than 6 points with positive current found.\n"
                            "Ensure the data contains a forward-bias sweep.")
            return
        
        volts = volts[pos_mask]
        currents_nA = currents_nA[pos_mask]
        curr_std_nA = curr_std_nA[pos_mask]

        sort_idx = np.argsort(volts)
        volts = volts[sort_idx]
        currents_nA = currents_nA[sort_idx]
        curr_std_nA = curr_std_nA[sort_idx]

        # ── Region Detection: Auto or Manual ────────────────────────────────
        rq_mode = self.rq_mode_var.get() if hasattr(self, 'rq_mode_var') else 'auto'

        if rq_mode == 'manual':
            def _parse_v(var, fallback):
                try:
                    return float(var.get())
                except (ValueError, AttributeError):
                    return fallback

            r1_vmin = _parse_v(self.rq_r1_vmin_var, volts[0])
            r1_vmax = _parse_v(self.rq_r1_vmax_var, volts[len(volts) // 4])
            r2_vmin = _parse_v(self.rq_r2_vmin_var, volts[3 * len(volts) // 4])
            r2_vmax = _parse_v(self.rq_r2_vmax_var, volts[-1])

            mask1 = (volts >= r1_vmin) & (volts <= r1_vmax)
            mask2 = (volts >= r2_vmin) & (volts <= r2_vmax)

            if np.sum(mask1) < 2 or np.sum(mask2) < 2:
                msg.showwarning("Rq Manual Mode", "Ensure both regions contain at least 2 points.")
                return

            V1, I1, I1_std = volts[mask1], currents_nA[mask1], curr_std_nA[mask1]
            V2, I2, I2_std = volts[mask2], currents_nA[mask2], curr_std_nA[mask2]

            knee_v_mid = (r1_vmax + r2_vmin) / 2.0
            knee_idx = int((np.abs(volts - knee_v_mid)).argmin())
            V_knee, I_knee = volts[knee_idx], currents_nA[knee_idx]
            curve_start = int((np.abs(volts - r1_vmax)).argmin())
            curve_end   = int((np.abs(volts - r2_vmin)).argmin())

            r1_label = f"{r1_vmin:.2f} V → {r1_vmax:.2f} V"
            r2_label = f"{r2_vmin:.2f} V → {r2_vmax:.2f} V"
            #anno_header = r"$\bf{Manual\ Rq\ Analysis}$" + "\n"
            anno_header = ""#r"$\bf{Rq\ Analysis}$" + "\n"

        else:
            knee_idx, curve_start, curve_end = self._find_knee_region(volts, currents_nA)

            idx1_end = max(3, curve_start + 1)
            V1, I1, I1_std = volts[:idx1_end], currents_nA[:idx1_end], curr_std_nA[:idx1_end]

            idx2_start = min(len(volts) - 3, curve_end)
            V2, I2, I2_std = volts[idx2_start:], currents_nA[idx2_start:], curr_std_nA[idx2_start:]

            V_knee, I_knee = volts[knee_idx], currents_nA[knee_idx]
            r1_label = f"start → {volts[curve_start]:.2f} V"
            r2_label = f"{volts[curve_end]:.2f} V → end"
            anno_header = r"$\bf{Automatic\ Rq\ Analysis}$" + "\n"

        def _linear_fit(I_arr, V_arr, I_std_arr=None):
            """
            ROOT-style Fit: I = m * V + c 
            Rq = 1/m. Minimizes Chi-Square directly via curve_fit.
            """
            def fit_func(v, m, c):
                return m * v + c

            has_std = I_std_arr is not None and np.any(I_std_arr > 0) and len(I_arr) > 2

            if has_std:
                floor = 0.05 * np.median(np.abs(I_arr))
                if not np.isfinite(floor) or floor <= 0: floor = 1e-6
                sigma_I = np.where((I_std_arr <= 0) | ~np.isfinite(I_std_arr), floor, I_std_arr)

                try:
                    popt, pcov = curve_fit(fit_func, V_arr, I_arr, sigma=sigma_I, absolute_sigma=True)
                    m, c = popt
                    m_err = float(np.sqrt(pcov[0, 0]))
                except Exception:
                    popt, pcov = curve_fit(fit_func, V_arr, I_arr)
                    m, c = popt
                    m_err = float(np.sqrt(pcov[0, 0])) if np.isfinite(pcov[0,0]) else 0.0
            else:
                popt, pcov = curve_fit(fit_func, V_arr, I_arr)
                m, c = popt
                m_err = float(np.sqrt(pcov[0, 0])) if np.isfinite(pcov[0,0]) else 0.0
                sigma_I = np.ones_like(I_arr)

            try:
                R_bias = float(self.rq_rbias_var.get())
                N = float(self.rq_ncells_var.get())
            except:
                R_bias, N = 0.0, 1.0

            # R_total is 1/slope (Ohms). 
            # Formula: Rq = N * (R_total - R_bias)
            R_total_ohm = (1.0 / m) * 1e9
            Rq_ohm = N * (R_total_ohm - R_bias)
            Rq_kohm = Rq_ohm / 1e3
            
            # Error propagation: 
            # sigma_Rtotal = sigma_m / m^2
            # sigma_Rq = N * sigma_Rtotal
            sigma_Rtotal = (m_err / (m**2)) * 1e9
            Rq_ohm_err = N * sigma_Rtotal
            Rq_kohm_err = Rq_ohm_err / 1e3

            I_pred = fit_func(V_arr, m, c)
            
            # Calculate exact Chi2 / NDF
            if has_std:
                chi2 = np.sum(((I_arr - I_pred) / sigma_I) ** 2)
                ndf = max(1, len(I_arr) - 2)
                chi2_ndf = chi2 / ndf
            else:
                chi2_ndf = 0.0

            return m, c, chi2_ndf, Rq_ohm, Rq_kohm, I_pred, Rq_ohm_err, Rq_kohm_err, m_err, R_total_ohm, sigma_Rtotal

        # Execute Chi-Square fits
        (m1, c1, chi2_1, Rq1_ohm, Rq1_kohm, I_pred1, Rq1_ohm_err, Rq1_kohm_err,
         m1_err, R1_total_ohm, R1_total_ohm_err) = _linear_fit(I1, V1, I1_std)
        (m2, c2, chi2_2, Rq2_ohm, Rq2_kohm, I_pred2, Rq2_ohm_err, Rq2_kohm_err,
         m2_err, R2_total_ohm, R2_total_ohm_err) = _linear_fit(I2, V2, I2_std)

        # ── Resistance display formatting
        def _fmt_resistance(val_ohm, err_ohm):
            unit_pref = self.rq_unit_var.get() if hasattr(self, 'rq_unit_var') else 'auto'
            try:
                prec = int(self.rq_precision_var.get()) if hasattr(self, 'rq_precision_var') else 3
            except (ValueError, TypeError):
                prec = 3
            prec = max(0, min(9, prec))

            if unit_pref == 'ohm':
                scale, sym = 1.0, '\u03a9'
            elif unit_pref == 'kohm':
                scale, sym = 1e3, 'k\u03a9'
            elif unit_pref == 'mohm':
                scale, sym = 1e6, 'M\u03a9'
            else:
                aval = abs(val_ohm)
                if aval >= 1e6:
                    scale, sym = 1e6, 'M\u03a9'
                elif aval >= 1e3:
                    scale, sym = 1e3, 'k\u03a9'
                else:
                    scale, sym = 1.0, '\u03a9'

            val_disp = val_ohm / scale
            err_disp = err_ohm / scale
            return f"{val_disp:.{prec}f} \u00b1 {err_disp:.{prec}f} {sym}", val_disp, err_disp, sym, prec

        rq1_disp_str, _, _, _, _rq_prec = _fmt_resistance(Rq1_ohm, Rq1_ohm_err)
        rtotal1_disp_str, _, _, _, _ = _fmt_resistance(R1_total_ohm, R1_total_ohm_err)
        rq2_disp_str, _, _, _, _ = _fmt_resistance(Rq2_ohm, Rq2_ohm_err)
        rtotal2_disp_str, _, _, _, _ = _fmt_resistance(R2_total_ohm, R2_total_ohm_err)

        print(f"[Rq] Mode: {rq_mode.upper()}  |  Knee at V={V_knee:.3f} V, I={I_knee:.3f} nA  (index {knee_idx}/{len(volts)-1})")
        print(f"[Rq] Region 1 ({r1_label}): slope={m1:.4e} nA/V  Rq={Rq1_kohm:.3f} +/- {Rq1_kohm_err:.3f} kOhm  Chi2/ndf={chi2_1:.4f}")
        print(f"[Rq] Region 2 ({r2_label}): slope={m2:.4e} nA/V  Rq={Rq2_kohm:.3f} +/- {Rq2_kohm_err:.3f} kOhm  Chi2/ndf={chi2_2:.4f}")

        region_choice = self.rq_region_display_var.get() if hasattr(self, 'rq_region_display_var') else 'both'
        show_r1 = region_choice in ('both', 'region1')
        show_r2 = region_choice in ('both', 'region2')

        # ── Plot on Analysis Tab (ax_main) ──────────────────────────────────
        self.fig_analysis.clf()
        ax_main = self.fig_analysis.add_subplot(111)
        self.fig_analysis.subplots_adjust(left=0.12, right=0.95, top=0.80, bottom=0.12)

        show_errbars = (not hasattr(self, 'show_errorbars_var')) or self.show_errorbars_var.get()
        capsize_val = int(self.errorbar_capsize_var.get()) if hasattr(self, 'errorbar_capsize_var') and self.errorbar_capsize_var.get().isdigit() else 3

        if show_errbars and np.any(curr_std_nA > 0):
            err_scale = self.get_errorbar_scale() if hasattr(self, 'get_errorbar_scale') else 1.0
            ax_main.errorbar(volts, currents_nA, yerr=curr_std_nA * err_scale, fmt='o', color='indigo',
                              ecolor='indigo', elinewidth=1.0, capsize=capsize_val, markersize=6,
                              alpha=0.6, label="Measured Data")
        else:
            ax_main.plot(volts, currents_nA, 'o', color='indigo',
                         markersize=6, alpha=0.6, linestyle='None', label="Measured Data")

        if show_r1:
            ax_main.plot(V1, I1, 's', color='steelblue', markersize=7, alpha=0.9, linestyle='None', label=f"Selected Linear Region") #Region 1 (Linear)->Selected Linear Region
            V1_line = np.linspace(max(0.0, V1.min() * 0.9), V1.max() * 1.05, 200)
            I1_line_pred = m1 * V1_line + c1
            ax_main.plot(V1_line, I1_line_pred, color='steelblue', linewidth=2.0, linestyle='--',
                         label=f"Fit R1: Rq={rq1_disp_str}")

        if show_r2:
            ax_main.plot(V2, I2, 'D', color='darkorange', markersize=7, alpha=0.9, linestyle='None', label=f"Selected Linear Region") #Region 2 (Linear)->Selected Linear Region
            V2_line = np.linspace(max(0.0, V2.min() * 0.9), V2.max() * 1.05, 200)
            I2_line_pred = m2 * V2_line + c2
            ax_main.plot(V2_line, I2_line_pred, color='darkorange', linewidth=2.0, linestyle='--',
                         label=f"Fit R2: Rq={rq2_disp_str}")

        if show_r1 and show_r2:
            ax_main.axvline(V_knee, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
            ax_main.plot(V_knee, I_knee, 'r^', markersize=10, zorder=5, label=f"Knee @ {V_knee:.2f} V")

        slope1_txt = f"  slope (m) = ({m1:.4e} \u00b1 {m1_err:.4e}) nA/V"
        r1_txt = (f"{slope1_txt}\n"
                   f"  $R_{{total}}$ = {rtotal1_disp_str}")#   $\\chi^2/ndf$ = {chi2_1:.2f}") #tanay

        anno = anno_header
        if show_r1 and show_r2:
            anno += f"Knee point: {V_knee:.3f} V\n\n"
        
        if show_r1:
            anno += r"$\bf{Region}$" + f"  ({r1_label})\n" + r1_txt + "\n" #tanay changed Region \2-> Region
            if show_r2: anno += "\n"
        if show_r2:
            slope2_txt = f"  slope (m) = ({m2:.4e} \u00b1 {m2_err:.4e}) nA/V"
            anno += (r"$\bf{Region}$" + f"  ({r2_label})\n" #tanay changed Region \2-> Region
                     f"{slope2_txt}\n"
                     f"  $R_{{total}}$ = {rtotal2_disp_str}")#,   $\\chi^2/ndf$ = {chi2_2:.2f}") #tanay

        ax_main.text(0.02, 0.97, anno, transform=ax_main.transAxes,
                     verticalalignment='top', fontsize=11,
                     bbox=dict(boxstyle="round", fc="white", alpha=0.92, ec="darkorange"), 
                     color="#2C3E50")

        eq_latex = r"$I = m \cdot V + c \ \ \Rightarrow\ \ R_{total} = \dfrac{1}{m}$"
        ax_main.text(0.40, 0.97, eq_latex, transform=ax_main.transAxes,
                     verticalalignment='top', fontsize=11, color='darkred', fontstyle='italic',
                     bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.92, ec="red"))

        try:
            rq_R_bias = float(self.rq_rbias_var.get())
        except (ValueError, AttributeError):
            rq_R_bias = 0.0
        try:
            rq_N_cells = float(self.rq_ncells_var.get())
        except (ValueError, AttributeError):
            rq_N_cells = 1.0

        rq_box = r"$\bf{Quenching\ Resistance}$" + "\n"
        rq_box += r"$R_q = N_{microcell} \times (R_{total} - R_{bias})$" + "\n"
        rq_box += f"  $N_{{microcell}}$ = {rq_N_cells:g}\n"
        rq_box += f"  $R_{{bias}}$ = {rq_R_bias:g} \u03a9\n"
        if show_r1:
            rq_box += f"  $R_q$  = {rq1_disp_str}\n"
        if show_r2:
            rq_box += f"  $R_q$  = {rq2_disp_str}\n"
        rq_box = rq_box.rstrip("\n")

        ax_main.text(0.98, 0.03, rq_box, transform=ax_main.transAxes,
                     verticalalignment='bottom', horizontalalignment='right', fontsize=11,
                     bbox=dict(boxstyle="round", fc="white", alpha=0.92, ec="steelblue"),
                     color="#2C3E50")

        ax_main.set_ylabel("Current (nA)", fontweight='bold', fontsize=13)
        ax_main.set_yscale(self.scale_var.get() or 'linear')
        
        title_suffix = "" #if (show_r1 and show_r2) else ("  (Region 1 only)" if show_r1 else "  (Region 2 only)") #tanay changed
        label_name = title if title else self.module_name.get()
        ax_main.set_title(f"Forward-Bias Quenching Resistance  –  {label_name}{title_suffix}",
                          fontsize=13, fontweight='bold', pad=45)        
        
        legend_ncol = 3 if (show_r1 and show_r2) else 2
        ax_main.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                       ncol=legend_ncol, mode="expand", borderaxespad=0., frameon=False, fontsize=10)
        ax_main.grid(True, which='both', linestyle='--', alpha=0.5)
        ax_main.set_xlabel("Bias Voltage (V)", fontweight='bold', fontsize=12)

        self.canvas_analysis.draw()

        # --- NEW: Overlay Rq result on the Measurement Tab ---
        try:
            # Check if we should overlay the text (always True for single, follows user setting for batch)
            should_overlay = True
            if getattr(self, 'batch_mode_active', False) and hasattr(self, '_batch_config'):
                should_overlay = self._batch_config.get("overlay_fit_on_iv", True)

            if should_overlay and hasattr(self, 'ax') and self.ax is not None:
                # Compile a compact summary for the Measurement Tab
                meas_summary = ""
                if show_r1:
                    meas_summary += f"Rq: {rq1_disp_str}\n"
                    #meas_summary += f"Optocoupler Resistance: {rq1_disp_str}\n"
                if show_r2:
                    meas_summary += f"Rq: {rq2_disp_str}"
                    #meas_summary += f"Optocoupler Resistance: {rq2_disp_str}"
                meas_summary = meas_summary.strip()
                
                if meas_summary:
                    # Dynamically shrink font size for small batch tiles
                    f_size = 7 if getattr(self, 'batch_mode_active', False) else 11
                    
                    # Annotate self.ax (Measurement Tab active subplot)
                    self.ax.text(0.02, 0.98, meas_summary, transform=self.ax.transAxes,
                                 verticalalignment='top', fontsize=f_size,
                                 bbox=dict(boxstyle="round", fc="white", alpha=0.85, ec="steelblue"),
                                 color="black", zorder=10)
                    
                    # Overlay the actual fit lines on the Measurement Tab so they match the text
                    if show_r1:
                        self.ax.plot(V1_line, I1_line_pred, color='steelblue', linewidth=1.5, linestyle='--')
                    if show_r2:
                        self.ax.plot(V2_line, I2_line_pred, color='darkorange', linewidth=1.5, linestyle='--')
                        
                    # Redraw the Measurement canvas so the user sees it immediately
                    if hasattr(self, 'figure_canvas') and self.figure_canvas:
                        self.figure_canvas.draw()
        except Exception as e:
            print(f"[Rq] Failed to overlay fit data on Measurement tab: {e}")
        # -----------------------------------------------------

        # --- Stash the fit results so callers (e.g. batch autosave) can
        # redraw this same linear Rq fit elsewhere without recomputing it,
        # and without confusing it with the reverse-bias breakdown fit. ---
        self.last_rq_result = {
            'show_r1': show_r1, 'show_r2': show_r2,
            'V1': V1, 'I1': I1, 'm1': m1, 'c1': c1,
            'rq1_disp_str': rq1_disp_str, 'r1_label': r1_label,
            'V2': V2, 'I2': I2, 'm2': m2, 'c2': c2,
            'rq2_disp_str': rq2_disp_str, 'r2_label': r2_label,
            'V_knee': V_knee, 'I_knee': I_knee,
        }

        self.plot_notebook.select(self.tab_analysis)

    def plot_analysis_results(self, volts, currents_nA, v_bd_deriv, popt, success, target_ax_iv=None, target_ax_prob=None, title=None, perr=None, current_std=None):
        is_batch = target_ax_iv is not None
        if perr is None:
            perr = np.zeros(6)

        if hasattr(self, 'show_geiger_var'): show_geiger = self.show_geiger_var.get()
        else: show_geiger = True

        if not is_batch:
            self.fig_analysis.clf()
            self.fig_analysis.subplots_adjust(left=0.10, right=0.95, top=0.88, bottom=0.10, hspace=0.1)

            if show_geiger:
                gs = GridSpec(2, 1, height_ratios=[3, 1])
                ax_iv = self.fig_analysis.add_subplot(gs[0])
                ax_prob = self.fig_analysis.add_subplot(gs[1], sharex=ax_iv)
            else:
                ax_iv = self.fig_analysis.add_subplot(111)
                ax_prob = None
        else:
            ax_iv = target_ax_iv
            ax_prob = target_ax_prob
            ax_iv.clear()
            if ax_prob: ax_prob.clear()

        star = mpath.Path.unit_regular_star(6)
        circle = mpath.Path.unit_circle()
        cut_star = mpath.Path(vertices=np.concatenate([circle.vertices, star.vertices[::-1, ...]]), codes=np.concatenate([circle.codes, star.codes]))

        markersize_val = 4 if is_batch else 10
        ax_iv.plot(volts, currents_nA, marker=cut_star, color='indigo', markersize=markersize_val, alpha=0.6, label="Measured Data", linestyle='None')

        # Current-uncertainty error bars on the measured points (same
        # self.curr_std_arr-derived data used for the live Measurement tab
        # and the Rq analysis plot). Skipped on the small batch-grid tiles,
        # where the whiskers would just be visual noise at that size.
        if not is_batch and current_std is not None and np.any(np.asarray(current_std) > 0):
            show_errbars = (not hasattr(self, 'show_errorbars_var')) or self.show_errorbars_var.get()
            if show_errbars:
                try:
                    flag_cap, capsize_val = self.is_number(self.errorbar_capsize_var.get()) if hasattr(self, 'errorbar_capsize_var') else (False, 4)
                    if not flag_cap or capsize_val < 0:
                        capsize_val = 4
                    err_scale = self.get_errorbar_scale() if hasattr(self, 'get_errorbar_scale') else 1.0
                    yerr_scaled = np.asarray(current_std) * err_scale
                    ax_iv.errorbar(volts, currents_nA, yerr=yerr_scaled, fmt='none',
                                   ecolor='indigo', elinewidth=1.0, capsize=capsize_val,
                                   capthick=1.0, alpha=0.5, zorder=1)
                except Exception as e:
                    print(f"[Breakdown] Could not draw error bars: {e}")

        if success:
            v_bd_fit = popt[0]
            y_val_nA = dinu_eq8_model(v_bd_fit, *popt) 
            idx = (np.abs(volts - v_bd_fit)).argmin()
            ax_iv.plot(v_bd_fit, y_val_nA, 'rx', markersize=markersize_val, markeredgewidth=2, label="Breakdown Point")

            if not is_batch:
                ax_iv.annotate(f"Breakdown Point: {v_bd_fit:.2f}V", xy=(v_bd_fit, y_val_nA), xytext=(v_bd_fit - (max(volts)*0.15), currents_nA[idx]-y_val_nA/2), color='red', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), bbox=dict(boxstyle="round", fc="white", alpha=0.7), fontsize=13)
            else:
                ax_iv.annotate(f"{v_bd_fit:.2f}V", xy=(v_bd_fit, y_val_nA), xytext=(v_bd_fit - (max(volts)*0.1), currents_nA[idx]), color='red', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), fontsize=7)

            #################################################################################
            overvol=float(self.set_ovv.get())
            y_val_nA_ov = dinu_eq8_model(v_bd_fit+overvol, *popt)

            ax_iv.plot(v_bd_fit+overvol, y_val_nA_ov, 'mP', markersize=markersize_val, markeredgewidth=2,label=f"Current at {overvol:0.2f} Overvoltage")
            
            if not is_batch:
                ax_iv.annotate(f"$V_{{bd}}+overvol$: {v_bd_fit+overvol:.2f} V\n I: {y_val_nA_ov:0.0f} nA", xy=(v_bd_fit+overvol, y_val_nA_ov), xytext=(v_bd_fit+overvol -7, y_val_nA_ov-0.7*y_val_nA_ov), color='m', fontweight='bold', arrowprops=dict(arrowstyle='->', color='red'), bbox=dict(boxstyle="round", fc="white", alpha=0.7), fontsize=13)
            #################################################################################

            v_smooth = np.linspace(min(volts), min(max(volts), popt[1]-0.1), 1000)
            i_fit_nA = dinu_eq8_model(v_smooth, *popt)
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
                
                if not is_batch:
                    ax_prob.set_ylabel("Geiger Prob.", fontweight='bold', color='blue', fontsize=13)
                    ax_prob.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)
                else:
                    ax_prob.set_ylabel("Prob.", fontweight='bold', color='blue', fontsize=7)
                    ax_prob.set_xlabel("V", fontweight='bold',fontsize=7)
                    
                ax_prob.set_ylim(-0.05, 1.1)
                ax_prob.grid(True, which='both', linestyle='--', alpha=0.5)
                
                if not is_batch:
                    formula_txt = r"$P_{Geiger} = 1 - e^{-p(V - V_{bd})}$"
                    ax_prob.text(0.02, 0.6, formula_txt, transform=ax_prob.transAxes,
                                 fontsize=11, color='darkblue', fontweight='bold',
                                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85, ec="blue"))
                    ax_prob.tick_params(axis='both', labelsize=10)
                else:
                    ax_prob.tick_params(axis='both', labelsize=6)

            vbd_err_txt = f" $\\pm$ {perr[0]:.2f}" if perr[0] > 0 else ""
            if self.show_dcr_var.get()==0:
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f}{vbd_err_txt} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}")
            else:
                if abs(self.C_ucell)>0:
                    DCR=popt[3]*1e-9/(self.C_ucell*1e3)
                else: DCR=0
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f}{vbd_err_txt} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}\n"+f"DCR : {DCR:0.3f} kHz")

            equation_latex = (r"$I_{tot} = I_{leak} + I_{aval}$" + "\n" + r"$I_{leak} = \exp(aV + b)$" + "\n" + r"$I_{aval} = A \cdot \Delta V \cdot (1 - e^{-p \Delta V}) \cdot \frac{V_{cr}-V_{bd}}{V_{cr}-V}$" + "\n")

            if not is_batch:
                ax_iv.text(0.33, 0.96, equation_para, transform=ax_iv.transAxes, verticalalignment='top', fontsize=15, bbox=dict(boxstyle="round", fc="white", alpha=0.90, ec="#27AE60"), color="#2C3E50")
                ax_iv.text(0.01, 0.92, equation_latex, transform=ax_iv.transAxes, verticalalignment='top', fontsize=15, bbox=dict(boxstyle="round", fc="white", alpha=0.90, ec="green"), color="black")
            else:
                summary = f"Vbd:{popt[0]:.1f}V DCR:{DCR:.1f}kHz" if self.show_dcr_var.get() and abs(self.C_ucell)>0 else f"Vbd:{popt[0]:.1f}V"
                ax_iv.text(0.02, 0.98, summary, transform=ax_iv.transAxes, verticalalignment='top', fontsize=6, bbox=dict(boxstyle="round", fc="white", alpha=0.8), color="black")

        ax_iv.set_yscale(self.scale_var.get() or 'linear')
        ax_iv.grid(True, which='both', linestyle='--', alpha=0.5)

        if not is_batch:
            ax_iv.set_ylabel("Current (nA)", fontweight='bold',fontsize=14)
            if not show_geiger: ax_iv.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)
            ax_iv.set_title(title if title else self.module_name.get(), pad=35)
            ax_iv.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=4, mode="expand", borderaxespad=0., frameon=False, fontsize=14)
            self.canvas_analysis.draw()
            if not self.batch_mode_active:
                self.plot_notebook.select(self.tab_analysis)
            print(f"Analysis Complete. Vbd: {v_bd_deriv:.2f}V")
        else:
            ax_iv.set_ylabel("I (nA)", fontsize=7)
            if not show_geiger: ax_iv.set_xlabel("V", fontsize=7)
            title_text = title if title else self.module_name.get()
            ax_iv.set_title(title_text, pad=3, fontsize=8, fontweight='bold')
            ax_iv.tick_params(axis='both', labelsize=6)
            ax_iv.legend(loc='lower right', fontsize=5)

    def _find_knee_region(self, volts, currents_nA):
        """
        Detects the non-linear curve of a diode by fitting asymptotic lines to the 
        Ohmic (tail) and Leakage (head) regions, finding their intersection, 
        and checking where the real data deviates from linearity.
        """
        n = len(volts)
        # Fallback if there are too few points to analyze safely
        if n < 10:
            mid = n // 2
            return mid, max(0, mid-1), min(n-1, mid+1)
            
        # 1. Fit the Ohmic region (assume the last 25% of data is fully linear)
        tail_count = max(4, int(0.25 * n))
        m_ohm, c_ohm = np.polyfit(volts[-tail_count:], currents_nA[-tail_count:], 1)
        
        # 2. Fit the Leakage region (assume the first 25% is fully flat/leakage)
        head_count = max(4, int(0.25 * n))
        m_leak, c_leak = np.polyfit(volts[:head_count], currents_nA[:head_count], 1)
        
        # 3. Find the theoretical Knee (Intersection of the two straight lines)
        if m_ohm == m_leak:
            v_knee = volts[n//2]
        else:
            v_knee = (c_leak - c_ohm) / (m_ohm - m_leak)
            
        # Find the actual data index closest to this theoretical intersection
        knee_idx = (np.abs(volts - v_knee)).argmin()
        # Force knee_idx to not be at an extreme edge
        knee_idx = max(3, min(knee_idx, n - 4))
        
        # 4. Define a tolerance (e.g., 5% of the total current range)
        tolerance = 0.05 * (np.max(currents_nA) - np.min(currents_nA))
        
        # 5. Find curve_end: Walk backwards from max voltage. 
        # When data deviates from the Ohmic trendline by more than the tolerance, the curve begins.
        ohmic_line = m_ohm * volts + c_ohm
        curve_end = knee_idx + 1
        for i in range(n - 1, knee_idx, -1):
            if np.abs(currents_nA[i] - ohmic_line[i]) > tolerance:
                curve_end = min(n - 3, i + 1)
                break
                
        # 6. Find curve_start: Walk forwards from min voltage.
        # When data deviates from the Leakage trendline, the knee is starting.
        leak_line = m_leak * volts + c_leak
        curve_start = knee_idx - 1
        for i in range(0, knee_idx):
            if np.abs(currents_nA[i] - leak_line[i]) > tolerance:
                curve_start = max(2, i - 1)
                break
                
        # 7. Failsafe logical bounds in case of highly unusual noise spikes
        if curve_start >= curve_end or curve_start >= knee_idx or curve_end <= knee_idx:
            # Fall back to a fixed 15% window around the intersection point
            window = max(2, int(0.15 * n))
            curve_start = max(2, knee_idx - window)
            curve_end = min(n - 3, knee_idx + window)
            
        return knee_idx, curve_start, curve_end        

    # ----------------------------------------------------------------------
    # SECTION: QUENCH-RESISTANCE (Rq) CONFIG POPUP
    # PURPOSE: Toplevel window for configuring how the forward-bias
    #          quench-resistance analysis picks its linear fit region
    #          (Auto vs Manual) and per-SiPM overrides.
    # ----------------------------------------------------------------------
    def _open_rq_config_popup(self):
        popup = Tk.Toplevel(self.window)
        popup.title("Rq Analysis Configuration")
        popup.configure(bg=self.colors['bg_sidebar'])
        popup.resizable(False, False)
        popup.grab_set()

        # ... (Keep existing Title and Mode selection code) ...
        # Center over main window
        #self.window.update_idletasks()
        #px = self.window.winfo_x() + self.window.winfo_width() // 2 - 200
        #py = self.window.winfo_y() + self.window.winfo_height() // 2 - 140
        #popup.geometry(f"400x300+{px}+{py}")

        ttk.Label(popup, text="Quench Resistance – Region Mode",
                  style='Header.TLabel', padding=(8, 6)).pack(fill=Tk.X)

        # ── Mode radio buttons ──────────────────────────────────────────────
        mode_lf = ttk.LabelFrame(popup, text="Detection Mode", style='Group.TLabelframe', padding=6)
        mode_lf.pack(fill=Tk.X, padx=10, pady=(6, 4))

        _mode_var = StringVar(value=self.rq_mode_var.get())

        Radiobutton(mode_lf, text='Auto  (algorithm finds linear regions automatically)',
                    variable=_mode_var, value='auto',
                    bg=self.colors['bg_sidebar'], fg='white',
                    selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar']).pack(anchor='w')
        Radiobutton(mode_lf, text='Manual  (specify voltage range for each region)',
                    variable=_mode_var, value='manual',
                    bg=self.colors['bg_sidebar'], fg='white',
                    selectcolor=self.colors['bg_sidebar'],
                    activebackground=self.colors['bg_sidebar']).pack(anchor='w')

        # ── Added: Bias Resistance & Microcell Count ────────────────────────
        corr_lf = ttk.LabelFrame(popup, text="Correction Parameters", style='Group.TLabelframe', padding=6)
        corr_lf.pack(fill=Tk.X, padx=10, pady=(6, 4))
        
        # Initialize variables if they don't exist
        if not hasattr(self, 'rq_rbias_var'): self.rq_rbias_var = StringVar(value="0")
        if not hasattr(self, 'rq_ncells_var'): self.rq_ncells_var = StringVar(value="1")

        ttk.Label(corr_lf, text="R_bias (\u03A9):", style='Sidebar.TLabel').grid(row=0, column=0, sticky='e')
        ttk.Entry(corr_lf, textvariable=self.rq_rbias_var, width=10).grid(row=0, column=1, padx=5)
        
        ttk.Label(corr_lf, text="N (Microcells):", style='Sidebar.TLabel').grid(row=0, column=2, sticky='e')
        ttk.Entry(corr_lf, textvariable=self.rq_ncells_var, width=10).grid(row=0, column=3, padx=5)

        # ── Added: Display Unit & Precision ─────────────────────────────────
        disp_lf = ttk.LabelFrame(popup, text="Display Options", style='Group.TLabelframe', padding=6)
        disp_lf.pack(fill=Tk.X, padx=10, pady=(6, 4))

        _unit_var = StringVar(value=self.rq_unit_var.get())
        ttk.Label(disp_lf, text="Rq / R_total Unit:", style='Sidebar.TLabel').grid(row=0, column=0, sticky='w')
        unit_row = ttk.Frame(disp_lf, style='Sidebar.TFrame')
        unit_row.grid(row=1, column=0, columnspan=4, sticky='w', pady=(2, 6))
        for _u_val, _u_txt in (('auto', 'Auto'), ('ohm', 'Ω'), ('kohm', 'kΩ'), ('mohm', 'MΩ')):
            Radiobutton(unit_row, text=_u_txt, variable=_unit_var, value=_u_val,
                        bg=self.colors['bg_sidebar'], fg='white',
                        selectcolor=self.colors['bg_sidebar'],
                        activebackground=self.colors['bg_sidebar']).pack(side=Tk.LEFT, padx=(0, 8))

        _prec_var = StringVar(value=self.rq_precision_var.get())
        ttk.Label(disp_lf, text="Decimal Precision:", style='Sidebar.TLabel').grid(row=2, column=0, sticky='w')
        prec_spin = ttk.Spinbox(disp_lf, from_=0, to=9, textvariable=_prec_var, width=5, wrap=False)
        prec_spin.grid(row=2, column=1, sticky='w', padx=(4, 0))

        # ... (Keep Manual Bounds and OK/Cancel buttons) ...
        manual_lf = ttk.LabelFrame(popup, text="Manual Region Bounds (V)", style='Group.TLabelframe', padding=6)
        manual_lf.pack(fill=Tk.X, padx=10, pady=(0, 6))

        # Temp vars so cancel doesn't alter stored values
        _r1min = StringVar(value=self.rq_r1_vmin_var.get())
        _r1max = StringVar(value=self.rq_r1_vmax_var.get())
        _r2min = StringVar(value=self.rq_r2_vmin_var.get())
        _r2max = StringVar(value=self.rq_r2_vmax_var.get())

        for c in range(4): manual_lf.columnconfigure(c, weight=1)

        ttk.Label(manual_lf, text="Region :", style='Sidebar.TLabel', #tanay change Region 1- Region
                  font=('Segoe UI', 9, 'bold')).grid(row=0, column=0, columnspan=4, sticky='w')
        ttk.Label(manual_lf, text="Vmin:", style='Sidebar.TLabel').grid(row=1, column=0, sticky='e', padx=(0,2))
        ttk.Entry(manual_lf, textvariable=_r1min, width=8).grid(row=1, column=1, sticky='ew', padx=2)
        ttk.Label(manual_lf, text="Vmax:", style='Sidebar.TLabel').grid(row=1, column=2, sticky='e', padx=(4,2))
        ttk.Entry(manual_lf, textvariable=_r1max, width=8).grid(row=1, column=3, sticky='ew', padx=2)

        ttk.Label(manual_lf, text="Region :", style='Sidebar.TLabel', #tanay change Region 2- Region
        
                  font=('Segoe UI', 9, 'bold')).grid(row=2, column=0, columnspan=4, sticky='w', pady=(6,0))
        ttk.Label(manual_lf, text="Vmin:", style='Sidebar.TLabel').grid(row=3, column=0, sticky='e', padx=(0,2))
        ttk.Entry(manual_lf, textvariable=_r2min, width=8).grid(row=3, column=1, sticky='ew', padx=2)
        ttk.Label(manual_lf, text="Vmax:", style='Sidebar.TLabel').grid(row=3, column=2, sticky='e', padx=(4,2))
        ttk.Entry(manual_lf, textvariable=_r2max, width=8).grid(row=3, column=3, sticky='ew', padx=2)

        # ── OK / Cancel ─────────────────────────────────────────────────────
        btn_row = ttk.Frame(popup, style='Sidebar.TFrame')
        btn_row.pack(fill=Tk.X, padx=10, pady=(0, 8))

        def _apply():
            self.rq_mode_var.set(_mode_var.get())
            self.rq_r1_vmin_var.set(_r1min.get())
            self.rq_r1_vmax_var.set(_r1max.get())
            self.rq_r2_vmin_var.set(_r2min.get())
            self.rq_r2_vmax_var.set(_r2max.get())
            self.rq_unit_var.set(_unit_var.get())
            try:
                _p = int(float(_prec_var.get()))
            except (ValueError, TypeError):
                _p = 3
            _p = max(0, min(9, _p))
            self.rq_precision_var.set(str(_p))
            mode_str = "Auto" if _mode_var.get() == 'auto' else (
                f"Manual  R1:[{_r1min.get()},{_r1max.get()}]  R2:[{_r2min.get()},{_r2max.get()}]")
            self.lbl_rq_mode.config(text=f"Mode: {mode_str}")
            popup.destroy()
            self.refresh_rq_analysis_if_visible()

        ttk.Button(btn_row, text="OK", style='Action.TButton', command=_apply).pack(
            side=Tk.RIGHT, padx=(4, 0))
        ttk.Button(btn_row, text="Cancel", style='Action.TButton',
                   command=popup.destroy).pack(side=Tk.RIGHT)

    def refresh_rq_analysis_if_visible(self):
        """Callback for the Region 1 / Region 2 / Both radio buttons. Just
        re-runs the Rq analysis with whatever data is already in self.xp/
        self.yp so the plot updates immediately when the display choice
        changes, without requiring a new sweep."""
        if self.show_rq_var.get() and self.xp and len(self.xp) >= 6:
            try:
                self.run_quench_resistance_analysis()
            except Exception as e:
                print(f"[Rq] Refresh on region-toggle failed: {e}")

    def _resolve_rq_params(self, sipm=None):
        """Return the (R_bias, N_microcell) string values to use for an Rq
        analysis. If a batch SiPM row dict is given and its per-channel
        R_bias / N (Microcells) boxes are non-empty, those override the
        main window's Forward Bias Rq settings for just that channel;
        blank boxes fall back to the main window values. Returns strings
        (not floats) so callers can drop them straight into the existing
        rq_rbias_var / rq_ncells_var StringVars."""
        main_rbias  = self.rq_rbias_var.get() if hasattr(self, 'rq_rbias_var') else "0"
        main_ncells = self.rq_ncells_var.get() if hasattr(self, 'rq_ncells_var') else "1"
        if not sipm:
            return main_rbias, main_ncells
        try:
            ch_rbias = sipm["rbias"].get().strip()
        except Exception:
            ch_rbias = ""
        try:
            ch_ncells = sipm["ncells"].get().strip()
        except Exception:
            ch_ncells = ""
        rbias  = ch_rbias if ch_rbias else main_rbias
        ncells = ch_ncells if ch_ncells else main_ncells
        return rbias, ncells

    def _apply_yscale_and_limits(self):
        """Apply the current Log/Linear scale and either auto-fit the Y
        range or pin it to the manual Y-Min/Y-Max fields. Centralized here
        so the live sweep loop, the warning/limit path, and the scale/range
        controls all behave the same way for both Positive and Negative IV.

        Note: relim()+autoscale_view() is always called first, even when
        the Y-range is manual. That call is what fits the X-axis (voltage)
        to the data as points come in -- skipping it leaves the X-axis
        stuck at its unresolved initial state and no data appears on
        screen. set_ylim() afterward just overrides the Y portion of what
        autoscale_view() picked."""
        if not (hasattr(self, 'ax') and self.ax):
            return
        current_scale = self.scale_var.get() or 'linear'
        self.ax.set_yscale(current_scale)
        self.ax.relim()
        self.ax.autoscale_view()
        if not self.auto_yscale_var.get():
            flag_min, y_min = self.is_number(self.ymin_var.get())
            flag_max, y_max = self.is_number(self.ymax_var.get())
            if flag_min and flag_max and y_max > y_min:
                if current_scale == 'log' and y_min <= 0:
                    msg.showwarning('warning', 'Y-Min must be > 0 for Log scale. Using 0.1 instead.')
                    y_min = 0.1
                self.ax.set_ylim(y_min, y_max)
            else:
                msg.showwarning('warning', 'Invalid Y-Min/Y-Max values, falling back to Auto Y-Range')
        if current_scale == 'log':
            # Final safety net against the "Data has no positive values,
            # and therefore cannot be log-scaled" crash. matplotlib's log
            # tick locator only raises later, deep inside the next
            # canvas.draw(), if the *current* view limits end up with a
            # non-positive bottom or top -- this can happen from stale
            # limits carried over from Linear scale (which commonly spans
            # negative current, e.g. Negative IV mode), from all-zero
            # data during autoscale, or from the manual Y-Min/Y-Max
            # validation above failing and leaving the old limits in
            # place untouched. Runs after both branches above so it
            # always has the final say. Force a small sane positive
            # range whenever the limits aren't both finite and positive.
            ylo, yhi = self.ax.get_ylim()
            if not (np.isfinite(ylo) and np.isfinite(yhi) and ylo > 0 and yhi > 0):
                self.ax.set_ylim(1e-3, 1)

    def apply_y_range(self):
        self.auto_yscale_var.set(False)
        self._apply_yscale_and_limits()
        self.refresh_errorbars()
        if self.figure_canvas: self.figure_canvas.draw()

    def change_scale(self):
        if hasattr(self, 'ax') and self.ax:
            current_scale = self.scale_var.get()
            # Refresh the line data too: log scale can't render negative
            # values, and a Negative IV run stores signed current in
            # self.yp, so flipping to log after the fact needs |I| re-applied
            # to the existing plot lines, not just the axis scale.
            if hasattr(self, 'yp') and self.yp:
                if current_scale == 'log':
                    y_plot = [abs(v) for v in self.yp]
                else:
                    y_plot = self.yp
                if hasattr(self, 'plot1') and self.plot1 is not None:
                    self.plot1.set_data(self.xp, y_plot)
                if hasattr(self, 'plot2') and self.plot2 is not None:
                    self.plot2.set_data(self.xp_ap, y_plot)
            self._apply_yscale_and_limits()
            self.refresh_errorbars()
            if self.figure_canvas: self.figure_canvas.draw()

    def get_errorbar_scale(self):
        """Returns the current error-bar Y-scale multiplier (self.errorbar_scale_var),
        applied to every yerr before plotting so the user can visually
        stretch or shrink whisker length without altering the underlying
        std-dev data. Falls back to 1.0 (no scaling) on bad/missing input,
        and clamps to >= 0 since a negative scale has no sensible meaning
        for an error-bar length."""
        if not hasattr(self, 'errorbar_scale_var'):
            return 1.0
        flag, scale_val = self.is_number(self.errorbar_scale_var.get())
        if not flag or scale_val < 0:
            return 1.0
        return scale_val

    def on_errorbar_scale_changed(self):
        """Callback for the Y-Scale (x) entry: re-draws whichever
        error-bar-bearing plot is currently populated with data, so the
        new scale takes effect immediately without needing a new sweep."""
        self.refresh_errorbars()
        mode = self.analysis_mode_var.get() if hasattr(self, 'analysis_mode_var') else ""
        if mode == "forward" and hasattr(self, 'show_rq_var') and self.show_rq_var.get():
            self.refresh_rq_analysis_if_visible()
        elif mode == "reverse" and hasattr(self, 'xp') and self.xp and len(self.xp) >= 5:
            try:
                self.run_breakdown_analysis()
            except Exception:
                pass

    def refresh_errorbars(self):
        """Redraw the (Measured I-V) error bars on the live Measurement
        plot from self.curr_std_arr (std-dev of the repeated readings
        taken at each point, per 'No. Meas Per Step'). Controlled by:
          - self.show_errorbars_var : show/hide toggle
          - self.errorbar_capsize_var : visual cap (whisker) length, in points
          - self.errorbar_scale_var  : multiplier on whisker length (Y-axis)
        Matplotlib error-bar artists can't be updated in place the way
        Line2D.set_data() can, so the old artist is removed and a fresh
        one is drawn each time new data comes in or a setting changes."""
        if not (hasattr(self, 'ax') and self.ax):
            return

        # Always clear the previous error-bar artist first; if the user
        # has hidden them or there's no data, we stop here with nothing drawn.
        if getattr(self, 'errbar_container', None) is not None:
            try:
                self.errbar_container.remove()
            except Exception:
                pass
            self.errbar_container = None

        if not self.show_errorbars_var.get():
            return
        if not (hasattr(self, 'xp') and self.xp and hasattr(self, 'curr_std_arr') and self.curr_std_arr):
            return
        if len(self.curr_std_arr) != len(self.xp):
            return  # arrays momentarily out of sync mid-update; skip this redraw

        try:
            flag_cap, capsize_val = self.is_number(self.errorbar_capsize_var.get())
            if not flag_cap or capsize_val < 0:
                capsize_val = 4
        except Exception:
            capsize_val = 4

        err_scale = self.get_errorbar_scale()
        y_vals = [abs(v) for v in self.yp] if self.scale_var.get() == 'log' else list(self.yp)
        yerr_scaled = [v * err_scale for v in self.curr_std_arr]

        try:
            self.errbar_container = self.ax.errorbar(
                self.xp, y_vals, yerr=yerr_scaled,
                fmt='none', ecolor='#3498DB', elinewidth=1.0,
                capsize=capsize_val, capthick=1.0, alpha=0.7, zorder=1)
        except Exception as e:
            print(f"[ErrorBars] Could not draw error bars: {e}")
            self.errbar_container = None

    # ------------------------------------------
    # 5.6 BATCH PROCESSING
    # ------------------------------------------
    # ----------------------------------------------------------------------
    # SECTION: BATCH-MODE CONFIGURATION WINDOW
    # PURPOSE: Builds the large Toplevel window used to configure a batch
    #          run: Arduino connection, per-board and per-SiPM channel
    #          tables, canvas/sequence options, and Confirm/Save/Cancel
    #          actions. This single method contains many self-contained
    #          sub-sections, each marked below.
    # ----------------------------------------------------------------------
    def open_batch_config(self):
        saved = getattr(self, "_batch_saved_state", None)

        win = Tk.Toplevel(self.window)
        win.title("Batch Mode Configuration")

        # ── Dynamically size/center the window to the actual screen ────────
        # Fixed "900x760" could be taller than a laptop's usable screen
        # height (after taskbar/menubars), which is exactly what was
        # clipping the bottom "Confirm / Save / Cancel" button row off
        # the visible window. Instead, size relative to the real screen
        # and always leave a safety margin for OS chrome.
        scr_w = win.winfo_screenwidth()
        scr_h = win.winfo_screenheight()
        win_w = min(900, int(scr_w * 0.92))
        win_h = min(760, int(scr_h * 0.88))
        win_h = max(win_h, 520)          # never go below a usable minimum
        pos_x = max(0, (scr_w - win_w) // 2)
        pos_y = max(0, (scr_h - win_h) // 2 - 20)
        win.geometry(f"{win_w}x{win_h}+{pos_x}+{pos_y}")
        win.minsize(700, 480)
        win.resizable(True, True)
        win.configure(bg="#1C2833")
        win.grab_set()           

        COLORS = {
            "bg":      "#1C2833",
            "panel":   "#2C3E50",
            "header":  "#8E44AD",
            "header2": "#1F618D",
            "accent":  "#9B59B6",
            "accent2": "#5DADE2",
            "text":    "#ECF0F1",
            "entry":   "#34495E",
            "success": "#27AE60",
            "danger":  "#E74C3C",
            "warn":    "#F39C12",
        }

        # ------------------------------------------------------------------
        # SUBSECTION: TITLE BAR
        # PURPOSE: Header strip at the top of the Batch Mode window.
        # ------------------------------------------------------------------
        title_bar = Tk.Frame(win, bg=COLORS["header"], height=40)
        title_bar.pack(fill=Tk.X)
        title_bar.pack_propagate(False)
        Tk.Label(title_bar,
                 text="⚙  Batch Mode  –  Boards (Relays) × SiPM Slots",
                 bg=COLORS["header"], fg="white",
                 font=("Segoe UI", 13, "bold")).pack(pady=8)

        # ── Bottom action bar (footer) ──────────────────────────────────
        # Built and pinned to the bottom of the window FIRST, right after
        # the title bar, using side=Tk.BOTTOM. Tkinter's pack manager
        # allocates space to widgets in the order they're packed, so
        # claiming the bottom strip now guarantees the status line and
        # the Confirm / Save / Cancel buttons always stay visible and
        # never get squeezed off-screen -- no matter how tall the
        # scrollable board/SiPM tables below end up being on a smaller
        # screen. Those scrollable sections (packed further down with
        # expand=True) simply shrink to whatever vertical space is left
        # and dynamically fit the window, with their own scrollbars
        # handling any overflow.
        # ------------------------------------------------------------------
        # SUBSECTION: FOOTER / ACTION BAR
        # PURPOSE: Bottom-pinned status line and Confirm / Save / Cancel /
        #          Reload buttons (packed first so they never get squeezed
        #          off-screen by the scrollable tables below).
        # ------------------------------------------------------------------
        footer_frame = Tk.Frame(win, bg=COLORS["bg"])
        footer_frame.pack(side=Tk.BOTTOM, fill=Tk.X)

        btn_bar = Tk.Frame(footer_frame, bg=COLORS["bg"])
        btn_bar.pack(side=Tk.BOTTOM, fill=Tk.X, padx=12, pady=6)

        # Reuse the same StringVar the main-GUI Batch Control label is bound
        # to (created in __init__) -- creating a new one here would replace
        # self._batch_status_var with a different object, silently
        # disconnecting the main-GUI label from all future .set() calls
        # (it would keep showing whatever text it had at that moment,
        # e.g. "No batch running", even while a batch is actively running).
        Tk.Label(footer_frame, textvariable=self._batch_status_var,
                 bg="#17202A", fg="#F0B27A",
                 font=("Consolas", 9), anchor="w").pack(
                     side=Tk.BOTTOM, fill=Tk.X, padx=12, pady=(2, 0))

        # ------------------------------------------------------------------
        # SUBSECTION: ARDUINO CONNECTION
        # PURPOSE: Serial port / baud rate selection for the Arduino used
        #          to switch board/SiPM relays during a batch run, plus
        #          the Terminal button below.
        # ------------------------------------------------------------------
        ard_frame = Tk.LabelFrame(win, text="  Arduino Connection  ",
                                  bg=COLORS["panel"], fg=COLORS["accent"],
                                  font=("Segoe UI", 10, "bold"))
        ard_frame.pack(fill=Tk.X, padx=12, pady=(8, 2))

        ports_found = self.batch_init_arduino()#[p.device for p in serial.tools.list_ports.comports()]
        if not ports_found:
            ports_found = ["(no port found)"]

        Tk.Label(ard_frame, text="Arduino Serial Port:", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=0, column=0, sticky="w", padx=8, pady=4)
        default_port = ports_found[0]
        if saved and saved.get("ard_port") in ports_found:
            default_port = saved["ard_port"]
        batch_ard_port_var = Tk.StringVar(value=default_port)
        ard_combo = ttk.Combobox(ard_frame, textvariable=batch_ard_port_var,
                                  values=ports_found, width=20, state="readonly")
        ard_combo.grid(row=0, column=1, padx=8, pady=4)

        Tk.Label(ard_frame, text="Baud Rate:", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=0, column=2, sticky="w", padx=8)
        default_baud = str(saved.get("baud", "9600")) if saved else "9600"
        batch_baud_var = Tk.StringVar(value=default_baud)
        ttk.Combobox(ard_frame, textvariable=batch_baud_var,
                     values=["9600", "115200", "57600"], width=10,
                     state="readonly").grid(row=0, column=3, padx=8)

        def _open_arduino_terminal():
            port = batch_ard_port_var.get()
            baud = int(batch_baud_var.get())
            # The 2-second status-light poller (_status_get_serial) keeps
            # its own connection open on this same port. A serial port can
            # only be held by one connection at a time, so if we don't
            # release the poller's connection first, opening term_ser below
            # will either fail outright or silently starve the poller for
            # as long as this Terminal window stays open -- which is why
            # the row lights get stuck gray/white and never turn green
            # while you're driving pins from here. Closing it first, and
            # pausing polling via _terminal_active, avoids the fight over
            # the port; the poller reconnects and lights resume updating
            # automatically a couple seconds after you close this window.
            self._terminal_active = True
            _status_close_serial()
            try:
                term_ser = serial.Serial(port, baud, timeout=1)
                time.sleep(1.5)
                while term_ser.in_waiting:
                    term_ser.readline()
            except Exception as e:
                self._terminal_active = False
                msg.showerror("Terminal", f"Cannot open port {port}:\n{e}",
                              parent=win)
                return

            # ----------------------------------------------------------
            # SUBSECTION: ARDUINO TERMINAL POPUP
            # PURPOSE: Raw serial console for sending manual commands to
            #          the Arduino while the Batch Mode window is open.
            # ----------------------------------------------------------
            term_win = Tk.Toplevel(win)
            term_win.title(f"Arduino Terminal  –  {port}  @  {baud}")
            term_win.geometry("640x420")
            term_win.configure(bg="#0D0D0D")
            term_win.grab_set()

            log_frame = Tk.Frame(term_win, bg="#0D0D0D")
            log_frame.pack(fill=Tk.BOTH, expand=True, padx=6, pady=(6, 2))

            log_sb = ttk.Scrollbar(log_frame)
            log_sb.pack(side=Tk.RIGHT, fill=Tk.Y)

            log_text = Tk.Text(
                log_frame, bg="#0D0D0D", fg="#00FF41",
                font=("Consolas", 10), wrap="word",
                yscrollcommand=log_sb.set,
                state="disabled"
            )
            log_text.pack(fill=Tk.BOTH, expand=True)
            log_sb.config(command=log_text.yview)

            def _log(line, color="#00FF41"):
                log_text.config(state="normal")
                log_text.insert(Tk.END, line + "\n", color)
                log_text.tag_config(color, foreground=color)
                log_text.see(Tk.END)
                log_text.config(state="disabled")

            _log(f"[Terminal] Connected to {port} @ {baud} baud", "#AAAAAA")
            _log(f"[Terminal] Example: Board 9 Update 12 1", "#AAAAAA")

            input_frame = Tk.Frame(term_win, bg="#1A1A1A")
            input_frame.pack(fill=Tk.X, padx=6, pady=(0, 4))

            cmd_var = Tk.StringVar()
            cmd_entry = Tk.Entry(
                input_frame, textvariable=cmd_var,
                bg="#1A1A1A", fg="#FFFF00",
                insertbackground="white",
                font=("Consolas", 11), relief="flat"
            )
            cmd_entry.pack(side=Tk.LEFT, fill=Tk.X, expand=True, padx=(4, 2), pady=4)
            cmd_entry.focus_set()

            def _send_cmd(event=None):
                raw = cmd_var.get().strip()
                if not raw:
                    return
                cmd_var.set("")
                _log(f">>> {raw}", "#FFFF00")
                try:
                    term_ser.write((raw + "\n").encode("utf-8"))
                    time.sleep(0.2)
                    reply_lines = []
                    while term_ser.in_waiting:
                        reply_lines.append(
                            term_ser.readline().decode("utf-8",
                                                       errors="replace").rstrip())
                    if reply_lines:
                        for rl in reply_lines:
                            _log(f"<<< {rl}", "#00FF41")
                    else:
                        _log("<<< (no reply)", "#888888")
                except Exception as exc:
                    _log(f"[ERROR] {exc}", "#FF4444")

            cmd_entry.bind("<Return>", _send_cmd)

            send_btn = Tk.Button(
                input_frame, text="Send ↵",
                bg=COLORS["success"], fg="white",
                font=("Segoe UI", 10, "bold"), relief="flat", padx=8,
                command=_send_cmd
            )
            send_btn.pack(side=Tk.LEFT, pady=4)

            def _poll():
                if not term_win.winfo_exists():
                    return
                try:
                    while term_ser.in_waiting:
                        line = term_ser.readline().decode(
                            "utf-8", errors="replace").rstrip()
                        if line:
                            _log(f"[RX] {line}", "#44DDFF")
                except Exception:
                    pass
                term_win.after(200, _poll)

            _poll()

            def _on_close():
                try:
                    term_ser.close()
                except Exception:
                    pass
                self._terminal_active = False
                term_win.destroy()

            term_win.protocol("WM_DELETE_WINDOW", _on_close)

        Tk.Button(ard_frame, text=" Arduino Terminal",
                  bg=COLORS["accent"], fg="white",
                  font=("Segoe UI", 10, "bold"), relief="flat", padx=8,
                  command=_open_arduino_terminal).grid(
                      row=0, column=4, padx=12, pady=4)

        # ── Live channel-status indicator infrastructure ────────────────────
        # A small light next to each row's "Del" button shows whether that
        # board/SiPM relay is currently ON (green), OFF (red), or unknown /
        # not reachable (gray). Every row's actual hardware pin state is
        # polled every 2 s for as long as this Batch Mode window stays
        # open, regardless of the "Use in I-V" checkbox -- that checkbox
        # only controls batch participation, not whether we bother
        # reading the relay's real state. Polling uses its own serial
        # connection (separate from the one a live batch run opens) and
        # simply reports "unknown" while a batch is actively running,
        # rather than fight it for the port.
        _status_ser_holder = {"ser": None}

        def _status_get_serial():
            ser = _status_ser_holder["ser"]
            if ser is not None:
                try:
                    if ser.is_open:
                        return ser
                except Exception:
                    pass
            try:
                ser = serial.Serial(batch_ard_port_var.get(),
                                     int(batch_baud_var.get()), timeout=0.5)
                time.sleep(0.3)
                while ser.in_waiting:
                    ser.readline()
                _status_ser_holder["ser"] = ser
                return ser
            except Exception:
                _status_ser_holder["ser"] = None
                return None

        def _status_close_serial():
            ser = _status_ser_holder["ser"]
            if ser is not None:
                try:
                    ser.close()
                except Exception:
                    pass
                _status_ser_holder["ser"] = None

        def _status_query_pin(slave_id, pin):
            """Ask the Arduino whether a given relay pin is currently ON.
            Returns True / False / None (unknown -- no connection, no
            reply, a live batch run is using the port already, or the
            firmware reported the pin as out of range / errored).

            Matches Master.ino's Read/R reply formats exactly:
              "Pin <pin> State: 1"            -- local (Board 0) ON
              "Pin <pin> State: 0"            -- local (Board 0) OFF
              "Pin <pin> State: OUT_OF_RANGE" -- local, bad pin number
              "Pin State (Board <id>): 1"     -- remote board ON
              "Pin State (Board <id>): 0"     -- remote board OFF
              "Error: Board <id> ..."         -- remote board unreachable
            Scanning only after the literal "State:" token (rather than
            the last 0/1 digit anywhere in the reply) avoids misreading
            a board address digit inside an error line as a pin state.
            """
            if self.batch_mode_active or getattr(self, "_terminal_active", False):
                return None
            ser = _status_get_serial()
            if ser is None:
                return None
            try:
                cmd = f"Board {slave_id} Read {pin}"
                ser.write((cmd.strip() + "\n").encode("utf-8"))
                time.sleep(0.12)
                reply = ""
                while ser.in_waiting:
                    reply += ser.readline().decode("utf-8", errors="replace")
                reply = reply.strip()
                if not reply:
                    return None
                # Remote-board reply: "Pin State (Board <id>): 1"
                m = re.findall(r"State\s*\(Board\s*\d+\):\s*([01])\b", reply)
                if m:
                    return m[-1] == "1"
                # Local-board (slave_id 0) reply: "Pin <pin> State: 1"
                # -- no "(Board <id>)" segment, so it needs its own match.
                m = re.findall(r"State:\s*([01])\b", reply)
                if m:
                    return m[-1] == "1"
                return None
            except Exception:
                _status_close_serial()
                return None

        def _status_set_pin(slave_id, pin, turn_on):
            """Send an Update command switching a relay pin ON/OFF via
            the status-polling serial connection (mirrors the
            'Board <id> Update <pin> <0|1>' command the live batch run
            itself sends -- see cmd_board_on/cmd_sipm_on below).
            Returns True if the Arduino replied at all (command
            accepted), False if it replied with an explicit error, or
            None if unreachable (no connection, no reply, or a live
            batch run is already using the port).
            """
            if self.batch_mode_active or getattr(self, "_terminal_active", False):
                return None
            ser = _status_get_serial()
            if ser is None:
                return None
            try:
                cmd = f"Board {slave_id} Update {pin} {1 if turn_on else 0}"
                ser.write((cmd.strip() + "\n").encode("utf-8"))
                time.sleep(0.15)
                reply = ""
                while ser.in_waiting:
                    reply += ser.readline().decode("utf-8", errors="replace")
                reply = reply.strip()
                if not reply:
                    return None
                if reply.lower().startswith("error"):
                    return False
                return True
            except Exception:
                _status_close_serial()
                return None

        def _style_toggle_btn(btn, var, connected=True):
            """Paint a channel's ON/OFF button to match its current
            (commanded) enabled state -- but only while the channel is
            actually reachable. When the most recent status check could
            not reach this board/pin (state=None, same condition that
            paints the status light gray), the button must follow the
            light: show OFF in gray and stop inviting clicks, since a
            command sent right now has nowhere to go anyway."""
            try:
                if not connected:
                    btn.config(text="OFF", bg="#7F8C8D", fg="white",
                               state=Tk.DISABLED)
                elif var.get():
                    btn.config(text="ON", bg=COLORS["success"], fg="white",
                               state=Tk.NORMAL)
                else:
                    btn.config(text="OFF", bg=COLORS["danger"], fg="white",
                               state=Tk.NORMAL)
            except Exception:
                pass

        def _set_light(row_data, state):
            """state: True=ON (green), False=OFF (red), None=unknown (gray)."""
            canvas = row_data.get("light_canvas")
            oval = row_data.get("light_oval")
            if canvas is None or oval is None:
                return
            if state is True:
                color = "#2ECC71"
            elif state is False:
                color = "#E74C3C"
            else:
                color = "#7F8C8D"
            try:
                canvas.itemconfig(oval, fill=color, outline=color)
            except Exception:
                pass

        def _make_light_widget(parent):
            cv = Tk.Canvas(parent, width=16, height=16,
                            bg=COLORS["bg"], highlightthickness=0)
            oval = cv.create_oval(3, 3, 13, 13, fill="#7F8C8D", outline="#7F8C8D")
            return cv, oval

        self._batch_status_ser_holder = _status_ser_holder

        self._batch_board_rows = []
        self._batch_sipm_rows  = []

        # ------------------------------------------------------------------
        # SUBSECTION: BOARD TABLE CONTROLS
        # PURPOSE: "Number of boards" spinbox and the force-check-all-now
        #          control for the board channel table below.
        # ------------------------------------------------------------------
        board_ctrl_frame = Tk.Frame(win, bg=COLORS["bg"])
        board_ctrl_frame.pack(fill=Tk.X, padx=12, pady=(8, 2))

        Tk.Label(board_ctrl_frame, text="① Boards / Relay Switches:",
                 bg=COLORS["bg"], fg=COLORS["text"],
                 font=("Segoe UI", 10, "bold")).pack(side=Tk.LEFT)

        saved_boards = (saved or {}).get("boards", [])
        n_boards_var = Tk.IntVar(value=len(saved_boards) if saved_boards else 2)
        n_boards_spin = Tk.Spinbox(board_ctrl_frame, from_=1, to=16,
                                    textvariable=n_boards_var, width=4,
                                    font=("Segoe UI", 11), bg=COLORS["entry"],
                                    fg=COLORS["text"],
                                    buttonbackground=COLORS["panel"])
        n_boards_spin.pack(side=Tk.LEFT, padx=8)

        Tk.Label(board_ctrl_frame,
                 text="(each board = one relay on a digital pin; "
                      "turned ON before its SiPM slots, OFF after)",
                 bg=COLORS["bg"], fg="#95A5A6", font=("Segoe UI", 9)).pack(
                     side=Tk.LEFT, padx=4)

        force_check_frame = Tk.Frame(win, bg=COLORS["bg"])
        force_check_frame.pack(fill=Tk.X, padx=12, pady=(0, 4))

        # ------------------------------------------------------------------
        # SUBSECTION: BOARD CHANNEL TABLE
        # PURPOSE: Scrollable table of board rows (label, slave ID, pin,
        #          enabled/participate flags, live ON/OFF status light,
        #          toggle and delete buttons) used to drive relays during
        #          a batch run.
        # ------------------------------------------------------------------
        board_table_outer = Tk.Frame(win, bg=COLORS["bg"])
        board_table_outer.pack(fill=Tk.X, padx=12, pady=2)

        # Row/header pixel sizes used to size the table to its actual
        # content (clamped between a usable minimum and a max) instead of
        # greedily expanding to fill the window -- that greedy expansion
        # was what starved the "Divided Canvas Settings" / "Sequence
        # Options" panels and the bottom buttons of vertical room.
        _ROW_H, _HDR_H = 28, 26
        _BOARD_MIN_H, _BOARD_MAX_H = 76, 170

        board_canvas = Tk.Canvas(board_table_outer, bg=COLORS["bg"],
                                  highlightthickness=0, height=_BOARD_MIN_H)
        board_scroll = ttk.Scrollbar(board_table_outer, orient="vertical",
                                      command=board_canvas.yview)
        board_canvas.configure(yscrollcommand=board_scroll.set)
        board_scroll.pack(side=Tk.RIGHT, fill=Tk.Y)
        board_canvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True)

        board_inner = Tk.Frame(board_canvas, bg=COLORS["bg"])
        board_win_id = board_canvas.create_window((0, 0), window=board_inner,
                                                    anchor="nw")
        board_inner.bind("<Configure>",
                          lambda e: board_canvas.configure(
                              scrollregion=board_canvas.bbox("all")))
        board_canvas.bind("<Configure>",
                          lambda e: board_canvas.itemconfig(
                              board_win_id, width=e.width))

        board_hdr_cols  = ["#", "Board Label", "Slave ID", "Digital Pin",
                           "Use in I-V", "On/Off", "Status", "Del"]
        board_hdr_widths = [3, 22, 10, 10, 12, 6, 6, 4]
        for col, (h, w) in enumerate(zip(board_hdr_cols, board_hdr_widths)):
            Tk.Label(board_inner, text=h, bg=COLORS["header"], fg="white",
                     font=("Segoe UI", 9, "bold"), width=w,
                     anchor="center", relief="flat").grid(
                         row=0, column=col, padx=1, pady=2, sticky="ew")
        # Let the "Board Label" column (index 1) soak up any extra
        # horizontal room when the window is widened; the narrow
        # numeric/checkbox/status columns stay at their natural size
        # instead of stretching into awkward gaps.
        for col in range(len(board_hdr_cols)):
            board_inner.grid_columnconfigure(col, weight=1 if col == 1 else 0)

        def _resize_board_canvas():
            n = len(self._batch_board_rows)
            h = min(_BOARD_MAX_H,
                    max(_BOARD_MIN_H, _HDR_H + _ROW_H * n))
            board_canvas.config(height=h)

        def _regrid_board_rows():
            for i, row_data in enumerate(self._batch_board_rows):
                r = i + 1
                row_data["num_var"].set(str(i + 1))
                for col, w in enumerate(row_data["widgets"]):
                    w.grid(row=r, column=col, padx=1, pady=1, sticky="ew")

        def _delete_board_row(row_data):
            if len(self._batch_board_rows) <= 1:
                msg.showwarning("Batch Mode",
                                 "At least one board is required.", parent=win)
                return
            for w in row_data["widgets"]:
                w.destroy()
            self._batch_board_rows.remove(row_data)
            _regrid_board_rows()
            _resize_board_canvas()
            n_boards_var.set(len(self._batch_board_rows))
            _refresh_board_choices()

        def _check_board_row_now(row_data):
            if not win.winfo_exists():
                return
            try:
                # Always ask the Arduino for the real pin state -- this
                # must NOT depend on the commanded on/off var, otherwise
                # a channel that's genuinely OFF (but still reachable)
                # incorrectly shows gray ("unknown") instead of red.
                state = _status_query_pin(row_data["slave_id"].get(),
                                           row_data["pin"].get())
            except (Tk.TclError, Exception):
                # A momentarily empty/invalid Spinbox field (e.g. mid-edit)
                # must not be allowed to propagate -- that would kill the
                # recurring win.after() poll loop for every row, not just
                # this one. Just report "unknown" for this pass instead.
                state = None
            _set_light(row_data, state)
            # The very first time we get a real reading for this row
            # (e.g. right after the window opens / the port is set),
            # sync the commanded ON/OFF button to match the actual
            # hardware state -- otherwise the button can show "ON"
            # (green) while the status light shows red, which looks
            # like a bug even though the two are tracking different
            # things (commanded vs. actual).
            if not row_data.get("_synced_once") and state is not None:
                row_data["_synced_once"] = True
                if row_data["enabled"].get() != state:
                    row_data["enabled"].set(state)
                    return  # the trace on "enabled" re-invokes this check
            btn = row_data.get("toggle_btn")
            if btn is not None:
                _style_toggle_btn(btn, row_data["enabled"],
                                   connected=(state is not None))

        def _build_board_rows(n):
            while len(self._batch_board_rows) > n:
                row_data = self._batch_board_rows.pop()
                for w in row_data["widgets"]:
                    w.destroy()
            while len(self._batch_board_rows) < n:
                idx = len(self._batch_board_rows)
                r = idx + 1

                if idx < len(saved_boards):
                    sb = saved_boards[idx]
                    init_label    = sb.get("label", f"Board_{idx+1}")
                    init_slave_id = sb.get("slave_id", 8)
                    init_pin      = sb.get("pin", idx + 2)
                    init_enabled  = sb.get("enabled", True)
                    # Older save files only had "enabled", which used to
                    # double as the batch-participation flag -- fall back
                    # to it so old configs still restore correctly.
                    init_participate = sb.get("participate", init_enabled)
                else:
                    init_label    = f"Board_{idx+1}"
                    init_slave_id = 8
                    init_pin      = idx + 2
                    init_enabled  = True
                    init_participate = True

                num_var      = Tk.StringVar(value=str(idx + 1))
                label_var    = Tk.StringVar(value=init_label)
                slave_id_var = Tk.IntVar(value=init_slave_id)
                pin_var      = Tk.IntVar(value=init_pin)
                enabled_var  = Tk.BooleanVar(value=init_enabled)
                participate_var = Tk.BooleanVar(value=init_participate)

                num_lbl = Tk.Label(board_inner, textvariable=num_var,
                                    bg=COLORS["panel"], fg=COLORS["text"],
                                    font=("Segoe UI", 9, "bold"), width=3,
                                    anchor="center")
                name_ent = Tk.Entry(board_inner, textvariable=label_var,
                                     bg=COLORS["entry"], fg=COLORS["text"],
                                     insertbackground="white", width=22,
                                     font=("Segoe UI", 9))
                slave_spin = Tk.Spinbox(board_inner, from_=0, to=127,
                                         textvariable=slave_id_var, width=10,
                                         bg=COLORS["entry"], fg=COLORS["text"],
                                         buttonbackground=COLORS["panel"],
                                         font=("Consolas", 9),
                                         justify="center")
                pin_spin = Tk.Spinbox(board_inner, from_=0, to=53,
                                       textvariable=pin_var, width=10,
                                       bg=COLORS["entry"], fg=COLORS["text"],
                                       buttonbackground=COLORS["panel"],
                                       font=("Consolas", 9),
                                       justify="center")
                def _toggle_board_channel(sid_v=slave_id_var, pin_v=pin_var,
                                           ena_v=enabled_var):
                    try:
                        sid = int(sid_v.get())
                        pin = int(pin_v.get())
                    except (ValueError, Tk.TclError):
                        msg.showwarning("Batch Mode",
                                         "Invalid Slave ID / Pin.", parent=win)
                        return
                    turn_on = not ena_v.get()
                    result = _status_set_pin(sid, pin, turn_on)
                    if result is None:
                        msg.showwarning(
                            "Batch Mode",
                            "Could not reach the Arduino to switch "
                            "this channel.", parent=win)
                        return
                    ena_v.set(turn_on)

                participate_chk = Tk.Checkbutton(
                    board_inner, variable=participate_var,
                    bg=COLORS["panel"], activebackground=COLORS["panel"],
                    selectcolor=COLORS["entry"], highlightthickness=0,
                    bd=0)

                toggle_btn = Tk.Button(board_inner, font=("Segoe UI", 8, "bold"),
                                        relief="flat", width=5, padx=2, pady=0,
                                        command=_toggle_board_channel)
                _style_toggle_btn(toggle_btn, enabled_var)
                light_cv, light_oval = _make_light_widget(board_inner)

                row_data = {
                    "num_var":     num_var,
                    "label":       label_var,
                    "slave_id":    slave_id_var,
                    "pin":         pin_var,
                    "participate": participate_var,
                    "enabled":     enabled_var,
                    "toggle_btn":  toggle_btn,
                    "light_canvas": light_cv,
                    "light_oval":   light_oval,
                    "_synced_once": False,
                }

                del_btn = Tk.Button(board_inner, text="✕",
                                     bg=COLORS["danger"], fg="white",
                                     font=("Segoe UI", 8, "bold"),
                                     relief="flat", padx=2, pady=0,
                                     command=lambda rd=row_data: _delete_board_row(rd))

                widgets = [num_lbl, name_ent, slave_spin, pin_spin,
                           participate_chk, toggle_btn, light_cv, del_btn]
                for col, w in enumerate(widgets):
                    w.grid(row=r, column=col, padx=1, pady=1, sticky="ew")

                row_data["widgets"] = widgets
                label_var.trace_add("write", lambda *_: _refresh_board_choices())
                enabled_var.trace_add(
                    "write", lambda *_, rd=row_data: _check_board_row_now(rd))
                self._batch_board_rows.append(row_data)
                _check_board_row_now(row_data)

        _build_board_rows(n_boards_var.get())
        _resize_board_canvas()

        def _on_n_boards_change(*_):
            try:
                n = int(n_boards_var.get())
                if 1 <= n <= 16:
                    _build_board_rows(n)
                    _resize_board_canvas()
                    _refresh_board_choices()
            except (ValueError, Tk.TclError):
                pass

        n_boards_var.trace_add("write", _on_n_boards_change)
        n_boards_spin.config(command=_on_n_boards_change)

        # ------------------------------------------------------------------
        # SUBSECTION: SiPM TABLE CONTROLS
        # PURPOSE: "Number of SiPMs" spinbox controlling the SiPM channel
        #          table below.
        # ------------------------------------------------------------------
        sipm_ctrl_frame = Tk.Frame(win, bg=COLORS["bg"])
        sipm_ctrl_frame.pack(fill=Tk.X, padx=12, pady=(10, 2))

        Tk.Label(sipm_ctrl_frame, text="② SiPM Slots (I-V measurements):",
                 bg=COLORS["bg"], fg=COLORS["text"],
                 font=("Segoe UI", 10, "bold")).pack(side=Tk.LEFT)

        Tk.Label(sipm_ctrl_frame,
                 text="  (leave R_bias / N blank to use the main window's "
                      "Forward Bias Rq settings for that channel)",
                 bg=COLORS["bg"], fg="#7F8C8D",
                 font=("Segoe UI", 8, "italic")).pack(side=Tk.LEFT)

        saved_sipms = (saved or {}).get("sipms", [])
        n_sipms_var = Tk.IntVar(value=len(saved_sipms) if saved_sipms else 2)
        n_sipms_spin = Tk.Spinbox(sipm_ctrl_frame, from_=1, to=64,
                                   textvariable=n_sipms_var, width=4,
                                   font=("Segoe UI", 11), bg=COLORS["entry"],
                                   fg=COLORS["text"],
                                   buttonbackground=COLORS["panel"])
        n_sipms_spin.pack(side=Tk.LEFT, padx=8)

        Tk.Label(sipm_ctrl_frame,
                 text="(each SiPM slot runs its own I-V + fit, while its "
                      "assigned board's relay is ON)",
                 bg=COLORS["bg"], fg="#95A5A6", font=("Segoe UI", 9)).pack(
                     side=Tk.LEFT, padx=4)

        # ------------------------------------------------------------------
        # SUBSECTION: SiPM CHANNEL TABLE
        # PURPOSE: Scrollable table of per-SiPM rows (label, board, pin,
        #          enabled/participate flags, per-channel Rq overrides,
        #          live status light) used during a batch run.
        # ------------------------------------------------------------------
        sipm_table_outer = Tk.Frame(win, bg=COLORS["bg"])
        sipm_table_outer.pack(fill=Tk.BOTH, expand=True, padx=12, pady=2)

        sipm_canvas = Tk.Canvas(sipm_table_outer, bg=COLORS["bg"],
                                 highlightthickness=0)
        sipm_scroll = ttk.Scrollbar(sipm_table_outer, orient="vertical",
                                     command=sipm_canvas.yview)
        sipm_canvas.configure(yscrollcommand=sipm_scroll.set)
        sipm_scroll.pack(side=Tk.RIGHT, fill=Tk.Y)
        sipm_canvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True)

        sipm_inner = Tk.Frame(sipm_canvas, bg=COLORS["bg"])
        sipm_win_id = sipm_canvas.create_window((0, 0), window=sipm_inner,
                                                  anchor="nw")
        sipm_inner.bind("<Configure>",
                        lambda e: sipm_canvas.configure(
                            scrollregion=sipm_canvas.bbox("all")))
        sipm_canvas.bind("<Configure>",
                         lambda e: sipm_canvas.itemconfig(
                             sipm_win_id, width=e.width))

        sipm_hdr_cols   = ["#", "SiPM Label", "Assign to Board",
                           "Digital Pin", "R_bias (\u03A9)", "N (\u03BCcell)",
                           "Use in I-V", "On/Off", "Status", "Del"]
        sipm_hdr_widths = [3, 22, 18, 10, 10, 8, 12, 6, 6, 4]
        for col, (h, w) in enumerate(zip(sipm_hdr_cols, sipm_hdr_widths)):
            Tk.Label(sipm_inner, text=h, bg=COLORS["header2"], fg="white",
                     font=("Segoe UI", 9, "bold"), width=w,
                     anchor="center", relief="flat").grid(
                         row=0, column=col, padx=1, pady=2, sticky="ew")
        # Same idea as the board table: give the "SiPM Label" column the
        # extra horizontal space on resize, everything else stays put.
        for col in range(len(sipm_hdr_cols)):
            sipm_inner.grid_columnconfigure(col, weight=1 if col == 1 else 0)

        def _refresh_board_choices(*_):
            labels = [r["label"].get() for r in self._batch_board_rows]
            for r in self._batch_sipm_rows:
                combo = r["board_combo"]
                combo["values"] = labels
                if r["board_var"].get() not in labels:
                    r["board_var"].set(labels[0] if labels else "")

        def _regrid_sipm_rows():
            for i, row_data in enumerate(self._batch_sipm_rows):
                r = i + 1
                row_data["num_var"].set(str(i + 1))
                for col, w in enumerate(row_data["widgets"]):
                    w.grid(row=r, column=col, padx=1, pady=1, sticky="ew")

        def _delete_sipm_row(row_data):
            if len(self._batch_sipm_rows) <= 1:
                msg.showwarning("Batch Mode",
                                 "At least one SiPM slot is required.", parent=win)
                return
            for w in row_data["widgets"]:
                w.destroy()
            self._batch_sipm_rows.remove(row_data)
            _regrid_sipm_rows()
            n_sipms_var.set(len(self._batch_sipm_rows))

        def _check_sipm_row_now(row_data):
            if not win.winfo_exists():
                return
            try:
                # Always ask the Arduino for the real pin state -- see
                # the comment in _check_board_row_now for why this must
                # not be gated on the commanded on/off var.
                board_label = row_data["board_var"].get().strip()
                slave_id = None
                for br in self._batch_board_rows:
                    if br["label"].get().strip() == board_label:
                        slave_id = br["slave_id"].get()
                        break
                if slave_id is not None:
                    state = _status_query_pin(slave_id, row_data["pin"].get())
                else:
                    state = None
            except (Tk.TclError, Exception):
                # Same reasoning as _check_board_row_now -- never let a
                # bad/empty field kill the recurring poll loop.
                state = None
            _set_light(row_data, state)
            # See the matching comment in _check_board_row_now: sync the
            # commanded ON/OFF button to the real state the first time we
            # get one, so it doesn't show "ON" next to a red light.
            if not row_data.get("_synced_once") and state is not None:
                row_data["_synced_once"] = True
                if row_data["enabled"].get() != state:
                    row_data["enabled"].set(state)
                    return  # the trace on "enabled" re-invokes this check
            btn = row_data.get("toggle_btn")
            if btn is not None:
                _style_toggle_btn(btn, row_data["enabled"],
                                   connected=(state is not None))

        def _build_sipm_rows(n):
            board_labels = [r["label"].get() for r in self._batch_board_rows]
            while len(self._batch_sipm_rows) > n:
                row_data = self._batch_sipm_rows.pop()
                for w in row_data["widgets"]:
                    w.destroy()
            while len(self._batch_sipm_rows) < n:
                idx = len(self._batch_sipm_rows)
                r = idx + 1

                if idx < len(saved_sipms):
                    ss = saved_sipms[idx]
                    init_label   = ss.get("label", f"SiPM_{idx+1}")
                    init_board   = ss.get("board_label", "")
                    init_pin     = ss.get("pin", idx + 2)
                    init_enabled = ss.get("enabled", True)
                    # Older save files only had "enabled", which used to
                    # double as the batch-participation flag -- fall back
                    # to it so old configs still restore correctly.
                    init_participate = ss.get("participate", init_enabled)
                    # Per-channel Rq correction overrides. Left blank
                    # (default) means "fall back to the main-window
                    # Forward Bias Rq settings" at analysis time -- see
                    # _resolve_rq_params() in run_batch_sequence.
                    init_rbias   = ss.get("rbias", "")
                    init_ncells  = ss.get("ncells", "")
                else:
                    init_label   = f"SiPM_{idx+1}"
                    init_board   = ""
                    init_pin     = idx + 2
                    init_enabled = True
                    init_participate = True
                    init_rbias   = ""
                    init_ncells  = ""

                if init_board not in board_labels:
                    init_board = (board_labels[idx % len(board_labels)]
                                   if board_labels else "")

                num_var     = Tk.StringVar(value=str(idx + 1))
                label_var   = Tk.StringVar(value=init_label)
                board_var   = Tk.StringVar(value=init_board)
                pin_var     = Tk.IntVar(value=init_pin)
                enabled_var = Tk.BooleanVar(value=init_enabled)
                participate_var = Tk.BooleanVar(value=init_participate)
                rbias_var   = Tk.StringVar(value=init_rbias)
                ncells_var  = Tk.StringVar(value=init_ncells)

                num_lbl = Tk.Label(sipm_inner, textvariable=num_var,
                                    bg=COLORS["panel"], fg=COLORS["text"],
                                    font=("Segoe UI", 9, "bold"), width=3,
                                    anchor="center")
                name_ent = Tk.Entry(sipm_inner, textvariable=label_var,
                                     bg=COLORS["entry"], fg=COLORS["text"],
                                     insertbackground="white", width=22,
                                     font=("Segoe UI", 9))
                board_combo = ttk.Combobox(sipm_inner, textvariable=board_var,
                                            values=board_labels, width=16,
                                            state="readonly")
                pin_spin = Tk.Spinbox(sipm_inner, from_=0, to=53,
                                       textvariable=pin_var, width=10,
                                       bg=COLORS["entry"], fg=COLORS["text"],
                                       buttonbackground=COLORS["panel"],
                                       font=("Consolas", 9),
                                       justify="center")

                # Per-channel Rq correction parameters. Left blank, the
                # per-SiPM analysis in run_batch_sequence falls back to
                # the main window's Forward Bias R_bias / N (Microcells)
                # values -- see _resolve_rq_params().
                rbias_ent = Tk.Entry(sipm_inner, textvariable=rbias_var,
                                      bg=COLORS["entry"], fg=COLORS["text"],
                                      insertbackground="white", width=10,
                                      font=("Consolas", 9), justify="center")
                ncells_ent = Tk.Entry(sipm_inner, textvariable=ncells_var,
                                       bg=COLORS["entry"], fg=COLORS["text"],
                                       insertbackground="white", width=8,
                                       font=("Consolas", 9), justify="center")
                def _toggle_sipm_channel(pin_v=pin_var, board_v=board_var,
                                          ena_v=enabled_var):
                    board_label = board_v.get().strip()
                    slave_id = None
                    for br in self._batch_board_rows:
                        if br["label"].get().strip() == board_label:
                            slave_id = br["slave_id"].get()
                            break
                    if slave_id is None:
                        msg.showwarning(
                            "Batch Mode",
                            "Assign this SiPM slot to a board first.",
                            parent=win)
                        return
                    try:
                        sid = int(slave_id)
                        pin = int(pin_v.get())
                    except (ValueError, Tk.TclError):
                        msg.showwarning("Batch Mode",
                                         "Invalid Slave ID / Pin.", parent=win)
                        return
                    turn_on = not ena_v.get()
                    result = _status_set_pin(sid, pin, turn_on)
                    if result is None:
                        msg.showwarning(
                            "Batch Mode",
                            "Could not reach the Arduino to switch "
                            "this channel.", parent=win)
                        return
                    ena_v.set(turn_on)

                participate_chk = Tk.Checkbutton(
                    sipm_inner, variable=participate_var,
                    bg=COLORS["panel"], activebackground=COLORS["panel"],
                    selectcolor=COLORS["entry"], highlightthickness=0,
                    bd=0)

                toggle_btn = Tk.Button(sipm_inner, font=("Segoe UI", 8, "bold"),
                                        relief="flat", width=5, padx=2, pady=0,
                                        command=_toggle_sipm_channel)
                _style_toggle_btn(toggle_btn, enabled_var)
                light_cv, light_oval = _make_light_widget(sipm_inner)

                row_data = {
                    "num_var":     num_var,
                    "label":       label_var,
                    "board_var":   board_var,
                    "board_combo": board_combo,
                    "pin":         pin_var,
                    "rbias":       rbias_var,
                    "ncells":      ncells_var,
                    "participate": participate_var,
                    "enabled":     enabled_var,
                    "toggle_btn":  toggle_btn,
                    "light_canvas": light_cv,
                    "light_oval":   light_oval,
                    "_synced_once": False,
                }

                del_btn = Tk.Button(sipm_inner, text="✕",
                                     bg=COLORS["danger"], fg="white",
                                     font=("Segoe UI", 8, "bold"),
                                     relief="flat", padx=2, pady=0,
                                     command=lambda rd=row_data: _delete_sipm_row(rd))

                widgets = [num_lbl, name_ent, board_combo, pin_spin,
                           rbias_ent, ncells_ent,
                           participate_chk, toggle_btn, light_cv, del_btn]
                for col, w in enumerate(widgets):
                    w.grid(row=r, column=col, padx=1, pady=1, sticky="ew")

                row_data["widgets"] = widgets
                enabled_var.trace_add(
                    "write", lambda *_, rd=row_data: _check_sipm_row_now(rd))
                board_var.trace_add(
                    "write", lambda *_, rd=row_data: _check_sipm_row_now(rd))
                self._batch_sipm_rows.append(row_data)
                _check_sipm_row_now(row_data)

        _build_sipm_rows(n_sipms_var.get())

        def _on_n_sipms_change(*_):
            try:
                n = int(n_sipms_var.get())
                if 1 <= n <= 64:
                    _build_sipm_rows(n)
            except (ValueError, Tk.TclError):
                pass

        n_sipms_var.trace_add("write", _on_n_sipms_change)
        n_sipms_spin.config(command=_on_n_sipms_change)

        # ── Continuous 2-second status polling while this window is open ───
        # Re-checks every enabled board/SiPM row's light on a 2 s cadence.
        # Stops on its own once the Batch Mode window is closed (the
        # win.winfo_exists() guard breaks the win.after() recursion).
        # ------------------------------------------------------------------
        # SUBSECTION: LIVE CHANNEL-STATUS POLLING
        # PURPOSE: Recurring 2-second poll that refreshes every board/SiPM
        #          row status light while this window stays open.
        # ------------------------------------------------------------------
        def _poll_all_channel_lights():
            if not win.winfo_exists():
                return
            try:
                for rd in list(self._batch_board_rows):
                    _check_board_row_now(rd)
                for rd in list(self._batch_sipm_rows):
                    _check_sipm_row_now(rd)
            except Exception:
                # Never let an unexpected error here break the recurring
                # after() chain below -- that would silently freeze every
                # row's light gray/unlit for the rest of the session.
                pass
            finally:
                win.after(2000, _poll_all_channel_lights)

        win.after(2000, _poll_all_channel_lights)

        def _force_check_all_now():
            """Manually triggered immediate re-check of every enabled
            board/SiPM row's light, independent of the 2 s auto-poll
            cadence. Also nudges the poller to reconnect if its serial
            connection had dropped."""
            if not win.winfo_exists():
                return
            _status_close_serial()
            try:
                for rd in list(self._batch_board_rows):
                    _check_board_row_now(rd)
                for rd in list(self._batch_sipm_rows):
                    _check_sipm_row_now(rd)
            except Exception:
                pass

        force_check_btn = Tk.Button(
            force_check_frame, text="🔄  Force Check All Channels Now",
            bg=COLORS["accent"], fg="white",
            font=("Segoe UI", 9, "bold"), relief="flat", padx=10, pady=4,
            command=_force_check_all_now)
        force_check_btn.pack(side=Tk.LEFT)

        # ------------------------------------------------------------------
        # SUBSECTION: DIVIDED CANVAS SETTINGS
        # PURPOSE: Configure how the plot canvas is subdivided for
        #          multi-channel batch plotting.
        # ------------------------------------------------------------------
        # ── Canvas Grid Options ────────────────────────────────────────────────
        grid_frame = Tk.LabelFrame(win, text="  Divided Canvas Settings  ",
                                   bg=COLORS["panel"], fg=COLORS["accent"],
                                   font=("Segoe UI", 10, "bold"))
        grid_frame.pack(fill=Tk.X, padx=12, pady=(8, 2))

        Tk.Label(grid_frame, text="Grid Rows (n):", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=0, column=0, sticky="w", padx=8, pady=4)
        default_rows = str((saved or {}).get("grid_rows", "2"))
        grid_rows_var = Tk.StringVar(value=default_rows)
        Tk.Entry(grid_frame, textvariable=grid_rows_var, width=6,
                 bg=COLORS["entry"], fg=COLORS["text"],
                 font=("Segoe UI", 10)).grid(row=0, column=1, padx=4)

        Tk.Label(grid_frame, text="Grid Columns (m):", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=0, column=2, sticky="w", padx=8)
        default_cols = str((saved or {}).get("grid_cols", "2"))
        grid_cols_var = Tk.StringVar(value=default_cols)
        Tk.Entry(grid_frame, textvariable=grid_cols_var, width=6,
                 bg=COLORS["entry"], fg=COLORS["text"],
                 font=("Segoe UI", 10)).grid(row=0, column=3, padx=4)
        
        Tk.Label(grid_frame, text="(Canvas refreshes after n × m plots)", 
                 bg=COLORS["panel"], fg="#95A5A6", font=("Segoe UI", 9)).grid(
                     row=0, column=4, sticky="w", padx=12)

        # ------------------------------------------------------------------
        # SUBSECTION: SEQUENCE OPTIONS
        # PURPOSE: Batch Run ID, fit-overlay, and other run-sequencing
        #          options for the batch.
        # ------------------------------------------------------------------
        # ── Sequence options ───────────────────────────────────────────────────
        options_frame = Tk.LabelFrame(win, text="  Sequence Options  ",
                                       bg=COLORS["panel"], fg=COLORS["accent"],
                                       font=("Segoe UI", 10, "bold"))
        options_frame.pack(fill=Tk.X, padx=12, pady=(8, 2))

        Tk.Label(options_frame,
                 text="Inter-board pause (s):", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=0, column=0, sticky="w", padx=8, pady=4)
        default_inter = str((saved or {}).get("inter_delay", "5"))
        inter_delay_var = Tk.StringVar(value=default_inter)
        Tk.Entry(options_frame, textvariable=inter_delay_var, width=6,
                 bg=COLORS["entry"], fg=COLORS["text"],
                 font=("Segoe UI", 10)).grid(row=0, column=1, padx=4)

        Tk.Label(options_frame,
                 text="Inter-channel pause (s):", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=0, column=2, sticky="w", padx=8, pady=4)
        default_inter_ch = str((saved or {}).get("inter_channel_delay", "0"))
        inter_channel_delay_var = Tk.StringVar(value=default_inter_ch)
        Tk.Entry(options_frame, textvariable=inter_channel_delay_var, width=6,
                 bg=COLORS["entry"], fg=COLORS["text"],
                 font=("Segoe UI", 10)).grid(row=0, column=3, padx=4)

        Tk.Label(options_frame,
                 text="Auto-save results:", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=0, column=4, sticky="w", padx=8)
        default_autosave = (saved or {}).get("autosave", True)
        autosave_var = Tk.BooleanVar(value=default_autosave)
        Tk.Checkbutton(options_frame, variable=autosave_var,
                        bg=COLORS["panel"], activebackground=COLORS["panel"],
                        selectcolor=COLORS["entry"]).grid(
                            row=0, column=5, padx=4)

        # --- NEW: Batch Run ID Input ---
        Tk.Label(options_frame, text="Batch Run ID:", bg=COLORS["panel"],
                 fg=COLORS["text"], font=("Segoe UI", 10)).grid(
                     row=1, column=0, sticky="w", padx=8, pady=(4, 4))
        
        default_batch_id = (saved or {}).get("batch_id", "")
        batch_id_var = Tk.StringVar(value=default_batch_id)
        Tk.Entry(options_frame, textvariable=batch_id_var, width=15,
                 bg=COLORS["entry"], fg=COLORS["text"],
                 font=("Segoe UI", 10)).grid(row=1, column=1, padx=4, pady=(4, 4))
        # -------------------------------

        # Overlay fit on I-V curve in Measurement Tab
        default_overlay_fit = (saved or {}).get("overlay_fit_on_iv", True)
        overlay_fit_var = Tk.BooleanVar(value=default_overlay_fit)
        # Keep fit_tab_var for backward-compat with _persist_state / _batch_config
        fit_tab_var = Tk.StringVar(value="Measurement Tab")
        Tk.Checkbutton(options_frame,
                        text="Overlay fit on I-V curve",
                        variable=overlay_fit_var,
                        bg=COLORS["panel"], fg=COLORS["text"],
                        activebackground=COLORS["panel"],
                        selectcolor=COLORS["entry"],
                        font=("Segoe UI", 10)).grid(
                            row=1, column=2, columnspan=4, sticky="w",
                            padx=8, pady=(4, 4))

        # NOTE: the status label and btn_bar frame are created earlier
        # (pinned to the window bottom right after the title bar) so the
        # action buttons always stay visible regardless of window height.

        # ------------------------------------------------------------------
        # SUBSECTION: STATE PERSISTENCE / CONFIRM / CANCEL / RELOAD
        # PURPOSE: Helper functions bound to the footer buttons -- save
        #          current form state, validate & confirm the batch setup,
        #          cancel without applying, or reload from batch.mac.
        # ------------------------------------------------------------------
        def _persist_state():
            self._batch_saved_state = {
                "ard_port": batch_ard_port_var.get(),
                "baud":     batch_baud_var.get(),
                "boards": [
                    {
                        "label":       r["label"].get(),
                        "slave_id":    r["slave_id"].get(),
                        "pin":         r["pin"].get(),
                        "enabled":     r["enabled"].get(),
                        "participate": r["participate"].get(),
                    }
                    for r in self._batch_board_rows
                ],
                "sipms": [
                    {
                        "label":       r["label"].get(),
                        "board_label": r["board_var"].get(),
                        "pin":         r["pin"].get(),
                        "rbias":       r["rbias"].get(),
                        "ncells":      r["ncells"].get(),
                        "enabled":     r["enabled"].get(),
                        "participate": r["participate"].get(),
                    }
                    for r in self._batch_sipm_rows
                ],
                "inter_delay":          inter_delay_var.get(),
                "inter_channel_delay":  inter_channel_delay_var.get(),
                "autosave":             autosave_var.get(),
                "overlay_fit_on_iv":    overlay_fit_var.get(),
                "fit_tab":              fit_tab_var.get(),
                "grid_rows":            grid_rows_var.get(),
                "grid_cols":            grid_cols_var.get(),
                "batch_id":             batch_id_var.get(),
            }

        def _confirm_batch_setup():
            boards = [r for r in self._batch_board_rows if r["participate"].get()]
            if not boards:
                msg.showwarning("Batch Mode",
                                 "No board is checked 'Use in I-V'!", parent=win)
                return
            sipms = [r for r in self._batch_sipm_rows if r["participate"].get()]
            if not sipms:
                msg.showwarning("Batch Mode",
                                 "No SiPM slot is checked 'Use in I-V'!", parent=win)
                return
            enabled_board_labels = {r["label"].get() for r in boards}
            if not any(r["board_var"].get() in enabled_board_labels for r in sipms):
                msg.showwarning(
                    "Batch Mode",
                    "None of the enabled SiPM slots are assigned to an "
                    "enabled board.", parent=win)
                return
            try:
                inter_d = float(inter_delay_var.get())
            except ValueError:
                msg.showerror("Batch Mode",
                               "Invalid numeric value in options.", parent=win)
                return
            try:
                inter_ch_d = float(inter_channel_delay_var.get())
            except ValueError:
                inter_ch_d = 0.0

            _persist_state()

            self._batch_config = {
                "boards":                boards,
                "sipms":                 sipms,
                "ard_port":              batch_ard_port_var.get(),
                "baud":                  int(batch_baud_var.get()),
                "inter_delay":           inter_d,
                "inter_channel_delay":   inter_ch_d,
                "autosave":              autosave_var.get(),
                "overlay_fit_on_iv":     overlay_fit_var.get(),
                "fit_tab":               fit_tab_var.get(),
                "grid_rows":             int(grid_rows_var.get() if grid_rows_var.get().isdigit() else 2),
                "grid_cols":             int(grid_cols_var.get() if grid_cols_var.get().isdigit() else 2),
                "batch_id":              batch_id_var.get().strip(),
                "status_var":            self._batch_status_var,
                "win":                   win,
            }
            # Batch is now configured & validated -- select Batch Run mode
            # on the main GUI so the user's next click of START TEST
            # launches it, then close this window.
            self.run_mode_var.set('batch')
            win.grab_release()
            _status_close_serial()
            win.destroy()
            msg.showinfo("Batch Mode",
                          "Batch setup confirmed.\n"
                          "Go to the main window and click START TEST "
                          "to begin the batch run.")

        def _cancel():
            _persist_state()
            _status_close_serial()
            win.destroy()

        Tk.Button(btn_bar, text="✔  Confirm Batch Setup",
                  bg=COLORS["success"], fg="white",
                  font=("Segoe UI", 11, "bold"), relief="flat", padx=10,
                  command=_confirm_batch_setup).pack(side=Tk.LEFT, padx=(0, 8))

        Tk.Button(btn_bar, text="💾  Save Settings",
                  bg=COLORS["accent2"], fg="white",
                  font=("Segoe UI", 10), relief="flat", padx=8,
                  command=_persist_state).pack(side=Tk.LEFT, padx=(0, 8))

        Tk.Button(btn_bar, text="✖  Cancel",
                  bg=COLORS["danger"], fg="white",
                  font=("Segoe UI", 10), relief="flat", padx=8,
                  command=_cancel).pack(side=Tk.LEFT)

        def reload_batch_mac():
            self.load_batch_macro("batch.mac")
            _cancel() # Close and immediately re-open window to apply loaded data
            self.open_batch_config()
            
        Tk.Button(btn_bar, text="Load batch",
                  bg="#34495E", fg="white",
                  font=("Segoe UI", 10, "bold"), relief="flat", padx=8,
                  command=reload_batch_mac).pack(side=Tk.RIGHT, padx=(0,8))

        win.protocol("WM_DELETE_WINDOW", _cancel)

        Tk.Label(btn_bar,
                 text="I-V params (Start/End/Step/Delay/Threshold) are taken from main window.",
                 bg=COLORS["bg"], fg="#7F8C8D",
                 font=("Segoe UI", 8)).pack(side=Tk.RIGHT, padx=4)

    def run_batch_sequence(self):
        cfg      = self._batch_config
        boards   = cfg["boards"]
        sipms    = cfg["sipms"]
        inter_d  = cfg["inter_delay"]
        inter_ch_d = cfg.get("inter_channel_delay", 0.0)
        overlay_fit = cfg.get("overlay_fit_on_iv", True)
        autosave = cfg["autosave"]
        status_v = cfg["status_var"]

        batch_ser = None
        try:
            batch_ser = serial.Serial(cfg["ard_port"], cfg["baud"], timeout=2)
            time.sleep(2)                          
            while batch_ser.in_waiting:            
                batch_ser.readline()
            status_v.set(f"Arduino connected on {cfg['ard_port']}")
            self.window.update_idletasks()
        except Exception as e:
            msg.showerror("Batch Mode", f"Cannot open Arduino port:\n{e}")
            return

        if self.search_flag == 0 or self.instrument is None:
            msg.showerror("Batch Mode",
                          "Keithley not connected!\n"
                          "Please connect the power supply first.")
            if batch_ser:
                batch_ser.close()
            return

        if self.RUN_IV_HV() == 0:
            if batch_ser:
                batch_ser.close()
            return

        try:
            _, _sv = self.is_number(self.start_voltage.get())
            _, _ev = self.is_number(self.end_voltage.get())
            _, _st = self.is_number(self.step_voltage.get())
            _, _dt = self.is_number(self.delay_time.get())
            n_steps = max(1, int(abs(_ev - _sv) / max(_st, 0.001)))
            sweep_timeout = max(120, n_steps * (_dt + 2) * 10)
        except Exception:
            sweep_timeout = 600        

        total_boards     = len(boards)
        total_sipms_done = 0
        self.batch_mode_active = True
        self.batch_pause_flag = 0
        self.batch_stop_flag = 0
        self.batch_pause_btn.config(text='BATCH PAUSE', bg='#E0E0E0')
        self._batch_status_var.set("Batch running…")

        # --- NEW: Create Mother Folder ---
        batch_id = cfg.get("batch_id", "")
        if not batch_id:
            import random
            batch_id = "Auto" + str(random.randint(1000, 9999))
            
        ts_mother = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
        self.current_mother_folder = f"./Results/Batch_Run_{batch_id}_{ts_mother}"
        os.makedirs(self.current_mother_folder, exist_ok=True)
        # ---------------------------------

        grid_rows = max(1, cfg.get("grid_rows", 2))
        grid_cols = max(1, cfg.get("grid_cols", 2))
        plots_per_page = grid_rows * grid_cols
        current_plot_count = 0

        try:
            for board_idx, board in enumerate(boards):
                self._batch_checkpoint()
                board_label = board["label"].get()

                try:
                    slave_id = int(board["slave_id"].get())
                    board_pin = int(board["pin"].get())
                except (ValueError, Tk.TclError):
                    status_v.set(f"[Board {board_idx+1}/{total_boards}] "
                                 f"Skipped {board_label}: invalid Slave ID / Pin")
                    self.window.update_idletasks()
                    continue

                board_sipms = [s for s in sipms
                               if s["board_var"].get() == board_label]
                if not board_sipms:
                    status_v.set(f"[Board {board_idx+1}/{total_boards}] "
                                 f"{board_label}: no SiPM slots assigned — skipping")
                    self.window.update_idletasks()
                    continue

                n_sipms_board = len(board_sipms)
                cmd_board_on  = f"Board {slave_id} Update {board_pin} 1"
                cmd_board_off = f"Board {slave_id} Update {board_pin} 0"

                status_v.set(f"[Board {board_idx+1}/{total_boards}] "
                             f"Switching ON → {board_label}  ({cmd_board_on})")
                self.window.update_idletasks()
                reply = self._batch_send_arduino(batch_ser, cmd_board_on)
                print(f"[Batch] Arduino reply to '{cmd_board_on}': {reply}")
                time.sleep(0.5)

                for sipm_idx, sipm in enumerate(board_sipms):
                    self._batch_checkpoint()
                    sipm_label = sipm["label"].get()
                    full_label = f"{board_label}_{sipm_label}"
                    total_sipms_done += 1

                    try:
                        sipm_pin = int(sipm["pin"].get())
                        cmd_sipm_on  = f"Board {slave_id} Update {sipm_pin} 1"
                        cmd_sipm_off = f"Board {slave_id} Update {sipm_pin} 0"
                        reply = self._batch_send_arduino(batch_ser, cmd_sipm_on)
                        print(f"[Batch] SiPM relay ON  '{cmd_sipm_on}': {reply}")
                        time.sleep(0.3)
                    except Exception as e:
                        print(f"[Batch] SiPM relay switch error: {e}")
                        cmd_sipm_off = None

                    self.module_name.set(full_label)
                    self.xp = [];  self.yp = [];  self.ypp = []
                    self.xp_ap = [];  self.temp_arr = []
                    self.humid_arr = [];  self.time_arr = []
                    self.curr_std_arr = []
                    self.warn_flag   = 0
                    self.stop_flag   = 0
                    self.pause_plot  = 0
                    self.run_flag    = 1      
                    self.sim_flag    = 0
                    self.run_index   = 0
                    self.run_init_flg = 0
                    self.legn_flag   = 0

                    flag1, current_th_num = self.is_number(self.current_th.get())
                    flag2, start_v  = self.is_number(self.start_voltage.get())
                    flag3, end_v    = self.is_number(self.end_voltage.get())
                    flag4, step_v   = self.is_number(self.step_voltage.get())
                    flag5, delay_t  = self.is_number(self.delay_time.get())
                    # Ramp-down step now comes straight from the main
                    # window's "Ramp Down (V)" field -- read fresh here so
                    # it's never stale, and never overwritten afterward
                    # (see batch_mode_active guard in start_process).
                    flag6, down_s   = self.is_number(self.down_step_voltage.get())

                    self.start_vol     = start_v
                    self.end_vol       = end_v
                    self.step_vol      = step_v
                    self.down_step_vol = down_s
                    self.time_delay    = delay_t
                    self.polarinit     = self.chk_polarity(end_v, start_v)
                    
                    # ── Setup / Refresh Canvas Grid ──────────────────────────
                    if current_plot_count % plots_per_page == 0:
                        # Measurement Grid Setup (unchanged — grid stays for Measurement tab)
                        self.figure.clf()
                        self.figure.subplots_adjust(left=0.10, right=0.88, top=0.92, bottom=0.15, hspace=0.7, wspace=0.4)
                        self.batch_axes = []
                        self.batch_axes2 = []
                        for i in range(plots_per_page):
                            ax = self.figure.add_subplot(grid_rows, grid_cols, i + 1)
                            ax2 = ax.twinx()
                            self.batch_axes.append(ax)
                            self.batch_axes2.append(ax2)
                        self.figure_canvas.draw()

                        # Ensure Measurement tab is visible at start of new grid
                        self.plot_notebook.select(self.tab_measure)
                        self.window.update_idletasks()

                    # Assign current measurement axis target based on progress
                    ax_index = current_plot_count % plots_per_page
                    self.ax = self.batch_axes[ax_index]
                    self.ax2 = self.batch_axes2[ax_index]
                    
                    self.ax.clear()
                    self.ax2.clear()
                    self.errbar_container = None  # was attached to a now-cleared/different axes

                    # Link plot lines strictly to the targeted subplot
                    self.plot1, = self.ax.plot([], [], 'o-', color='#3498DB', markersize=3, label="Measured I-V Data")
                    self.plot2, = self.ax.plot([], [], 'x', color='#E74C3C', markersize=3, label="Validation of applied vs measured voltage")
                    self.plot3, = self.ax.plot([], [], 'b', linestyle='None', label=None)
                    self.plot4, = self.ax2.plot([], [], 'ro', linestyle='None', label=None)
                    self.plot5, = self.ax2.plot([], [], 'bd', markersize=2, label="Temp")
                    self.plot6, = self.ax2.plot([], [], 'ms', markersize=2, label="Humidity")

                    # Subplot UI Formatting
                    self.ax.set_title(full_label, fontsize=9, pad=3, fontweight='bold')
                    self.ax.set_xlabel('Voltage (V)', color='green', fontsize=8)
                    # Match single-mode's live I-V y-label style (see
                    # plot_VI_graph): red, bold, "Current / in nA". Uses a
                    # normal set_ylabel (two lines via \n) rather than
                    # multicolor_ylabel's fixed-offset AnchoredOffsetbox,
                    # which on these small grid tiles was overshooting the
                    # subplot's left edge and getting clipped off entirely.
                    self.ax.set_ylabel("Current\nin nA", fontsize=8, fontweight='bold', color='red', labelpad=4)
                    self.ax.set_yscale(self.scale_var.get() or 'linear')
                    if not self.auto_yscale_var.get():
                        flag_min, y_min = self.is_number(self.ymin_var.get())
                        flag_max, y_max = self.is_number(self.ymax_var.get())
                        if flag_min and flag_max and y_max > y_min:
                            if self.scale_var.get() == 'log' and y_min <= 0:
                                y_min = 0.1
                            self.ax.set_ylim(y_min, y_max)
                    self.ax.tick_params(axis='both', labelsize=8)
                    self.ax2.tick_params(axis='both', labelsize=8)
                    
                    self.ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5, color='gray')
                    self.ax.set_facecolor('white')

                    if self.var.get() == 0:
                        self.ax2.set_visible(False)
                    else:
                        self.ax2.set_visible(True)
                        self.ax2.yaxis.set_label_position("right")
                        self.ax2.yaxis.tick_right()
                        #self.ax2.set_ylabel('Temp (\u00B0C) / Humidity (%)', color='m', fontsize=7, fontweight='bold', rotation=270, labelpad=15)
                        self.ax2.text(1.08, 0.52, 'Temp (\u00B0C)', color='b', 
              fontsize=7, fontweight='bold', rotation=270, 
              transform=self.ax2.transAxes, ha='left', va='bottom')

# 3. Place the second part of the label (Humidity)
# y=0.48 starts it slightly below the center.
                        self.ax2.text(1.08, 0.48, ' / Humidity (%)', color='m', 
              fontsize=7, fontweight='bold', rotation=270, 
              transform=self.ax2.transAxes, ha='left', va='top')


                    h1, l1 = self.ax.get_legend_handles_labels()
                    h2, l2 = self.ax2.get_legend_handles_labels()
                    if h1 or h2:
                        self.ax.legend(h1+h2, l1+l2, bbox_to_anchor=(0.5, -0.30), loc='upper center', ncol=3, fontsize=6, framealpha=0.9)

                    self.curr_th = current_th_num * 1e-6
                    self.set_current_threshold(self.curr_th)

                    self.show_green_light()

                    status_v.set(
                        f"[Board {board_idx+1}/{total_boards}] {board_label} → "
                        f"SiPM {sipm_idx+1}/{n_sipms_board}: {sipm_label}  "
                        f"Running I-V ({start_v} V → {end_v} V, step {step_v} V)")
                    self.window.update_idletasks()

                    self.start_process()

                    t_start = time.time()
                    paused_elapsed = 0.0
                    while (self.run_flag == 1 or self.awaiting_reconnect) and self.stop_flag == 0:
                        if self.batch_pause_flag or self.awaiting_reconnect:
                            # Freeze the timeout clock while paused -- the
                            # sweep itself is already frozen (pause_plot=1
                            # stops auto_run_process from scheduling its
                            # next step), so don't let a long pause look
                            # like a sweep timeout. A disconnect
                            # (awaiting_reconnect) also has to hold here,
                            # not just batch_pause_flag -- otherwise
                            # run_flag==0 (set by handle_disconnect_
                            # during_sweep) would immediately fall through
                            # to ramp-down/next-channel instead of waiting
                            # for the user to reconnect and press
                            # RESUME/BATCH RESUME.
                            p0 = time.time()
                            while (self.batch_pause_flag or self.awaiting_reconnect) and self.stop_flag == 0:
                                self.window.update()
                                time.sleep(0.1)
                            paused_elapsed += time.time() - p0
                            continue
                        self.window.update()
                        time.sleep(0.05)
                        if time.time() - t_start - paused_elapsed > sweep_timeout:
                            print(f"[Batch] Sweep timeout for {full_label}; "
                                  f"forcing stop.")
                            self.run_flag  = 0
                            self.stop_flag = 1
                            break

                    status_v.set(
                        f"[Board {board_idx+1}/{total_boards}] {board_label} → "
                        f"SiPM {sipm_idx+1}/{n_sipms_board}: {sipm_label}  "
                        f"Ramping down…")
                    self.window.update_idletasks()
                    
                    try:
                        self.ramp_down_complete = False
                        print('Ramp Down Step from batch:: ', down_s)
                        # NOTE: ramp_down_zero must run unconditionally here.
                        # Previously this was gated on `self.warn_flag == 0`,
                        # but warn_flag gets set to 1 whenever a sweep ends
                        # via the current-limit warning path (very common on
                        # forward I-V, which routinely hits the current
                        # threshold). With the gate, ramp_down_zero was never
                        # called in that case, ramp_down_complete never
                        # turned True, and the wait loop below would hang for
                        # up to an hour -- batch_mode_active stayed True the
                        # whole time, which is why a new batch run reported
                        # "already in progress" and the next channel never
                        # started.
                        self.ramp_down_zero(down_s, delay_t)
                        t_ramp_start = time.time()
                        while not self.ramp_down_complete:
                            self.window.update()
                            time.sleep(0.05)
                            if self.batch_stop_flag:
                                # User pressed BATCH STOP while ramping down
                                # -- don't sit here for up to an hour, bail
                                # out as soon as the flag is seen.
                                print(f"[Batch] Stop requested during ramp-down for {full_label}; breaking wait.")
                                break
                            if time.time() - t_ramp_start > 3600:   
                                   print(f"[Batch] Ramp-down timeout for {full_label}; forcing continue.")
                                   break
                    except Exception as e:
                            print(f"[Batch] Ramp-down error: {e}")    

                    # Channel is now safely ramped down -- this is a safe
                    # point to honor a batch-stop request and unwind out
                    # of the whole sequence (never mid-sweep, mid-ramp).
                    if self.batch_stop_flag:
                        # Mirror the normal end-of-channel/end-of-board
                        # relay cleanup so we don't leave the SiPM or board
                        # relay energized when bailing out early.
                        try:
                            if cmd_sipm_off:
                                reply = self._batch_send_arduino(batch_ser, cmd_sipm_off)
                                print(f"[Batch] SiPM relay OFF '{cmd_sipm_off}' (stop): {reply}")
                            reply = self._batch_send_arduino(batch_ser, cmd_board_off)
                            print(f"[Batch] Board relay OFF '{cmd_board_off}' (stop): {reply}")
                        except Exception as e:
                            print(f"[Batch] Relay cleanup on stop error: {e}")
                        raise BatchStopRequested()

                    # Batch Analysis Trigger (Using Targeted Axes)
                    fit_success = False
                    fit_popt = None
                    # Clear any Rq fit left over from a previous SiPM so a
                    # reverse-mode (or skipped) run never picks up a stale
                    # forward-mode result in the autosave step below.
                    self.last_rq_result = None
                    mode = self.analysis_mode_var.get()
                    if mode == "reverse" and self.calc_vbd_var.get() and len(self.xp) >= 5:
                        try:
                            volts = np.array(self.xp)
                            currents_nA = np.array(self.yp) * 1e9 if self.sim_flag == 1 else np.array(self.yp)
                            current_std = None
                            if hasattr(self, 'curr_std_arr') and self.curr_std_arr and len(self.curr_std_arr) == len(volts):
                                current_std = np.array(self.curr_std_arr, dtype=float)
                            v_bd_deriv = find_vbd_derivative(volts, currents_nA)
                            popt, success, perr = optimize_fit(volts, currents_nA, v_bd_deriv, user_params=self.user_fit_params, current_std=current_std)
                            fit_popt, fit_success = popt, success

                            if overlay_fit and success:
                                v_bd_fit = popt[0]
                                v_smooth = np.linspace(min(volts), min(max(volts), popt[1]-0.1), 500)
                                i_fit_nA = dinu_eq8_model(v_smooth, *popt)

                                # Overlay the fit on the Measurement Tab I-V curve
                                self.ax.plot(v_smooth, i_fit_nA, 'g--', linewidth=1.5, label="Fit Model")
                                self.ax.axvline(v_bd_fit, color='red', linestyle='--', alpha=0.7)
                                y_val_nA = dinu_eq8_model(v_bd_fit, *popt)
                                self.ax.plot(v_bd_fit, y_val_nA, 'rx', markersize=6)

                                # Add summary box
                                summary = f"Vbd:{popt[0]:.1f}V"
                                if self.show_dcr_var.get() and abs(self.C_ucell) > 0:
                                    DCR = popt[3]*1e-9 / (self.C_ucell*1e3)
                                    summary += f" DCR:{DCR:.1f}kHz"

                                self.ax.text(0.02, 0.98, summary, transform=self.ax.transAxes,
                                             verticalalignment='top', fontsize=6,
                                             bbox=dict(boxstyle="round", fc="white", alpha=0.8),
                                             color="black")
                                # Keep the same legend style/position used for
                                # this subplot (bbox_to_anchor centered strip
                                # below the axes) -- just re-pull the handles
                                # so the newly added "Fit Model" line is
                                # folded in, rather than replacing it with a
                                # differently-styled lower-right legend.
                                h1, l1 = self.ax.get_legend_handles_labels()
                                h2, l2 = self.ax2.get_legend_handles_labels()
                                if h1 or h2:
                                    self.ax.legend(h1 + h2, l1 + l2, bbox_to_anchor=(0.5, -0.30),
                                                   loc='upper center', ncol=3, fontsize=6, framealpha=0.9)
                                self.figure_canvas.draw()

                            # Render a single full-size analysis plot for THIS SiPM into the
                            # Analysis tab (replaces previous SiPM's plot; no grid here).
                            
                            # --- FIX: Briefly select the Analysis Tab so text scales properly ---
                            self.plot_notebook.select(self.tab_analysis)
                            self.window.update_idletasks()
                            
                            self.run_breakdown_analysis(title=full_label)
                            self.canvas_analysis.draw()
                        except Exception as e:
                            print(f"[Batch] Breakdown analysis error: {e}")

                    elif mode == "forward" and self.show_rq_var.get() and len(self.xp) >= 6:
                        # Quenching-resistance analysis for forward-bias batch
                        # runs. This branch was previously missing entirely --
                        # only the reverse/breakdown-voltage branch above ran
                        # in batch mode, so ticking "Quench R" in Forward mode
                        # never produced an Rq fit during a batch sequence.
                        try:
                            self.plot_notebook.select(self.tab_analysis)
                            self.window.update_idletasks()

                            # Per-channel R_bias / N (Microcells) override:
                            # blank boxes in the batch table fall back to
                            # the main window's Forward Bias Rq settings.
                            # Swap the shared StringVars in for the
                            # duration of this channel's fit, then restore
                            # them so the main window is left untouched.
                            _orig_rbias  = self.rq_rbias_var.get()
                            _orig_ncells = self.rq_ncells_var.get()
                            _ch_rbias, _ch_ncells = self._resolve_rq_params(sipm)
                            self.rq_rbias_var.set(_ch_rbias)
                            self.rq_ncells_var.set(_ch_ncells)
                            try:
                                self.run_quench_resistance_analysis(title=full_label)
                                fit_success = getattr(self, 'last_rq_result', None) is not None
                            finally:
                                self.rq_rbias_var.set(_orig_rbias)
                                self.rq_ncells_var.set(_orig_ncells)
                            self.canvas_analysis.draw()
                        except Exception as e:
                            print(f"[Batch] Quench resistance (Rq) analysis error: {e}")

                    # ── Save results for this SiPM (always, regardless of fit tab) ──
                    if len(self.xp) > 0:
                        if autosave:
                            self._batch_autosave(full_label, ax_index=ax_index,
                                                  fit_popt=fit_popt, fit_success=fit_success,
                                                  analysis_mode=mode,
                                                  rq_result=getattr(self, 'last_rq_result', None))
                        else:
                            try:
                                _std_b = list(self.curr_std_arr) if hasattr(self, 'curr_std_arr') and len(self.curr_std_arr) == len(self.xp) else [0.0] * len(self.xp)
                                pd.DataFrame({
                                    "VOLTS":         self.xp,
                                    "CURRNT_NAMP":   self.yp,
                                    "CURR_STD_NAMP": _std_b,
                                    "TEMP_DEGC":     self.temp_arr,
                                    "RH_PRCNT":      self.humid_arr,
                                    "TIME":          self.time_arr,
                                }).to_csv("temp.csv", index=False)
                            except Exception:
                                pass
                                
                    # --- FIX: Return to the Measurement Tab after saving is complete ---
                    self.plot_notebook.select(self.tab_measure)
                    self.window.update_idletasks()
                    # -------------------------------------------------------------------


                    # ── Post-Processing Tab Switch and Wait ─────────────────
                    if cmd_sipm_off:
                        try:
                            reply = self._batch_send_arduino(batch_ser, cmd_sipm_off)
                            print(f"[Batch] SiPM relay OFF '{cmd_sipm_off}': {reply}")
                            time.sleep(0.3)
                        except Exception as e:
                            print(f"[Batch] SiPM relay OFF error: {e}")

                    current_plot_count += 1  # Add this line back right here

                    # ── Inter-channel pause (between SiPMs on same board) ──────
                    remaining_sipms = len(board_sipms) - (sipm_idx + 1)
                    if remaining_sipms > 0 and inter_ch_d > 0:
                        status_v.set(f"[Board {board_idx+1}/{total_boards}] "
                                     f"Inter-channel pause {inter_ch_d} s…")
                        self.window.update_idletasks()
                        deadline_ch = time.time() + inter_ch_d
                        while time.time() < deadline_ch:
                            self._batch_checkpoint()
                            self.window.update()
                            time.sleep(0.1)
                status_v.set(f"[Board {board_idx+1}/{total_boards}] "
                             f"Switching OFF → {board_label}  ({cmd_board_off})")
                self.window.update_idletasks()
                reply = self._batch_send_arduino(batch_ser, cmd_board_off)
                print(f"[Batch] Arduino reply to '{cmd_board_off}': {reply}")
                time.sleep(0.3)

                if board_idx < total_boards - 1 and inter_d > 0:
                    status_v.set(f"[Board {board_idx+1}/{total_boards}] "
                                 f"Waiting {inter_d} s before next board…")
                    self.window.update_idletasks()
                    deadline = time.time() + inter_d
                    while time.time() < deadline:
                        self._batch_checkpoint()
                        self.window.update()
                        time.sleep(0.1)

        except BatchStopRequested:
            status_v.set("Batch stopped by user.")
            self.window.update_idletasks()
            self._batch_status_var.set("Batch stopped by user")
        finally:
            self.batch_mode_active = False
            if batch_ser and batch_ser.is_open:
                batch_ser.close()                    
            self.batch_pause_btn.config(text='BATCH PAUSE', bg='#E0E0E0')

        if self.batch_stop_flag:
            self.batch_stop_flag = 0
            self.show_yellow_light()   # stopped mid-run → amber
            try:
                cfg["win"].deiconify()
                cfg["win"].lift()
            except Exception:
                pass
            msg.showinfo("Batch Mode",
                         f"Batch stopped by user.\n{total_sipms_done} SiPM "
                         f"measurement(s) completed before stopping.")
            return

        self.show_complete_light()   # all boards done → teal complete
        self._batch_status_var.set(
            f"Batch complete — {total_sipms_done} SiPM(s)")
        status_v.set(f"  Batch complete  –  {total_sipms_done} SiPM "
                     f"measurement(s) across {total_boards} board(s).")

        try:
            cfg["win"].deiconify()
            cfg["win"].lift()
        except Exception:
            pass

        msg.showinfo("Batch Mode",
                     f"Batch run complete!\n{total_sipms_done} SiPM "
                     f"measurement(s) across {total_boards} board(s).")
        # User pressed OK on the "Batch complete" pop-up → drop the status
        # from teal "Complete" down to amber "Paused / Hold".
        self.show_yellow_light()

    def _batch_autosave(self, label, ax_index=None, fit_popt=None, fit_success=False,
                        analysis_mode=None, rq_result=None):
        try:
            safe_label = re.sub(r'[^A-Za-z0-9_\-]', '_', label)
            ts = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
            
            # --- NEW: Route to Mother Folder if it exists ---
            if hasattr(self, 'current_mother_folder') and self.current_mother_folder:
                directory = f"{self.current_mother_folder}/SiPM_{safe_label}_{ts}"
            else:
                directory = f"./Results/Batch_{safe_label}_{ts}"
            # ------------------------------------------------
            
            os.makedirs(directory, exist_ok=True)
            base = f"{directory}/{ts}_{safe_label}"

            # 1. Per-SiPM data log (always, one file per SiPM)
            std_as = list(self.curr_std_arr) if hasattr(self, 'curr_std_arr') and len(self.curr_std_arr) == len(self.xp) else [0.0] * len(self.xp)
            alldata = pd.DataFrame({
                "VOLTS":           self.xp,
                "CURRNT_NAMP":     self.yp,
                "CURR_STD_NAMP":   std_as,
                "TEMP_DEGC":       self.temp_arr,
                "RH_PRCNT":        self.humid_arr,
                "TIME":            self.time_arr,
            })
            alldata.to_csv(base + "_Result_Log.csv", index=False)

            # 2. Measurement composite grid (overview of all channels on this page)
            self.figure.savefig(base + "_Composite_IV_Grid.png")
            try:
                # fig_analysis now holds a single full-size rich plot for THIS SiPM only
                self.fig_analysis.savefig(base + "_Analysis_SinglePlot.png")
            except Exception:
                pass

            # 3. Separate single I-V plot for THIS SiPM only, redrawn standalone
            #    (independent of the composite grid layout / other subplots)
            try:
                single_fig = plt.Figure(figsize=(6, 4.5), dpi=150)
                single_ax = single_fig.add_subplot(111)
                volts = np.array(self.xp)
                currents = np.array(self.yp)

                single_ax.plot(volts, currents, 'o-', color='#3498DB',
                                markersize=4, label="Measured Data")

                if analysis_mode == "forward" and rq_result is not None:
                    # Forward-bias quenching-resistance run: draw the same
                    # linear Rq fit(s) shown on the Analysis tab, instead of
                    # the reverse-bias breakdown model below. Using the
                    # wrong model here was the original bug -- this single
                    # PNG used to always assume dinu_eq8_model (Vbd fit)
                    # regardless of which analysis actually produced the
                    # numbers shown in the GUI.
                    try:
                        r = rq_result
                        if r.get('show_r1'):
                            V1, m1, c1 = r['V1'], r['m1'], r['c1']
                            V1_line = np.linspace(max(0.0, V1.min() * 0.9), V1.max() * 1.05, 200)
                            single_ax.plot(V1_line, m1 * V1_line + c1, color='steelblue',
                                            linewidth=1.5, linestyle='--',
                                            label=f"Fit R1: Rq={r['rq1_disp_str']}")
                        if r.get('show_r2'):
                            V2, m2, c2 = r['V2'], r['m2'], r['c2']
                            V2_line = np.linspace(max(0.0, V2.min() * 0.9), V2.max() * 1.05, 200)
                            single_ax.plot(V2_line, m2 * V2_line + c2, color='darkorange',
                                            linewidth=1.5, linestyle='--',
                                            label=f"Fit R2: Rq={r['rq2_disp_str']}")
                        if r.get('show_r1') and r.get('show_r2'):
                            single_ax.axvline(r['V_knee'], color='red', linestyle=':', alpha=0.7)

                        summary_lines = []
                        if r.get('show_r1'):
                            summary_lines.append(f"Rq: {r['rq1_disp_str']}")
                        if r.get('show_r2'):
                            summary_lines.append(f"Rq: {r['rq2_disp_str']}")
                        summary = "\n".join(summary_lines)
                        if summary:
                            single_ax.text(0.02, 0.98, summary, transform=single_ax.transAxes,
                                            verticalalignment='top', fontsize=9,
                                            bbox=dict(boxstyle="round", fc="white", alpha=0.85),
                                            color="black")
                    except Exception as e:
                        print(f"[Batch] Single-plot Rq overlay error: {e}")

                elif fit_success and fit_popt is not None:
                    try:
                        currents_nA = currents * 1e9 if self.sim_flag == 1 else currents
                        v_bd_fit = fit_popt[0]
                        v_smooth = np.linspace(min(volts), min(max(volts), fit_popt[1] - 0.1), 500)
                        i_fit_nA = dinu_eq8_model(v_smooth, *fit_popt)
                        y_val_nA = dinu_eq8_model(v_bd_fit, *fit_popt)

                        single_ax.plot(v_smooth, i_fit_nA, 'g--', linewidth=1.5, label="Fit Model")
                        single_ax.axvline(v_bd_fit, color='red', linestyle='--', alpha=0.7)
                        single_ax.plot(v_bd_fit, y_val_nA, 'rx', markersize=8, label="Breakdown Point")

                        summary = f"Vbd: {v_bd_fit:.2f} V"
                        if self.show_dcr_var.get() and abs(self.C_ucell) > 0:
                            DCR = fit_popt[3] * 1e-9 / (self.C_ucell * 1e3)
                            summary += f"\nDCR: {DCR:.2f} kHz"
                        single_ax.text(0.02, 0.98, summary, transform=single_ax.transAxes,
                                        verticalalignment='top', fontsize=9,
                                        bbox=dict(boxstyle="round", fc="white", alpha=0.85),
                                        color="black")
                    except Exception as e:
                        print(f"[Batch] Single-plot fit overlay error: {e}")

                single_ax.set_title(label, fontsize=11, fontweight='bold')
                single_ax.set_xlabel('Voltage (V)', color='green', fontsize=10)
                self.multicolor_ylabel(single_ax, ('Current', 'in nA'), ('r', 'r'),
                                        axis='y', size=10, weight='bold', xx=-0.12, yy=0.5)
                single_ax.set_yscale(self.scale_var.get() or 'linear')
                single_ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
                single_ax.legend(loc='lower right', fontsize=8)
                single_fig.tight_layout()
                single_fig.savefig(base + "_IV_Single.png")
                plt.close(single_fig)
            except Exception as e:
                print(f"[Batch] Single-plot save error: {e}")

            print(f"[Batch] Saved → {directory}")
        except Exception as e:
            print(f"[Batch] Auto-save error: {e}")

    # ------------------------------------------
    # 5.7 POST PROCESSING TAB
    # ------------------------------------------
    # ----------------------------------------------------------------------
    # SECTION: POST-PROCESS TAB LAYOUT
    # PURPOSE: Builds the "Post Process" tab: a left control panel
    #          (file selection, voltage-range sliders, display toggles,
    #          action buttons) and a right-hand plot area for reviewing
    #          previously saved (or live) CSV result files.
    # ----------------------------------------------------------------------
    def _setup_post_process_tab(self):
        try:
            import screeninfo
            screen = screeninfo.get_monitors()[0]
            width = screen.width
            height = screen.height
        except:
            width = 1920
            height = 1080

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

        # ------------------------------------------------------------------
        # SUBSECTION: LEFT CONTROL PANEL (scrollable)
        # PURPOSE: Scrollable container for all Post-Process controls
        #          (File / Voltage / Display / Actions sections below).
        # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # SUBSECTION: RIGHT PLOT AREA
        # PURPOSE: Container where post_plot() renders the selected CSV
        #          data.
        # ------------------------------------------------------------------
        post_right = Frame(self.tab3, bg=COLORS["bg_right"])
        post_right.pack(side=Tk.RIGHT, fill=Tk.BOTH, expand=True)

        header_frame = Frame(post_left_inner, bg=COLORS["header"], height=30)
        header_frame.pack(side=Tk.TOP, fill=Tk.X)
        header_frame.pack_propagate(False)

        Label(header_frame, text="Controls Panel",
            bg=COLORS["header"], fg="white",
            font=("Arial", 10, "bold")).pack(pady=6)

        Monitor_frame3 = Frame(post_left_inner)
        Monitor_frame3.pack(side=Tk.TOP, anchor="nw", fill=Tk.X)
        Monitor_frame3.grid_rowconfigure(0, weight=1)

        # ------------------------------------------------------------------
        # SUBSECTION: FILE SECTION
        # PURPOSE: Pick the result CSV file to load/plot.
        # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # SUBSECTION: VOLTAGE SECTION
        # PURPOSE: Start/End voltage-range sliders, range text box, and
        #          plot-title entry.
        # ------------------------------------------------------------------
        voltage_section = Frame(Monitor_frame3, bg=COLORS["panel"], bd=1, relief="groove")
        voltage_section.pack(fill=Tk.X, pady=(0, 4), padx=1)

        Label(voltage_section, text="Voltage", bg=COLORS["panel"], fg=COLORS["warning"],
            font=("Arial", 10, "bold")).pack(anchor="w", padx=4, pady=(2, 1))

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

        # ------------------------------------------------------------------
        # SUBSECTION: DISPLAY SECTION
        # PURPOSE: Log-scale, Temp/Humidity, Breakdown-Voltage, Geiger and
        #          DCR display toggles, plus current-unit selection.
        # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # SUBSECTION: ACTIONS SECTION
        # PURPOSE: Save plot, load Live Data, and Exit buttons.
        # ------------------------------------------------------------------
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

        self.image_label_vi = None

        if getattr(sys, 'frozen', False):
            base_path = sys._MEIPASS
        else:
            base_path = os.path.abspath(".")

        image_path_k = os.path.join(base_path, 'light_files', 'keithley.png')
        try:
            image_keithley = Image.open(image_path_k)
            resized_image_keithley = image_keithley.resize((500, 215))
            photo_keithley = ImageTk.PhotoImage(resized_image_keithley)
            self.image_label_keithley = ttk.Label(self.keithley_img_frame, image=photo_keithley, style='Panel.TLabel')
            self.image_label_keithley.image = photo_keithley
            self.image_label_keithley.pack(anchor='center', expand=True)
        except:
             self.image_label_keithley = ttk.Label(self.keithley_img_frame, text="Keithley Device Ready", font=('Segoe UI', 24, 'bold'), foreground='#BDC3C7')
             self.image_label_keithley.pack(anchor='center', expand=True)

        # ------------------------------------------------------------------
        # SUBSECTION: KEITHLEY IMAGE / PLOT PLACEHOLDER
        # PURPOSE: Placeholder frame stacked under the Measurement tab's
        #          device image, later used by post_plot() to draw the
        #          rendered figure canvas on top.
        # ------------------------------------------------------------------
        self.plot_frame = Frame(post_right, bg="white")
        self.plot_frame.place(relx=0, rely=0, relwidth=1, relheight=1)
        self.plot_frame.lower()

    def post_plot(self, log_file, voltage_start=None, voltage_end=None):
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

        if self.post_canvas:
            self.post_canvas.get_tk_widget().destroy()
            self.post_canvas = None

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

        scale_factor = self.CURRENT_SCALE[self.current_unit_var.get()]
        current = (current_nA * 1e-9) * scale_factor

        self.voltage_min, self.voltage_max = voltage.min(), voltage.max()

        if voltage_start is not None and voltage_end is not None:
            mask = (voltage >= voltage_start) & (voltage <= voltage_end)
            voltage = voltage[mask]
            current = current[mask]
            current_nA = current_nA[mask]
            temperature = temperature[mask]
            humidity = humidity[mask]
        else:
            self.x_start_var.set(self.voltage_min)
            self.x_end_var.set(self.voltage_max)

        show_bd = self.breakdown_voltage_var.get()
        show_geiger = self.giger_prob_var.get()

        v_bd_fit = None
        popt = None
        perr = np.zeros(6)
        fit_success = False

        if show_bd or show_geiger:
            try:
                v_bd_deriv = find_vbd_derivative(voltage, current_nA)
                popt, fit_success, perr = optimize_fit(voltage, current_nA, v_bd_deriv,
                                                user_params=self.user_fit_params if hasattr(self, 'user_fit_params') else None)

                if fit_success:
                    v_bd_fit = popt[0]
            except Exception as e:
                print(f"Breakdown analysis failed: {e}")
                fit_success = False

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

        star = mpath.Path.unit_regular_star(6)
        circle = mpath.Path.unit_circle()
        cut_star = mpath.Path(vertices=np.concatenate([circle.vertices, star.vertices[::-1, ...]]),
                            codes=np.concatenate([circle.codes, star.codes]))

        ax1.plot(voltage, current, marker=cut_star, color='indigo', markersize=8,
                alpha=0.6, label='Measured Data', linestyle='None')

        if show_bd and fit_success and v_bd_fit is not None:
            ax1.axvline(
                v_bd_fit,
                color='red',
                linestyle='--',
                linewidth=2,
                alpha=0.7,
                label=rf'$Breakdown\,Voltage \, V_{{bd}} = {v_bd_fit:.2f}\,\mathrm{{V}}$'
            )

            v_smooth = np.linspace(min(voltage), min(max(voltage), popt[1]-0.1), 500)
            i_fit_nA = dinu_eq8_model(v_smooth, *popt)
            i_fit_scaled = (i_fit_nA * 1e-9) * scale_factor
            ax1.plot(v_smooth, i_fit_scaled, 'g--', linewidth=2, label='Fit Model')

            y_val_nA = dinu_eq8_model(v_bd_fit, *popt)
            y_val_scaled = (y_val_nA * 1e-9) * scale_factor
            ax1.plot(v_bd_fit, y_val_scaled, 'rx', markersize=12,
                    markeredgewidth=3, label='Breakdown Point')

            overvol=float(self.set_ovv.get())
            y_val_nA_ov = dinu_eq8_model(v_bd_fit+overvol, *popt)* 1e-9* scale_factor
            text_pos=max(voltage)-3 

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
            0.5, 1.18,
            self.set_title.get(),
            transform=ax1.transAxes,
            ha="center",
            va="bottom",
            fontsize=16,
            fontweight="bold",
            clip_on=False
        )

        ax1.grid(True, alpha=0.3, linestyle='--')

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

            formula_txt = r"$P_{Geiger} = 1 - e^{-p(V - V_{bd})}$"
            greiger=ax_prob.text(0.3, 0.8, formula_txt, transform=ax_prob.transAxes,
                        fontsize=12, color='darkblue', fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.9, ec="blue"),
                        ha='right', va='top')
            vbd_err_txt = f" $\\pm$ {perr[0]:.2f}" if perr[0] > 0 else ""
            if self.show_dcr_var.get()==0:
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f}{vbd_err_txt} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}")
            else:
                if abs(self.C_ucell)>0:
                    DCR=popt[3]*1e-9/(self.C_ucell*1e3)
                else: DCR=0
                equation_para = (r"$\bf{Fit\ Parameters:}$" + "\n" + f"Breakdown ($V_{{bd}}$): {popt[0]:.2f}{vbd_err_txt} V\n" + f"Critical ($V_{{cr}}$): {popt[1]:.2f} V\n" + f"Geiger Shape ($p$): {popt[2]:.2f}\n" + f"Amplitude ($A$): {popt[3]:.2e}\n" + f"Leak Slope ($a$): {popt[4]:.2e}\n" + f"Leak Offset ($b$): {popt[5]:.2e}\n"+f"DCR : {DCR:0.3f} kHz")

            equation_latex = (r"$I_{tot} = I_{leak} + I_{aval}$" + "\n" + r"$I_{leak} = \exp(aV + b)$" + "\n" + r"$I_{aval} = A \cdot \Delta V \cdot (1 - e^{-p \Delta V}) \cdot \frac{V_{cr}-V_{bd}}{V_{cr}-V}$" + "\n")

            eqn=ax1.text(0.01, 0.95,equation_latex,transform=ax1.transAxes,va='top',fontsize=13,zorder=5,bbox=dict(boxstyle="round", fc="white", alpha=0.92, ec="green"))

            para=ax1.text(0.35, 0.95, equation_para,transform=ax1.transAxes,va='top',fontsize=13,zorder=4,bbox=dict(boxstyle="round", fc="white", alpha=0.85, ec="#27AE60"))

        if not show_geiger: ax1.set_xlabel("Bias Voltage (V)", fontweight='bold',fontsize=14)
        ax1.grid(True, which='both', linestyle='--', alpha=0.5)
        ax1.legend(
            bbox_to_anchor=(0., 1.00, 1., .102),
            loc='lower left',
            ncol=4,
            mode="expand",
            borderaxespad=0.,
            frameon=False,
            fontsize=14
        )

        self.post_canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.post_canvas.draw()
        self.post_canvas.get_tk_widget().pack(fill='both', expand=True)

    def update_voltage_range_from_sliders(self, val):
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

            data = pd.read_csv(file_path)
            v_min = float(data['VOLTS'].min())
            v_max = float(data['VOLTS'].max())

            self.voltage_min = v_min
            self.voltage_max = v_max

            self.voltage_start_slider.config(from_=v_min, to=v_max)
            self.voltage_end_slider.config(from_=v_min, to=v_max)

            self.x_start_var.set(v_min)
            self.x_end_var.set(v_max)

            self.voltage_range_text.delete("1.0", Tk.END)
            self.voltage_range_text.insert(Tk.END, f"{v_min:.2f}, {v_max:.2f}")

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

        data = pd.read_csv(file_path)
        v_min = float(data['VOLTS'].min())
        v_max = float(data['VOLTS'].max())

        self.voltage_min = v_min
        self.voltage_max = v_max

        self.voltage_start_slider.config(from_=v_min, to=v_max)
        self.voltage_end_slider.config(from_=v_min, to=v_max)

        self.x_start_var.set(v_min)
        self.x_end_var.set(v_max)

        self.voltage_range_text.delete("1.0", Tk.END)
        self.voltage_range_text.insert(Tk.END, f"{v_min:.2f}, {v_max:.2f}")

        self.post_plot(file_path, v_min, v_max)

    def apply_voltage_range_from_text(self, event=None):
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


# ==========================================
# 6. APPLICATION ENTRY POINT
# ==========================================
if __name__ == "__main__":
    app = KeithleyGUI()
    app.run()
