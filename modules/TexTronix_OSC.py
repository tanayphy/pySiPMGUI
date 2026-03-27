'''
 MIT License

Copyright (c) 2024 Tanay Dey

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import tkinter as Tk
from tkinter import Toplevel, Frame, Button, Label, Entry, Radiobutton, StringVar, BooleanVar, Checkbutton, GROOVE, VERTICAL, HORIZONTAL, RIGHT, LEFT, BOTTOM, TOP, BOTH, X, Y, Scrollbar, ttk, messagebox as msg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import mplhep as hep
import numpy as np
import math
import time
import random
import threading
import queue
import re
import os
import csv
from datetime import datetime
from itertools import zip_longest
import pyvisa

# --- GLOBAL CONFIGURATION ---
try:
    plt.style.use(hep.style.ROOT)
except:
    plt.style.use('default')

# --- 1. MOCK SCOPE CLASS ---
class MockScope:
    def write(self, cmd): 
        if "SAVe:WAVEform" in cmd:
            time.sleep(0.1) 
        pass
        
    def query(self, cmd):
        time.sleep(0.002) 
        if "*OPC?" in cmd: return '1'
        if "VALue?" in cmd:
            if "MEAS1" in cmd: 
                return str(random.choice([0.2, 0.5, 1.2, 1.5, 2.0]))
            return str(random.gauss(10, 2))
        return "0"

# --- MAIN GUI CLASS ---
class OscilloscopeGUI:
    def __init__(self):
        # Initialize Main Window
        self.window = Tk.Tk()
        self.window.title("TekTronix Controller - Dr. Tanay Dey")
        # Default size, but we'll center it
        self.center_window(self.window, 400, 750)
        
        # Global Variables (Converted to Instance Attributes)
        self.scope = None        
        self.data_queue = queue.Queue()
        
        # Build the GUI
        self.setup_gui()
        
        # Protocol Handlers
        self.window.protocol("WM_DELETE_WINDOW", self.exits)

    def run(self):
        self.window.mainloop()

    # --- HELPER: CENTER WINDOW ---
    def center_window(self, win, width, height):
        """Centers a window on the screen."""
        screen_width = win.winfo_screenwidth()
        screen_height = win.winfo_screenheight()
        x = int((screen_width / 2) - (width / 2))
        y = int((screen_height / 2) - (height / 2))
        win.geometry(f'{width}x{height}+{x}+{y}')

    # --- STYLE CONFIGURATION ---
    def _configure_styles(self):
        """Defines professional styles for the GUI."""
        self.style = ttk.Style()
        self.style.theme_use('clam')

        # Color Palette
        self.colors = {
            'bg_main': '#F4F6F9',        # Light Gray/Blueish
            'bg_sidebar': '#2C3E50',     # Dark Blue/Gray
            'fg_sidebar': '#ECF0F1',     # White text
            'accent': '#2980B9',         # Professional Blue
            'accent_hover': '#3498DB',
            'warning': '#E74C3C',        # Red
            'success': '#27AE60',        # Green
            'text_dark': '#2C3E50'
        }

        # Frame Styles
        self.style.configure('Main.TFrame', background=self.colors['bg_main'])
        self.style.configure('Sidebar.TFrame', background=self.colors['bg_sidebar'])
        
        # Label Styles
        self.style.configure('Sidebar.TLabel', background=self.colors['bg_sidebar'], foreground=self.colors['fg_sidebar'], font=('Segoe UI', 10))
        self.style.configure('Header.TLabel', background=self.colors['bg_sidebar'], foreground=self.colors['fg_sidebar'], font=('Segoe UI', 12, 'bold'))
        self.style.configure('Panel.TLabel', background=self.colors['bg_main'], foreground=self.colors['text_dark'], font=('Segoe UI', 10))
        
        # Labelframe Styles
        self.style.configure('Group.TLabelframe', background=self.colors['bg_sidebar'], foreground=self.colors['fg_sidebar'])
        self.style.configure('Group.TLabelframe.Label', background=self.colors['bg_sidebar'], foreground='#BDC3C7', font=('Segoe UI', 9, 'bold'))
        
        self.style.configure('PanelGroup.TLabelframe', background=self.colors['bg_main'], foreground=self.colors['text_dark'])
        self.style.configure('PanelGroup.TLabelframe.Label', background=self.colors['bg_main'], foreground='#7f8c8d', font=('Segoe UI', 9, 'bold'))

        # Button Styles
        self.style.configure('Action.TButton', font=('Segoe UI', 10, 'bold'), padding=5)
        self.style.map('Action.TButton', background=[('active', self.colors['accent_hover'])], foreground=[('active', 'black')])

        # Treeview
        self.style.configure("Treeview", font=('Segoe UI', 9), rowheight=25)
        self.style.configure("Treeview.Heading", font=('Segoe UI', 9, 'bold'))


    # --- 6. SETUP MAIN GUI ---
# --- 6. SETUP MAIN GUI ---
    def setup_gui(self):
        # Apply your professional styles
        self._configure_styles()
        
        # Create a main container frame
        main_frame = ttk.Frame(self.window, style='Main.TFrame')
        main_frame.pack(fill=BOTH, expand=True)
        
        # Title Label
        lbl_title = ttk.Label(
            main_frame, 
            text="TekTronix Controller", 
            font=('Segoe UI', 18, 'bold'), 
            background=self.colors['bg_main'], 
            foreground=self.colors['accent']
        )
        lbl_title.pack(pady=(150, 10))
        
        # Under Construction Warning
        lbl_status = ttk.Label(
            main_frame, 
            text="🚧 App is Under Construction 🚧", 
            font=('Segoe UI', 14, 'bold'), 
            background=self.colors['bg_main'], 
            foreground=self.colors['warning']
        )
        lbl_status.pack(pady=10)
        
        # Subtitle/Description
        lbl_desc = ttk.Label(
            main_frame, 
            text="Core features are currently being developed.\nPlease check back later for the full release.", 
            font=('Segoe UI', 10), 
            justify=Tk.CENTER, 
            background=self.colors['bg_main'],
            foreground=self.colors['text_dark']
        )
        lbl_desc.pack(pady=20)
        
        # Exit Button
        btn_exit = ttk.Button(
            main_frame, 
            text="Exit Application", 
            style='Action.TButton', 
            command=self.exits
        )
        btn_exit.pack(pady=30)        
    def exits(self, event=None):
        try: 
            if self.stop_event: self.stop_event.set()
        except: pass
        self.window.destroy()
        os._exit(0)

if __name__ == "__main__":
    app = OscilloscopeGUI()
    app.run()
