import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import subprocess
import sys
import os
from PIL import Image, ImageTk  # Requires: pip install pillow

class MasterLauncher:
    def __init__(self, root):
        self.root = root
        self.root.title("SINP Lab Control Center")
        self.root.geometry("500x450")  # Increased height slightly for the exit button
        
        # Configure Styles
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure('Main.TFrame', background='#2C3E50')
        self.style.configure('Header.TLabel', background='#2C3E50', foreground='white', font=('Segoe UI', 16, 'bold'))
        self.style.configure('Action.TButton', font=('Segoe UI', 12, 'bold'), padding=10)
        # Style for the Exit button (Red)
        self.style.configure('Exit.TButton', font=('Segoe UI', 12, 'bold'), padding=10, background='#E74C3C', foreground='white')
        self.style.map('Exit.TButton', background=[('active', '#C0392B')])
        
        self.main_frame = ttk.Frame(self.root, style='Main.TFrame')
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        # Header
        header = ttk.Label(self.main_frame, text="MASTER CONTROL DASHBOARD", style='Header.TLabel')
        header.pack(pady=(40, 30))

        # Buttons Container
        btn_frame = ttk.Frame(self.main_frame, style='Main.TFrame')
        btn_frame.pack(fill=tk.X, padx=50)

        # Button 1: Keithley Controller
        self.btn_keithley = ttk.Button(btn_frame, text="Launch Keithley 2410 Controller", 
                                     style='Action.TButton', 
                                     command=self.launch_keithley)
        self.btn_keithley.pack(fill=tk.X, pady=10)

        # Button 2: Oscilloscope Controller
        self.btn_scope = ttk.Button(btn_frame, text="Launch TekTronix Oscilloscope", 
                                  style='Action.TButton', 
                                  command=self.launch_scope)
        self.btn_scope.pack(fill=tk.X, pady=10)

        # --- NEW EXIT BUTTON ---
        self.btn_exit = ttk.Button(btn_frame, text="EXIT ALL APP", 
                                 style='Exit.TButton', 
                                 command=self.on_close)
        self.btn_exit.pack(fill=tk.X, pady=(30, 10))

        # Status Bar
        self.status_var = tk.StringVar(value="System Ready")
        self.status_bar = tk.Label(self.main_frame, textvariable=self.status_var, 
                                 bd=1, relief=tk.SUNKEN, anchor=tk.W, bg="#34495E", fg="white")
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

        # Keep track of processes
        self.processes = []

    def launch_app(self, script_name):
        """Helper to launch a python script safely."""
        try:
            # Check if file exists
            if not os.path.exists(script_name):
                messagebox.showerror("Error", f"File not found: {script_name}\nMake sure it is in the same folder.")
                return

            # Determine python interpreter (uses the same one running this launcher)
            interpreter = sys.executable
            
            # Launch as a subprocess
            self.status_var.set(f"Launching {script_name}...")
            # Use 'pythonw' logic on windows if needed, or just interpreter
            proc = subprocess.Popen([interpreter, script_name])
            self.processes.append(proc)
            
        except Exception as e:
            messagebox.showerror("Launch Error", f"Failed to launch {script_name}:\n{str(e)}")
            self.status_var.set("Error launching application")

    def launch_keithley(self):
        # MAKE SURE THIS MATCHES YOUR FILENAME EXACTLY
        self.launch_app("modules/keythley_2400_PS_Fit.py") 

    def launch_scope(self):
        # MAKE SURE THIS MATCHES YOUR FILENAME EXACTLY
        self.launch_app("modules/TexTronix_OSC.py") 

    def on_close(self):
        """Cleanup child processes when master closes."""
        # Loop through all launched processes and terminate them
        for proc in self.processes:
            try:
                if proc.poll() is None:  # Check if the process is still running
                    proc.terminate()     # Send termination signal
            except Exception:
                pass # Process might have already closed
        
        self.root.destroy() # Close the main launcher window

if __name__ == "__main__":
    root = tk.Tk()
    app = MasterLauncher(root)
    # Ensure clicking the 'X' button also runs the cleanup logic
    root.protocol("WM_DELETE_WINDOW", app.on_close)
    root.mainloop()
