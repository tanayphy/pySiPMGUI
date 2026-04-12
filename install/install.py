import subprocess
import sys

def install_requirements():
    """
    Installs the required external libraries using pip.
    Built-in Python libraries are omitted as they require no installation.
    """
    # List of required external packages
    packages = [
        "pyvisa",      # for import pyvisa
        "matplotlib",  # for import matplotlib
        "numpy",       # for import numpy
        "requests",    # for import requests
        "Pillow",      # for from PIL import Image
        "screeninfo",  # for import screeninfo
        "pyserial",    # for import serial (Package is pyserial, not serial)
        "pandas",      # for import pandas
        "scipy"        # for from scipy.optimize import curve_fit
    ]

    print("=========================================")
    print("Initializing SINP Lab Dependency Setup...")
    print("=========================================\n")

    for package in packages:
        try:
            print(f"Attempting to install: {package}")
            # Run pip install securely using the current python executable
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
            print(f"✅ Successfully installed {package}\n")
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to install {package}. Error: {e}\n")
        except Exception as e:
            print(f"❌ An unexpected error occurred while installing {package}: {e}\n")

    print("=========================================")
    print("Setup complete! You can now run your main application.")
    print("=========================================")

if __name__ == "__main__":
    install_requirements()