
# PySiPMGUI

**A Universal Python-Based Software for Photodetector I-V Quality Assurance**

PySiPMGUI is an open-source, non-blocking Graphical User Interface (GUI) developed in Python 3. It is designed for the automated I-V characterization of Silicon Photomultipliers (SiPMs), providing a free and platform-independent alternative to commercial software for the particle and astroparticle physics communities.

## 🚀 Key Features
* **Automated I-V Sweeps:** Controls PyVISA-compatible source-measure units (e.g., Keithley 2400 series).
* **Equipment Safety:** Features recursive safe voltage ramping and hardware/software current compliance interlocks.
* **Environmental Logging:** Interfaces with Arduino (e.g., SHT30 sensor) to log real-time temperature and humidity.
* **Physics-Based Analysis:** Fits I-V curves to extract Breakdown Voltage (V_BD), Geiger probability, and estimates the Dark Count Rate (DCR).
* **Responsive GUI:** Built with Tkinter, allowing real-time plotting and instant PAUSE/STOP controls.

## ⚙️ Installation

Python 3 is required. Install the necessary dependencies using pip:

```bash
pip install -r requirements.txt
pip install pyvisa pyserial numpy scipy matplotlib pwlf
