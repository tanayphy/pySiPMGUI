
######################################################################################################################################################################################################
Copyright (c) 2024 Tanay Dey
Creative Commons Attribution-NonCommercial 4.0

Title:
PySiPMGUI: A Universal Python-Based Software for Photodetector I-V Quality Assurance: From Underground Dark Matter Searches to Astroparticle Cherenkov Cameras

Author: Dr. Tanay Dey

Co Authors: Suraj Shaw, Ritabrata Banerjee, Pratik Majumdar, Satyaki
Bhattacharya

https://doi.org/10.48550/arXiv.2603.24781

https://github.com/tanayphy/pySiPMGUI.git

#####################################################################################################################################################################################################

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
python3 install/install.py 

or

pip install -r requirements.txt
 
```



# PySiPMGUI User Manual

Welcome to PySiPMGUI! This tool allows you to safely and automatically characterize Silicon Photomultipliers (SiPMs) by measuring their Current-Voltage (I-V) curves, monitoring environmental conditions, and extracting key parameters like Breakdown Voltage (**V_BD**) and Dark Count Rate (**DCR**).

## 1. Hardware Setup
Before launching the software, ensure your physical setup is properly connected:
1. **Place the SiPM** inside a light-tight dark box.
2. **Connect the SMU** (e.g., Keithley 2400 series) to the SiPM biasing circuit.
3. **Connect the SMU to the PC** via your USB-to-GPIB or direct USB cable.
4. *(Optional)* **Connect the Arduino** (with the SHT30 sensor inside the dark box) to your PC via USB.

## 2. Launching the Software
Open your terminal or command prompt, navigate to the folder containing your script, and run:
`python main.py`
*(Ensure all dependencies from your `requirements.txt` are installed before running).*

## 3. Interface Overview
The GUI is divided into two main sections:
* **Left Control Panel:** Contains all the input fields for connecting devices, setting voltage parameters, and triggering analysis.
* **Center Plotting Area:** Displays the live I-V curve and environmental data as the test runs. It includes tabs for **Measurement** and **Post Process/Analysis Result**.

## 4. Running an I-V Measurement
Follow these steps to safely bias your SiPM and capture data:

1. **Connect the SMU:**
   * Find the **Instrument Connection** section.
   * Enter the correct VISA address for your power supply (e.g., `GPIB0::24::INSTR`).
   * Click **Connect**.
2. **Set Sweep Parameters:**
   * **Polarity:** Select Positive or Negative biasing.
   * **Start / End:** Enter your starting and ending voltages (e.g., 0 to 30 V).
   * **Step:** Define the voltage increment (e.g., 0.2 V).
   * **Delay:** Set the wait time between applying voltage and measuring current (e.g., 0.5 s) to allow the circuit to stabilize.
   * **Compliance:** Set a safe maximum current limit (e.g., 10000 nA) to prevent damaging the SiPM.
3. **Connect Environment Sensor (Optional):**
   * Check the **Connect Arduino** box.
   * Select the appropriate COM port from the dropdown menu to log temperature and humidity.
4. **Start the Scan:**
   * Click the green **START TEST** button.
   * The software will begin safely ramping the voltage. You will see the data points appear on the graph in real-time.
   * *Emergency controls:* You can click **PAUSE** to temporarily halt the sweep, or **STOP** to immediately trigger a safe ramp-down to 0 V.

## 5. Analyzing the Data
Once the measurement finishes (or if you load a previously saved `.csv` file), you can analyze the I-V curve to extract the SiPM's characteristics:

1. Go to the **Analysis** section on the left panel.
2. Check the boxes for what you want to calculate: **Breakdown V**, **Geiger Prob**, or **DCR**.
3. Click **Fit Params**. 
   * *Note:* If you leave the parameter fields blank, the software will automatically estimate the initial breakdown voltage by finding the maximum derivative of the log-current curve.
4. Switch to the **Analysis Result** tab in the center area to view the fitted curve, the calculated breakdown voltage, and the estimated Dark Count Rate.

## 6. Saving Your Work
* **Data Logging:** During the test, the software holds the data in memory. Once finished, you will be prompted to save the results.
* **Formats:** Data is saved as a `.csv` file, and you can also save images of your plots directly from the GUI for your records.

## 7. Troubleshooting
* **"Instrument Connection Failed":** Verify that your SMU is powered on, the cables are secure, and the correct VISA address is entered.
* **Test stops abruptly:** You likely hit the current compliance limit. The software has initiated a safe ramp-down to protect your device. Try increasing the limit slightly if you are sure the SiPM can handle it, or check for light leaks in your dark box.
* **Arduino not reading:** Ensure the correct COM port is selected and no other program (like the Arduino IDE Serial Monitor) is currently using that port.
