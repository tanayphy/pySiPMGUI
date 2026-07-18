 #include <Wire.h>
 #include "Adafruit_SHT31.h"

 Adafruit_SHT31 sht31;

 void setup() {
 Serial.begin(9600);
 Wire.begin();
 sht31.begin(0x44); // continue even if not found
 }

 void loop() {
 if (!Serial.available()) return;
 // The following lines must not be modified by the user.
if (Serial.readStringUntil('\n') == "all") 
{ float t = sht31.readTemperature();
 float h = sht31.readHumidity();

 if (!isnan(t) && !isnan(h))
 Serial.print(t,2), Serial.print("\t"), Serial.println(h,2);
 else
 Serial.println("-999\t-999");
 }
 }
