#include <Wire.h>

const int SLAVE_ADDRESS = 8; // CHANGE THIS FOR S2, S3, etc.

volatile char cmdBuffer[32]; 
volatile bool newData = false;

int minPin = 2;
int maxPin = 53;

// Holds the result of the most recent "Read [PIN]" command, sent back
// to the Master the next time it does Wire.requestFrom() on us.
// Defaults to 0 until a Read command has actually been processed.
volatile byte lastReadValue = 0;

void setup() {
  Serial.begin(9600);           // Optional: plug this board into a PC via
                                 // its own USB port (separate from the I2C
                                 // bus) to watch exactly what it receives
                                 // and writes -- handy for confirming
                                 // whether a given pin's command actually
                                 // arrives/executes, independent of the
                                 // Master board.
  Wire.begin(SLAVE_ADDRESS);    
  Wire.onReceive(receiveEvent); 
  Wire.onRequest(requestEvent);
}

void loop() {
  if (newData) {
    processCommand(String((char*)cmdBuffer));
    newData = false; 
  }
}

void receiveEvent(int howMany) {
  int i = 0;
  while (Wire.available() && i < 31) {
    cmdBuffer[i] = Wire.read();
    i++;
  }
  cmdBuffer[i] = '\0'; 
  newData = true;
}

// Called automatically by the Wire library when the Master does
// Wire.requestFrom(SLAVE_ADDRESS, 1) after sending a "Read [PIN]"
// command. Just hands back whatever processCommand() last measured.
void requestEvent() {
  Wire.write(lastReadValue);
}

void processCommand(String input) {
  input.trim();
  String inputUpper = input;
  inputUpper.toUpperCase();

  // ── UPDATE / U : set a pin ON or OFF ────────────────────────────────
  String args = "";
  if (inputUpper.startsWith("UPDATE ")) {
    args = input.substring(7);
  } else if (inputUpper.startsWith("U ")) {
    args = input.substring(2);
  }
  
  // If args is no longer empty, we know we received a valid update command
  if (args.length() > 0) {
    args.trim();
    int spaceIndex = args.indexOf(' ');
    
    if (spaceIndex > 0) {
      int pin = args.substring(0, spaceIndex).toInt();
      int state = args.substring(spaceIndex + 1).toInt();
      
      if (pin >= minPin && pin <= maxPin) {
        pinMode(pin, OUTPUT);
        digitalWrite(pin, state > 0 ? HIGH : LOW);
        Serial.print("UPDATE pin "); Serial.print(pin);
        Serial.print(" -> "); Serial.println(state > 0 ? "HIGH" : "LOW");
      } else {
        Serial.print("UPDATE pin "); Serial.print(pin);
        Serial.println(" REJECTED: out of range");
      }
    }
    return;
  }

  // ── READ / R : measure a pin's current ON/OFF state ─────────────────
  // Stores the result in lastReadValue; the Master picks it up via a
  // separate I2C requestFrom() -> requestEvent() round trip, since a
  // slave can't push data to the Master on its own.
  String readArgs = "";
  if (inputUpper.startsWith("READ ")) {
    readArgs = input.substring(5);
  } else if (inputUpper.startsWith("R ")) {
    readArgs = input.substring(2);
  }

  if (readArgs.length() > 0) {
    readArgs.trim();
    int pin = readArgs.toInt();

    if (pin >= minPin && pin <= maxPin) {
      // Relay pins are driven as OUTPUT; digitalRead() on an OUTPUT pin
      // reads back the output latch, i.e. the last commanded state.
      pinMode(pin, OUTPUT);
      lastReadValue = digitalRead(pin) ? 1 : 0;
      Serial.print("READ pin "); Serial.print(pin);
      Serial.print(" -> "); Serial.println(lastReadValue);
    } else {
      lastReadValue = 0;
      Serial.print("READ pin "); Serial.print(pin);
      Serial.println(" REJECTED: out of range, reporting 0");
    }
  }
}
