#include <Wire.h>
const int MASTER_ADDRESS = 0;

int minPin = 2;
int maxPin = 53;
void setup() {
  Serial.begin(9600);
  Wire.begin(); 
  
  Serial.println("Master Ready.");
  Serial.println("Formats allowed:");
  Serial.println("  Board [ID] Update [PIN] [STATE]  (e.g., Board 8 Update 13 1)");
  Serial.println("  B [ID] U [PIN] [STATE]           (e.g., B 8 U 13 1)");
  Serial.println("  Board [ID] Read [PIN]            (e.g., Board 8 Read 13)");
  Serial.println("  B [ID] R [PIN]                   (e.g., B 8 R 13)");
}

void loop() {
  if (Serial.available() > 0) {
    String input = Serial.readStringUntil('\n');
    input.trim();

    if (input.length() > 0) {
      String inputUpper = input;
      inputUpper.toUpperCase();

      // Chop off "BOARD " or "B " if the user typed it
      if (inputUpper.startsWith("BOARD ")) {
        input = input.substring(6);
      } else if (inputUpper.startsWith("B ")) {
        input = input.substring(2);
      }
      
      input.trim(); // Clean up any remaining spaces

      // Now the string just looks like "8 UPDATE 13 1" or "8 U 13 1"
      // or "8 READ 13" or "8 R 13"
      int spaceIndex = input.indexOf(' ');
      
      if (spaceIndex > 0) {
        int slaveAddress = input.substring(0, spaceIndex).toInt();
        String command = input.substring(spaceIndex + 1);

        Serial.println(command);

        // Figure out up front whether this is a Read request, since a
        // Read needs a follow-up I2C requestFrom() to fetch the pin
        // state back from a remote (non-Master) board.
        String cmdUpper = command;
        cmdUpper.trim();
        cmdUpper.toUpperCase();
        bool isReadCmd = cmdUpper.startsWith("READ ") || cmdUpper.startsWith("R ");

        if (slaveAddress==MASTER_ADDRESS)
        {
          processCommand(command);
          }
        else
        {
        Wire.beginTransmission(slaveAddress);
        Wire.write(command.c_str()); 
        byte error = Wire.endTransmission();
        
        if (error == 0) {
          Serial.print("Success: Sent '");
          Serial.print(command);
          Serial.print("' to Board ");
          Serial.println(slaveAddress);

          if (isReadCmd) {
            // Give the slave's main loop a moment to process the
            // buffered command and update its pin-state byte before we
            // ask for it back.
            delay(15);  // give the slave's loop() + any relay/level-shifter
                        // settling time a bit more margin before we ask
                        // for the pin state back
            byte reply = 0;
            byte n = Wire.requestFrom(slaveAddress, (uint8_t)1);
            if (n >= 1) {
              reply = Wire.read();
              Serial.print("Pin State (Board ");
              Serial.print(slaveAddress);
              Serial.print("): ");
              Serial.println(reply);
            } else {
              Serial.print("Error: Board ");
              Serial.print(slaveAddress);
              Serial.println(" did not reply to Read request.");
            }
          }
        } else {
          Serial.print("Error: Board ");
          Serial.print(slaveAddress);
          Serial.println(" not responding. Check wiring.");
        }
      }
      } else {
        Serial.println("Error: Invalid format.");
      }
    }
  }
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
      }
    }
    return;
  }

  // ── READ / R : report a pin's current ON/OFF state over Serial ─────
  // Only reached for the local board (MASTER_ADDRESS); remote boards
  // report their state back over I2C via requestEvent() in Slave.ino.
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
      int state = digitalRead(pin);
      Serial.print("Pin ");
      Serial.print(pin);
      Serial.print(" State: ");
      Serial.println(state);
    } else {
      // Deliberately avoid ending this line in a stray '0' or '1' digit
      // -- the PC side just scans for the last 0/1 digit in the reply
      // to decide ON/OFF, so an out-of-range message must not end in
      // one of those digits.
      Serial.print("Pin ");
      Serial.print(pin);
      Serial.println(" State: OUT_OF_RANGE");
    }
  }
}
