import os
import PyInstaller
print "Creating exe"
os.system("pyinstaller.exe --onefile --windowed \
--icon=icons/confetti.ico \
--exclude-module=matplotlib \
--exclude-module=tkinter \
--exclude-module=zmq \
--exclude-module=twisted \
ConfettiGUI.py")
print "All Done"