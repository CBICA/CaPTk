#! /usr/bin/env python
import os
import sys
import ctypes
from  ctypes import CFUNCTYPE, c_int,c_char,c_char_p
import numpy
# Update 7/2/2020: Modified ConfettiCore.py to accommodate nuitka packaging. 
class ConfettiApi():
    def __init__(self):
        try:
            self.libraryName = os.path.normpath(self.getCoreModulePath())
            self.CoreModule = ctypes.cdll.LoadLibrary(self.libraryName)
            print "Successfully loaded the ConfettiCore binary."
            self.listner= None
            self.txtCallback_ref = None
        except Exception as e:
            self.CoreModule= None
            print("Failed to load ConfettiCore binary. Exception: ")
            print(e)
    def getCoreModulePath(self):
        if os.name is 'nt':
            libraryName = "ConfettiCore.dll"
        else:
            libraryName = "libConfettiCore.so"
        if (getattr(sys, 'frozen', False)):
            application_path = os.path.dirname(sys.executable)
        elif __file__:
            print "Using __file__"
            if os.name is 'nt':
                application_path = self.attemptWindowsLongPathFix(os.path.dirname(__file__))
            else:
                application_path = os.path.dirname(__file__)
        baseDir=os.path.join(application_path, '..' ) # use to search directory above dist
		# match if dll found in either this or one-higher directory
        for root, dirs, files in os.walk(application_path):
            for file in files:
                if libraryName in file:
                    return (os.path.join(root, file))
        for root, dirs, files in os.walk(baseDir):
            for file in files:
                if libraryName in file:
                    return (os.path.join(root, file))

    def attemptWindowsLongPathFix(self, in_path):
        # Do not call this on non-Windows platforms.
        try:
            path = unicode(in_path)
        except Exception as e:
            print(e) 
        GetLongPathName = ctypes.windll.kernel32.GetLongPathNameW
        buffer = ctypes.create_unicode_buffer(GetLongPathName(path, 0, 0))
        GetLongPathName(path, buffer, len(buffer))
        print buffer.value
        return buffer.value

    def registerCppCallback(self):
        if self.CoreModule is None:
            return

        TXT_CALLBACK =CFUNCTYPE(c_int, c_char_p)
        self.CoreModule.registerCallback.restype = None
        self.CoreModule.registerCallback.argtypes = [TXT_CALLBACK]
        self.txtCallback_ref = TXT_CALLBACK(self.pyTxtCallbackFunc)
        self.CoreModule.registerCallback(self.txtCallback_ref)
        self.listner = None
    def pyTxtCallbackFunc(self,txt):
        if self.CoreModule is None:
            return
        if txt is not None and self.listner is not None:
            self.listner("CORE:"+ str(txt))
        else:
            print("CORE:"+ str(txt))
        return 0
    def setListner(self, listner):
        if self.CoreModule is None:
            return
        self.listner= listner
        self.listner("Listening to core")


    def genConnectivitySig(self,tdiCsvFile, signatureFileName,fiberFile):
        if self.CoreModule is None:
            return
        res = self.CoreModule.genConnectivitySig(c_char_p(tdiCsvFile),c_char_p(signatureFileName),c_char_p(fiberFile))

    def genAdaptiveCluster(self,signatureFile, tempRoot, outFile):
        if self.CoreModule is None:
            return
        self.CoreModule.genAdaptiveCluster.restype = c_int
        self.CoreModule.genAdaptiveCluster.argtypes = [c_char_p,c_char_p,c_char_p]
        res = self.CoreModule.genAdaptiveCluster(c_char_p(signatureFile),c_char_p(tempRoot),c_char_p(outFile))

    def extractTracts(self,outDIR, tempRoot, inputIDFile, inputNamesFile, inputLabelFile, inputBfloatFile):
        if self.CoreModule is None:
            return
        self.CoreModule.extractTracts.restype = c_int
        self.CoreModule.extractTracts.argtypes = [c_char_p,c_char_p,c_char_p,c_char_p,c_char_p,c_char_p]
        res = self.CoreModule.extractTracts(c_char_p(outDIR),c_char_p(tempRoot),c_char_p(inputIDFile),
                                            c_char_p(inputNamesFile),c_char_p(inputLabelFile),c_char_p(inputBfloatFile))

    def getValue(self,varName):
        if self.CoreModule is None:
            return
        self.CoreModule.getValue.restype = c_char_p
        self.CoreModule.getValue.argtypes = [c_char_p]
        res = self.CoreModule.getValue(c_char_p(varName))
        return res

