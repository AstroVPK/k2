import math as math
import numpy as np
import urllib, urllib2
import os as os
import sys as sys
import subprocess
import argparse
import matplotlib.pyplot as plt
import pdb

try:
	import libcarma as libcarma
except ImportError:
	print 'libcarma is not setup. Setup libcarma by sourcing bin/setup.sh'
	sys.exit(1)

class k2LC(libcarma.basicLC):

	sap = ['sap', 'raw', 'uncal', 'un-cal', 'uncalibrated', 'un-calibrated']
	pdcsap = ['pdcsap', 'mast', 'cal', 'calib', 'calibrated']
	k2sff = ['k2sff', 'vj', 'vanderburg', 'vanderburgjohnson', 'vanderburg-johnson']
	k2sc = ['k2sc', 'aigrain']
	k2varcat = ['k2varcat', 'armstrong']

	def _getCanonicalFileName(self, name, campaign, processing):
		fileName = ''
		if processing in self.sap or processing in self.pdcsap:
			fileName = ''.join(['ktwo', name, '-', campaign, '_llc.dat'])
		elif processing in self.k2sff or processing in self.k2sc or processing in self.k2varcat:
			if processing in self.k2sff:
				fileName = ''.join(['hlsp_k2sff_k2_lightcurve_' ,name, '-', campaign, '_kepler_v1_llc.dat'])
			elif processing in self.k2sc:
				fileName = ''.join(['hlsp_k2sc_k2_llc_', name, '-', campaign, '_kepler_v1_lc.dat'])
			elif processing in self.k2varcat:
				fileName = ''.join(['hlsp_k2varcat_k2_lightcurve_', name, '-', campaign, '_kepler_v2_llc.dat'])
			else:
				raise ValueError('Unrecognized k2LC type')
		else:
			raise ValueError('Unrecognized k2LC type')
		return fileName

	def _getMAST(self, name, campaign, path, goid, gopi):
		baseURL = 'http://archive.stsci.edu/pub/k2/lightcurves'
		recordFile = 'k2List.dat'

		fileName = self._getCanonicalFileName(name, campaign, 'mast')
		fileNameFits = ''.join([fileName[0:-3], 'fits'])
		filePath = os.path.join(path, fileName)
		filePathFits = ''.join([filePath[0:-3], 'fits'])

		if not os.path.isfile(filePathFits):
			recordFilePath = os.path.join(path, recordFile)
			with open(recordFilePath, 'a') as record:
				record.write('%s %s %s %s\n'%(name, campaign, goid, gopi))
			camp = ''.join(['c', str(int(campaign[1:]))])
			name1Dir = ''.join([name[0:4], '00000'])
			name2Dir = ''.join([name[4:6], '000'])
			fullURL = '/'.join([baseURL, camp, name1Dir, name2Dir, fileNameFits])
			result = urllib.urlretrieve(fullURL, filePathFits)
		if not os.path.isfile(filePath):
			subprocess.call(['topcat', '-stilts', 'tcopy',  'in=%s'%(filePathFits), 'ofmt=ascii', 'out=%s'%(filePath)])

	def _getHLSP(self, name, campaign, path):
		baseURL = 'http://archive.stsci.edu/missions/hlsp'

		fileName = self._getCanonicalFileName(name, campaign, 'k2sff')
		fileNameFits = ''.join([fileName[0:-3], 'fits'])
		filePath = os.path.join(path, fileName)
		filePathFits = os.path.join(path, fileNameFits)
		if not os.path.isfile(filePathFits):
			name1Dir = ''.join([name[0:4], '00000'])
			name2Dir = name[4:]
			fullURL = '/'.join([baseURL, 'k2sff', campaign, name1Dir, name2Dir, fileNameFits])
			try:
				ret = urllib2.urlopen(fullURL)
			except urllib2.HTTPError:
				pass
			else:
				result = urllib.urlretrieve(fullURL, filePathFits)
		if not os.path.isfile(filePath) and os.path.isfile(filePathFits):
			subprocess.call(['topcat', '-stilts', 'tcopy',  'in=%s'%(filePathFits), 'ofmt=ascii', 'out=%s'%(filePath)])

		fileName = self._getCanonicalFileName(name, campaign, 'k2sc')
		fileNameFits = ''.join([fileName[0:-3], 'fits'])
		filePath = os.path.join(path, fileName)
		filePathFits = os.path.join(path, fileNameFits)
		if not os.path.isfile(filePathFits):
			name1Dir = ''.join([name[0:4], '00000'])
			fullURL = '/'.join([baseURL, 'k2sc', campaign, name1Dir, fileNameFits])
			try:
				ret = urllib2.urlopen(fullURL)
			except urllib2.HTTPError:
				pass
			else:
				result = urllib.urlretrieve(fullURL, filePathFits)
		if not os.path.isfile(filePath) and os.path.isfile(filePathFits):
			subprocess.call(['topcat', '-stilts', 'tcopy',  'in=%s'%(filePathFits), 'ofmt=ascii', 'out=%s'%(filePath)])

		fileName = self._getCanonicalFileName(name, campaign, 'k2varcat')
		fileNameFits = ''.join([fileName[0:-3], 'fits'])
		filePath = os.path.join(path, fileName)
		filePathFits = os.path.join(path, fileNameFits)
		if not os.path.isfile(filePathFits):
			name1Dir = ''.join([name[0:4], '00000'])
			name2Dir = ''.join([name[4:6], '000'])
			fullURL = '/'.join([baseURL, 'k2varcat', campaign, name1Dir, name2Dir, fileNameFits])
			try:
				ret = urllib2.urlopen(fullURL)
			except urllib2.HTTPError:
				pass
			else:
				result = urllib.urlretrieve(fullURL, filePathFits)
		if not os.path.isfile(filePath) and os.path.isfile(filePathFits):
			subprocess.call(['topcat', '-stilts', 'tcopy',  'in=%s'%(filePathFits), 'ofmt=ascii', 'out=%s'%(filePath)])

	def _readMAST(self, name, campaign, path, processing):
		fileName = self._getCanonicalFileName(name, campaign, processing)
		filePath = os.path.join(path, fileName)
		with open(filePath,'r') as k2File:
			allLines = k2File.readlines()
		self._numCadences = len(allLines) - 1
		startT = -1.0
		lineNum = 1
		while startT == -1.0:
			words = allLines[lineNum].split()
			nextWords = allLines[lineNum + 1].split()
			if words[0] != '""' and nextWords[0] != '""':
				startT = float(words[0])
				dt = float(nextWords[0]) - float(words[0])
			else:
				lineNum += 1
		self.startT = startT
		self._dt = dt ## Increment between epochs.
		self.cadence = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.t = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.x = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.y = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.yerr = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.mask = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E']) ## Numpy array of mask values.
		for i in xrange(self.numCadences):
			words = allLines[i +1].split()
			self.cadence[i] = int(words[2])
			if words[9] == '0':
				self.t[i] = float(words[0]) - self.startT
				if processing in self.sap:
					try:
						self.y[i] = float(words[3])
						self.yerr[i] = float(words[4])
						self.mask[i] = 1.0
					except ValueError:
						self.y[i] = 0.0
						self.yerr[i] = math.sqrt(sys.float_info[0])
						self.mask[i] = 0.0
				elif processing  in self.pdcsap:
					try:
						self.y[i] = float(words[7])
						self.yerr[i] = float(words[8])
						self.mask[i] = 1.0
					except ValueError:
						self.y[i] = 0.0
						self.yerr[i] = math.sqrt(sys.float_info[0])
						self.mask[i] = 0.0
				else:
					raise ValueError('Unrecognized k2LC type')
			else:
				if words[0] != '""':
					self.t[i] = float(words[0]) - self.startT
				else:
					self.t[i] = self.t[i - 1] + self.dt
				self.yerr[i] = math.sqrt(sys.float_info[0])
				self.mask[i] = 0.0
		self._dt = float(np.nanmedian(self.t[1:] - self.t[:-1])) ## Increment between epochs.
		self._T = self.t[-1] - self.t[0] ## Total duration of the light curve.

	def _readK2SFF(self, name, campaign, path, processing):
		fileNameMAST = self._getCanonicalFileName(name, campaign, 'mast')
		filePathMAST = os.path.join(path, fileNameMAST)
		with open(filePathMAST, 'r') as k2FileMAST:
			allLinesMAST = k2FileMAST.readlines()
		self._numCadences = len(allLinesMAST) - 1
		startT = -1.0
		lineNum = 1
		while startT == -1.0:
			words = allLinesMAST[lineNum].split()
			nextWords = allLinesMAST[lineNum + 1].split()
			if words[0] != '""' and nextWords[0] != '""':
				startT = float(words[0])
				dt = float(nextWords[0]) - float(words[0])
			else:
				lineNum += 1
		self.startT = startT
		self._dt = dt ## Increment between epochs.
		self.cadence = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.t = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.x = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.y = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.yerr = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.mask = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E']) ## Numpy array of mask values.
		for i in xrange(self.numCadences):
			words = allLinesMAST[i + 1].split()
			self.cadence[i] = int(words[2])
			if words[9] == '0':
				self.t[i] = float(words[0]) - self.startT
			else:
				if words[0] != '""':
					self.t[i] = float(words[0]) - self.startT
				else:
					self.t[i] = self.t[i - 1] + self.dt
		self._dt = float(np.nanmedian(self.t[1:] - self.t[:-1])) ## Increment between epochs.
		self._T = self.t[-1] - self.t[0] ## Total duration of the light curve.

		fileName = self._getCanonicalFileName(name, campaign, 'k2sff')
		filePath = os.path.join(path, fileName)
		with open(filePath,'r') as k2File:
			allLines = k2File.readlines()
		for line in allLines[1:]:
			words = line.rstrip('\n').split()
			cadNum = int(words[5])
			index = np.where(self.cadence == cadNum)[0][0]
			self.t[index] = float(words[0]) - self.startT
			self.y[index] = float(words[2])
			self.mask[index] = 1.0

		valSum = 0.0
		countSum = 0.0
		for i in xrange(self.numCadences - 1):
			valSum += self.mask[i + 1]*self.mask[i]*math.pow((self.y[i + 1] - self.y[i]), 2.0)
			countSum += self.mask[i + 1]*self.mask[i]
		noise = math.sqrt(valSum/countSum)
		for i in xrange(self.numCadences):
			if self.mask[i] == 1.0:
				self.yerr[i] = noise

	def _readK2SC(self, name, campaign, path, processing):
		fileNameMAST = self._getCanonicalFileName(name, campaign, 'mast')
		filePathMAST = os.path.join(path, fileNameMAST)
		with open(filePathMAST, 'r') as k2FileMAST:
			allLinesMAST = k2FileMAST.readlines()
		self._numCadences = len(allLinesMAST) - 1
		startT = -1.0
		lineNum = 1
		while startT == -1.0:
			words = allLinesMAST[lineNum].split()
			nextWords = allLinesMAST[lineNum + 1].split()
			if words[0] != '""' and nextWords[0] != '""':
				startT = float(words[0])
				dt = float(nextWords[0]) - float(words[0])
			else:
				lineNum += 1
		self.startT = startT
		self._dt = dt ## Increment between epochs.
		self.cadence = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.t = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.x = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.y = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.yerr = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.mask = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E']) ## Numpy array of mask values.
		for i in xrange(self.numCadences):
			words = allLinesMAST[i + 1].split()
			self.cadence[i] = int(words[2])
			if words[9] == '0':
				self.t[i] = float(words[0]) - self.startT
			else:
				if words[0] != '""':
					self.t[i] = float(words[0]) - self.startT
				else:
					self.t[i] = self.t[i - 1] + self.dt
		self._dt = float(np.nanmedian(self.t[1:] - self.t[:-1])) ## Increment between epochs.
		self._T = self.t[-1] - self.t[0] ## Total duration of the light curve.

		fileName = self._getCanonicalFileName(name, campaign, 'k2sc')
		filePath = os.path.join(path, fileName)
		with open(filePath,'r') as k2File:
			allLines = k2File.readlines()
		for line in allLines[1:]:
			words = line.rstrip('\n').split()
			if int(words[7]) == 0: 
				time = float(words[0]) - self.startT
				index = np.where(self.t == time)[0][0]
				self.y[index] = float(words[8])
				self.yerr[index] = float(words[6])
				self.mask[index] = 1.0

	def _readK2VARCAT(self, name, campaign, path, processing):
		fileNameMAST = self._getCanonicalFileName(name, campaign, 'mast')
		filePathMAST = os.path.join(path, fileNameMAST)
		with open(filePathMAST, 'r') as k2FileMAST:
			allLinesMAST = k2FileMAST.readlines()
		self._numCadences = len(allLinesMAST) - 1
		startT = -1.0
		lineNum = 1
		while startT == -1.0:
			words = allLinesMAST[lineNum].split()
			nextWords = allLinesMAST[lineNum + 1].split()
			if words[0] != '""' and nextWords[0] != '""':
				startT = float(words[0])
				dt = float(nextWords[0]) - float(words[0])
			else:
				lineNum += 1
		self.startT = startT
		self._dt = dt ## Increment between epochs.
		self.cadence = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.t = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.x = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.y = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.yerr = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E'])
		self.mask = np.require(np.zeros(self.numCadences), requirements=['F', 'A', 'W', 'O', 'E']) ## Numpy array of mask values.
		for i in xrange(self.numCadences):
			words = allLinesMAST[i + 1].split()
			self.cadence[i] = int(words[2])
			if words[9] == '0':
				self.t[i] = float(words[0]) - self.startT
			else:
				if words[0] != '""':
					self.t[i] = float(words[0]) - self.startT
				else:
					self.t[i] = self.t[i - 1] + self.dt
		self._dt = float(np.nanmedian(self.t[1:] - self.t[:-1])) ## Increment between epochs.
		self._T = self.t[-1] - self.t[0] ## Total duration of the light curve.

		fileName = self._getCanonicalFileName(name, campaign, 'k2varcat')
		filePath = os.path.join(path, fileName)
		try:
			k2File = open(filePath,'r')
		except IOError:
			pass
		else:
			allLines = k2File.readlines()
			for line in allLines[1:]:
				words = line.rstrip('\n').split()
				time = float(words[0]) - self.startT
				index = np.where(self.t == time)[0][0]
				self.y[index] = float(words[3])
				self.yerr[index] = float(words[4])
				self.mask[index] = 1.0
			k2File.close()

	def read(self, name, band = None, path = None, **kwargs):
		processing = kwargs.get('processing', 'sap').lower()
		campaign = kwargs.get('campaign', 'c05').lower()
		fileName = self._getCanonicalFileName(name, campaign, processing)
		goid = kwargs.get('goid', '').lower()
		gopi = kwargs.get('gopi', '').lower()
		if path is None:
			try:
				path = os.environ['K2DATADIR']
			except KeyError:
				raise KeyError('Environment variable "K2DATADIR" not set! Please set "K2DATADIR" to where all K2 data should live first...')
		filePath = os.path.join(path, fileName)

		self._computedCadenceNum = -1
		self._tolIR = 1.0e-3
		self._fracIntrinsicVar = 0.0
		self._fracNoiseToSignal = 0.0
		self._maxSigma = 2.0
		self._minTimescale = 2.0
		self._maxTimescale = 0.5
		self._pSim = 0
		self._qSim = 0
		self._pComp = 0
		self._qComp = 0
		self._isSmoothed = False ## Has the LC been smoothed?
		self._dtSmooth = 0.0
		self._isRegular = True
		self.XSim = np.require(np.zeros(self.pSim), requirements=['F', 'A', 'W', 'O', 'E']) ## State of light curve at last timestamp
		self.PSim = np.require(np.zeros(self.pSim*self.pSim), requirements=['F', 'A', 'W', 'O', 'E']) ## Uncertainty in state of light curve at last timestamp.
		self.XComp = np.require(np.zeros(self.pComp), requirements=['F', 'A', 'W', 'O', 'E']) ## State of light curve at last timestamp
		self.PComp = np.require(np.zeros(self.pComp*self.pComp), requirements=['F', 'A', 'W', 'O', 'E']) ## Uncertainty in state of light curve at last timestamp.
		self._name = str(name) ## The name of the light curve (usually the object's name).
		self._band = str(r'Kep') ## The name of the photometric band (eg. HSC-I or SDSS-g etc..).
		self._xunit = r'$d$' ## Unit in which time is measured (eg. s, sec, seconds etc...).
		self._yunit = r'who the f*** knows?' ## Unit in which the flux is measured (eg Wm^{-2} etc...).

		self._getMAST(name, campaign, path, goid, gopi)
		self._getHLSP(name, campaign, path)

		if processing in self.sap or processing in self.pdcsap:
			self._readMAST(name, campaign, path, processing)
		elif processing in self.k2sff:
			self._readK2SFF(name, campaign, path, processing)
		elif processing in self.k2sc:
			self._readK2SC(name, campaign, path, processing)
		elif processing in self.k2varcat:
			self._readK2VARCAT(name, campaign, path, processing)

		count = int(np.sum(self.mask))
		y_meanSum = 0.0
		yerr_meanSum = 0.0
		for i in xrange(self.numCadences):
			y_meanSum += self.mask[i]*self.y[i]
			yerr_meanSum += self.mask[i]*self.yerr[i]
		if count > 0.0:
			self._mean = y_meanSum/count
			self._meanerr = yerr_meanSum/count
		else:
			self._mean = 0.0
			self._meanerr = 0.0
		y_stdSum = 0.0
		yerr_stdSum = 0.0
		for i in xrange(self.numCadences):
			y_stdSum += math.pow(self.mask[i]*self.y[i] - self._mean, 2.0)
			yerr_stdSum += math.pow(self.mask[i]*self.yerr[i] - self._meanerr, 2.0)
		if count > 0.0:
			self._std = math.sqrt(y_stdSum/count)
			self._stderr = math.sqrt(yerr_stdSum/count)
		else:
			self._std = 0.0
			self._stderr = 0.0

	def write(self, name, path = None, **kwrags):
		pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-id', '--ID', type = str, default = '205905563', help = r'EPIC ID')
	parser.add_argument('-p', '--processing', type = str, default = 'sap', help = r'sap/pdcsap/k2sff/k2sc/k2varcat etc...')
	parser.add_argument('-c', '--campaign', type = str, default = 'c03', help = r'Campaign')
	parser.add_argument('-goid', '--goID', type = str, default = '', help = r'Guest Observer ID')
	parser.add_argument('-gopi', '--goPI', type = str, default = '', help = r'Guest Observer PI')
	args = parser.parse_args()

	LC = k2LC(name = args.ID, band = 'Kep', processing = args.processing, campaign = args.campaign, goid = args.goID, gopi = args.goPI)

	LC.plot()
	LC.plotacf()
	LC.plotsf()
	plt.show(False)