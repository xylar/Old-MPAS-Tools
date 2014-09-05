import numpy
from netCDF4 import Dataset
from optparse import OptionParser

parser = OptionParser()           

options, args = parser.parse_args()


print args

fileName = args[0]
ncFile = Dataset(fileName,'r')
timeIndex = len(ncFile.dimensions['Time'])-1
surfPressure =  ncFile.variables['pressure'][timeIndex,:,0]
ncFile.close()

fileName = args[1]
ncFile = Dataset(fileName,'r+')

dT_dS = -5.73e-2
dT_dp=-7.53e-8
T0=9.39e-2
gamma=0.0001
Cp = 3974.0
L = 334000.0

timeIndex = len(ncFile.dimensions['Time'])-1

tempFlux = ncFile.variables['surfaceTemperatureFlux']
massFlux = ncFile.variables['surfaceMassFlux']
saltFlux = ncFile.variables['surfaceSalinityFlux']

seaSurfacePressure = ncFile.variables['seaSurfacePressure'][timeIndex,:]
iceMask = numpy.array(seaSurfacePressure > 0.0,float)
surfSalinity =  ncFile.variables['salinity'][timeIndex,:,0]
surfTemperature =  ncFile.variables['temperature'][timeIndex,:,0]
#surfPressure =  ncFile.variables['pressure'][timeIndex,:,0]

freezingTemperature = T0 + dT_dS*surfSalinity + dT_dp*surfPressure

temperatureFlux = -iceMask*gamma*(surfTemperature-freezingTemperature) 
freshwaterFlux = -temperatureFlux*(Cp/L)
temperatureFlux += freshwaterFlux*freezingTemperature
print numpy.amin(freshwaterFlux), numpy.amax(freshwaterFlux)
print numpy.amin(temperatureFlux), numpy.amax(temperatureFlux)

tempFlux[timeIndex,:] = temperatureFlux
massFlux[timeIndex,:] = freshwaterFlux
saltFlux[timeIndex,:] = 0.0

ncFile.close()