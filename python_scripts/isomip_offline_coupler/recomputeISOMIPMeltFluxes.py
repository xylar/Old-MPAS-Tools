import numpy
from netCDF4 import Dataset
from optparse import OptionParser

parser = OptionParser()           

options, args = parser.parse_args()


print args

fileName = args[0]
ncFile = Dataset(fileName,'r')
timeIndex = len(ncFile.dimensions['Time'])-1
boundaryLayerPressure =  ncFile.variables['pressure'][timeIndex,:,0]
ncFile.close()

fileName = args[1]
ncFile = Dataset(fileName,'r+')

# From Holland and Jenkins 1999
dT_dS = -5.73e-2
dT_dp=-7.53e-8
T0=9.39e-2

# From Hunter (2006)
gammaT=1e-4
#gammaS=1e-3*gammaT
Cp = 3974.0
L = 334000.0

timeIndex = len(ncFile.dimensions['Time'])-1

temperatureFluxVar = ncFile.variables['surfaceTemperatureFlux']
massFluxVar = ncFile.variables['surfaceMassFlux']
salinityFluxVar = ncFile.variables['surfaceSalinityFlux']

seaSurfacePressure = ncFile.variables['seaSurfacePressure'][timeIndex,:]
iceMask = numpy.array(seaSurfacePressure > 0.0,float)
boundaryLayerSalinity =  ncFile.variables['salinity'][timeIndex,:,0]
boundaryLayerTemperature =  ncFile.variables['temperature'][timeIndex,:,0]
#surfPressure =  ncFile.variables['pressure'][timeIndex,:,0]

interfaceTemperature = T0 + dT_dS*boundaryLayerSalinity \
  + dT_dp*boundaryLayerPressure

# using (3) and (4) from Hunter (2006) 
# or (7) from Jenkins et al. (2001) if gamma constant 
# and no heat flux into ice
# mass flux is in m/s
massFlux = -iceMask*gammaT*(Cp/L)*(interfaceTemperature-boundaryLayerTemperature)


# interface salinity would found using (6) in Jenkins et al. (2001)
# interfaceSalinity = gammaS*boundaryLayerSalinity/(massFlux+gammaS)
# However, this simly leads to the (salinity*thickness) flux equal to zero,
# consistent with no source of salt from the ice.

# gammaS drops out because we are not using interfaceSalinity in the
# computation of interfaceTemperature


# Using (13) from Jenkins et al. (2001)
# temp flux is in deg C*m/s
temperatureFlux = massFlux*interfaceTemperature + gammaT*(interfaceTemperature-boundaryLayerTemperature)
# salinity flux is in PSU*m/s
salinityFlux = 0.0
# The following is true, but trivial 
# salinityFlux = massFlux*interfaceSalinity + gammaS*(interfaceSalinity-boundaryLayerSalinity)

print "min/max fw flux:", numpy.amin(massFlux), numpy.amax(massFlux)
print "min/max T flux: ", numpy.amin(temperatureFlux), numpy.amax(temperatureFlux)
print "min/max S flux: ", numpy.amin(salinityFlux), numpy.amax(salinityFlux)

temperatureFluxVar[timeIndex,:] = temperatureFlux
massFluxVar[timeIndex,:] = massFlux
salinityFluxVar[timeIndex,:] = salinityFlux

ncFile.close()