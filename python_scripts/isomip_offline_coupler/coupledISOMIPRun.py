import numpy
import os
import os.path
from optparse import OptionParser
import subprocess

def getPrevDate(currentDate,coupleHours):
  #daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
  prevDate = numpy.array(currentDate)
  prevDate[3] -= options.coupleHours
  while(prevDate[3] < 0):
    prevDate[2] -= 1
    prevDate[3] += 24
    while(prevDate[2] < 1):
      prevDate[1] -= 1
      prevDate[2] += 30 #daysInMonth[numpy.mod(prevDate[1]-1,12)]
      while(prevDate[1] < 1):
        prevDate[0] -= 1
        prevDate[1] += 12
  return tuple(prevDate)
  
def readCurrentDate():
  filePointer = open('restart_timestamp','r')
  content = filePointer.read().splitlines()
  filePointer.close()
  dateString = content[0].lstrip()
  year = int(dateString[0:4])
  month = int(dateString[5:7])
  day = int(dateString[8:10])
  hour = int(dateString[11:13])
  return (year,month,day,hour)
  
def dateToStepIndex(date):
  totalHours = date[3] + 24*((date[2]-1) + 30*((date[1]-1) + 12*date[0]))
  stepIndex = totalHours/coupleHours
  return stepIndex

def namelistFromTemplate(init):
  infile = open('namelist.ocean_forward.template','r')
  outfile = open('namelist.ocean_forward', 'w')
  replacements = {}
  if(init):
    replacements['@doRestart'] = '.false.'
    replacements['@startTime'] = '"0000-01-01_00:00:00"'
  else:
    replacements['@doRestart'] = '.true.'
    replacements['@startTime'] = '"file"'
  replacements['@coupleTime'] = '"0000_%02i:00:00"'%coupleHours
  
  for line in infile:
      for src, target in replacements.iteritems():
          line = line.replace(src, target)
      outfile.write(line)
  infile.close()
  outfile.close()
  
parser = OptionParser()           
parser.add_option("--folder", type="string", default=".", dest="folder")
parser.add_option("--init", action="store_true", dest="init")
parser.add_option("--coupleHours", type="int", default=6, dest="coupleHours")
parser.add_option("--coupleStepCount", type="int", default=120, dest="coupleStepCount")
parser.add_option("--keepOutputFrequency", type="int", default=40, dest="keepOutputFrequency")
parser.add_option("--mpasCommand", type="string", default="mpirun -n 8 ./ocean_forward_model", dest="mpasCommand")
parser.add_option("--pythonCommand", type="string", default="python recomputeISOMIPMeltFluxes.py", dest="pythonCommand")
parser.add_option("--ncoCommandPrefix", type="string", default="", dest="ncoCommandPrefix")

options, args = parser.parse_args()

os.chdir(options.folder)

coupleHours = options.coupleHours

if(options.init):
  namelistFromTemplate(True)

  args = options.mpasCommand.split()
  status = subprocess.call(args)
  
  if status != 0:
    print "ocean_forward_model failed! Exiting."
    exit(status)
    
  # use the first step of the output coupling step output file
  # to create an output file for all kept steps
  inFileName = "output.0000-01-01_00.00.00.nc"
  outFileName = "output.nc"
  args = options.ncoCommandPrefix.split()
  args.extend(["ncrcat","-O","-d","Time,0",inFileName,outFileName])
  status = subprocess.call(args)
  
  if status != 0:
    print "ncrcat failed! Exiting."
    exit(status)

  exit(0)

makeTemplate = True

for stepIndex in range(options.coupleStepCount):

  if(makeTemplate):
    namelistFromTemplate(False)
    makeTemplate = False
 
  dates = []
  dates.append(readCurrentDate())
  dates.append(getPrevDate(dates[0],coupleHours))
    
  
  outputFile = "output.%04i-%02i-%02i_%02i.00.00.nc"%dates[1]
  stepIndex = dateToStepIndex(dates[0])
  if(numpy.mod(stepIndex,options.keepOutputFrequency) == 0 and os.path.exists(outputFile)):
    # append the last time step from output onto the original output file
    print "adding last time step from", outputFile, "to output.nc"
    args = options.ncoCommandPrefix.split()
    args.extend(["ncrcat","--record_append","-d","Time,1",outputFile,"output.nc"])
    status = subprocess.call(args)
    
    if status != 0:
      print "ncrcat failed! Exiting."
      exit(status)

  if(os.path.exists(outputFile)):
    print "deleting unneeded file", outputFile
    os.remove(outputFile)
    
  # do we need to remove the previous output/restart files?
  stepIndex = dateToStepIndex(dates[1])
  if(numpy.mod(stepIndex,options.keepOutputFrequency) != 0):
    restartFile = "restart.%04i-%02i-%02i_%02i.00.00.nc"%dates[1]
    if(os.path.exists(restartFile)):
      print "deleting unneeded file", restartFile
      os.remove(restartFile)
      
    
  outputFile = "output.nc"
  restartFile = "restart.%04i-%02i-%02i_%02i.00.00.nc"%dates[0]
  for writeFile in [outputFile, restartFile]:
    args = options.pythonCommand.split()
    args.extend([outputFile, writeFile])
    status = subprocess.call(args)
    
    if status != 0:
      print "recomputeISOMIPMeltFluxes.py failed! Exiting."
      exit(status)
  
  args = options.mpasCommand.split()
  status = subprocess.call(args)
  
  if status != 0:
    print "ocean_forward_model failed! Exiting."
    exit(status)
    
