import numpy
import os
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
parser.add_option("--folder", type="string", default="isomip_spherical", dest="folder")
parser.add_option("--init", action="store_true", dest="init")
parser.add_option("--coupleHours", type="int", default=6, dest="coupleHours")
parser.add_option("--coupleStepCount", type="int", default=10, dest="coupleStepCount")

options, args = parser.parse_args()

codePath = os.getcwd()

os.chdir(options.folder)

coupleHours = options.coupleHours

init = options.init
makeTemplate = True

for stepIndex in range(options.coupleStepCount):

  if(makeTemplate):
    namelistFromTemplate(init)
    if(not init):
      makeTemplate = False
 
  args = ["mpirun","-n","8","./ocean_forward_model"]
  status = subprocess.call(args)
  
  if status != 0:
    print "ocean_forward_model failed! Exiting."
    exit(status)
    
  currentDate = readCurrentDate()
  prevDate = getPrevDate(currentDate,coupleHours)
  outputFile = "output.%04i-%02i-%02i_%02i.00.00.nc"%prevDate
  restartFile = "restart.%04i-%02i-%02i_%02i.00.00.nc"%currentDate
  
  for writeFile in [outputFile, restartFile]:
    args = ["python", "%s/recomputeISOMIPMeltFluxes.py"%codePath,
            outputFile,writeFile]
    status = subprocess.call(args)
    
    if status != 0:
      print "recomputeISOMIPMeltFluxes.py failed! Exiting."
      exit(status)
  
  init = False

