from abaqus import session
from odbAccess import *
from abaqusConstants import *
from textRepr import *

coordFile = open('coordsA58.txt','w+')
timeFile = open('timeA58.txt','w+')

odbFileName = 'sph_axiM_RAW4T1_58.odb'

a = openOdb(odbFileName)
# Get the nodeSet
b = a.rootAssembly.instances['PART-1-1'].nodeSets['EDGESET']
c = a.steps['Step-1']

#print d.values[2].data[0]
for val in c.frames:

	C01 = val.fieldOutputs['COORD'].getSubset(region=b).values[0].data[0] 
    
	#-----------------------------------------
	coordFile.write( '%14.7e\n' % (C01))
	# coordFile.write('\n')

    
for val in c.frames:

    T01 = val.frameValue
    
	#-----------------------------------------
    timeFile.write( '%14.7e\n' % (T01))
    # timeFile.write('\n')


coordFile.close()
timeFile.close()

print('Job Done')
	
#d = c.frames[-1].fieldOutputs['U']
#e = d.values.data
