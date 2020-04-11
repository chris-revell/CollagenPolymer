import numpy as np
import os
from sys import argv


data           = np.genfromtxt("output/output.txt",delimiter=", ")
conditions     = np.genfromtxt("output/conditions.txt")
Nmonomers      = 2*int(conditions[0,1])
tmax           = conditions[1,1]
outputInterval = conditions[2,1]

for i in range(int(tmax/outputInterval)):
    os.system("clear")
    print("Rendering {}/{}".format(i,int(tmax/outputInterval)))
    outfile = open("output/povrayTmp{:03d}.pov".format(i),"w")
    outfile.write("#include \"colors.inc\"\n")
    outfile.write("camera {\n" )
    outfile.write("  sky <0,0,1>           \n")
    outfile.write("  direction <-1,0,0>      \n")
    outfile.write("  right <-4/3,0,0>      \n")
    outfile.write("  location <0,-8,18> \n" )
    outfile.write("  look_at <0,0,0>     \n" )
    outfile.write("  angle 15      \n")
    outfile.write("}\n")
    outfile.write("global_settings { ambient_light White }\n")
    outfile.write("light_source {\n" )
    outfile.write("  <10,-10,20>   \n")
    outfile.write("  color White*2 \n")
    outfile.write("}\n")
    outfile.write("background { color White }\n" )
    subdata = data[i*Nmonomers:(i+1)*Nmonomers,:]
    for j in range(Nmonomers):
        outfile.write("sphere{{<{},{},{}>,{} texture{{pigment{{color Green}}}}}}\n".format(subdata[j,0],subdata[j,1],subdata[j,2],0.1))
    outfile.close()
    os.system("povray output/povrayTmp{:03d}.pov > /dev/null 2>&1".format(i))
    os.system("rm output/povrayTmp{:03d}.pov".format(i))

os.system("convert -delay 10 -loop 0 output/povrayTmp*.png output/animated.gif")
os.system("rm output/*.png")
