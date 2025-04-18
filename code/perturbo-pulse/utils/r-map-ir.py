ifc_index_file = "rset.dat"
def magnitude(vec):
    sum = 0
    for i in vec:
        sum += i**2
    return sum

# read all R vectors
R = set()
R1 = []
R2 = [] 
with open(ifc_index_file, 'r') as lines:
    for line in lines:

        line = line.rstrip()
        line_data = line.split()
        
        rvec1 = tuple([float(x) for x in line_data[1:4]])
        rvec2 = tuple([float(x) for x in line_data[4:]])
        R1.append(rvec1)
        R2.append(rvec2)
        
        if rvec1 not in R:
            R.add(rvec1)
        if rvec2 not in R:
            R.add(rvec2)

# sort by magnitude of R
R = list(R)
R.sort(key=lambda x: magnitude(x))

print("Size of R grid = " + str(len(R)))

#output the index function 
file = open('rset_index.dat','w')
for i in range(len(R1)):
    file.write(str(R.index(R1[i]))+'    '+str(R.index(R2[i]))+'\n' )
file.close()

#output the R function 
file = open('rset_irr.dat','w')
file.write(str(len(R))+'\n')
for r in R:
    file.write(str(r[0])+'    '+str(r[1])+'    '+str(r[2])+'\n' )
file.close()
