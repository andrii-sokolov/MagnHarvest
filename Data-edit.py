file = []

z = -0.003

for i in range(42):
    file.append(open("./Data/Mechanical_resonance_shift"+format(z,"6f")+".dat","r"))
    z+=0.0001
data_file = []

temp_data = []

for f in file:
    temp_data.append([])
    for line in f:
        temp_data[-1].append(line.split())

z = 0.0042

for t in temp_data:
    for l in t:
        data_file.append([z,float(l[0]),-float(l[5])])
    z-=0.0001
    data_file.append([" "," "," "])

out_file = open("3ddata.dat","w")

for l in data_file:
    out_file.write(str(l[0])+"\t"+str(l[1])+"\t"+str(l[2])+"\n")

out_file.close()

out_2 = open("2ddata.dat","w")

z = 0.0042
for i in range(len(temp_data)):
    out_2.write(str(z)+"\t"+str(-float(temp_data[i][22][5]))+"\n")
    z-=0.0001

out_2.close()
