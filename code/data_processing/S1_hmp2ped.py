# python $scripts/get_soybean_map.py $genomefile/219.hmp.txt 219_snp.map 219_snp.ped

import sys

n=0
string1 = []
for line in open(sys.argv[1],"r"):
    if n==0:
        pass
    else:
        l = line.strip().split("\t")
        id=l[0]
        if l[2][0]=="0":
            chr = l[2][-1]
        else:
            chr = l[2]
        pos = l[3]
        string1.append(chr+"\t"+id+"\t"+"0"+"\t"+pos+"\n")
    n+=1

w1 = open(sys.argv[2],"w")
w1.writelines(string1)
w1.close()

m=0
ll = ["R","S","K","M","W","Y"]
YY = {"R":"A/G","S":"G/C","K":"G/T","M":"A/C","W":"A/T","Y":"T/C"}
dir = {}
for line in open(sys.argv[1],"r"):
    l = line.strip().split("\t")
    if m==0:
        list1 = l[11:]
    else:
        for i in range(len(list1)):
            if list1[i] not in dir.keys():
                if l[11+i] not in ll:
                    if l[11+i]=="N":
                        dir[list1[i]] = ["0 0"]
                    else:
                        dir[list1[i]] = [l[11+i]+" "+l[11+i]]
                else:
                    dir[list1[i]] = [YY[l[11+i]].split("/")[0]+" "+YY[l[11+i]].split("/")[1]]
            else:
                if l[11+i] not in ll:
                    if l[11+i]=="N":
                        dir[list1[i]].append("0 0")
                    else:
                        dir[list1[i]].append(l[11+i]+" "+l[11+i])
                else:
                    dir[list1[i]].append(YY[l[11+i]].split("/")[0]+" "+YY[l[11+i]].split("/")[1])
                
    m+=1

sorted_dir = sorted(dir.items(),key=lambda item:item[0],reverse=False)

# print sorted_dir

string2 = []
for k,v in sorted_dir:
    string2.append(k+" "+k+" "+"0 0 0 -9 "+" ".join(v)+"\n")

w2 = open(sys.argv[3],"w")
w2.writelines(string2)
w2.close()
        
    