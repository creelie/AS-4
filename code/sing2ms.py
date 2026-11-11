import sys,re
for name in sys.argv[1:]:
    s=open(name+".sing").read()
    n=int(re.search(r"x\(1\.\.(\d+)\)",s).group(1))
    body=s.split("ideal I=",1)[1].split(";",1)[0]
    body=re.sub(r"x\((\d+)\)",r"x\1",body)
    for ch,fn in (("0",name+"_0.ms"),):
        open(fn,"w").write(",".join(f"x{i}" for i in range(1,n+1))+"\n"+ch+"\n"+body.strip()+"\n")
