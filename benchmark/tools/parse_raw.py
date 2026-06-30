import sys, re
tag, ws, path = sys.argv[1], sys.argv[2], sys.argv[3]
lines = open(path).readlines()
done = [i for i,l in enumerate(lines) if "grid-done" in l]
if   len(done) >= 2: start, end = done[-2]+1, done[-1]
elif len(done) == 1: start, end = 0, done[-1]
else:                start, end = 0, len(lines)
door=arm=setn=cn=None; rows=0
for line in lines[start:end]:
    mc = re.search(r"%cell\s+%([a-z0-9-]+)\s+%([a-z0-9-]+)\s+(\d+)\s+([\d.]+)", line)
    if mc: door,arm,setn,cn = mc.group(1),mc.group(2),mc.group(3),mc.group(4).replace('.',''); continue
    mt = re.search(r"took\s+\S+/([\d.]+)", line)
    if mt and door is not None:
        print(f"{ws},{tag},{door},{arm},{cn},{setn},{mt.group(1).replace('.','')},ok"); rows+=1; door=None
sys.stderr.write(f"  {tag}: {rows} rows\n")
