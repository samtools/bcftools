#!/usr/bin/python

import glob, gzip, csv, sys, os, copy, re
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

def usage(msg=None):
    if msg==None:
        print 'Usage: plot.py [OPTIONS] <dir>'
        print 'Options:'
        print '   -H, --highlight +group1,-group2       Highlight calls shared within group1 but not present in group2'
        print '   -i, --interactive                     Run interactively'
        print '   -l, --min-length <num>                Filter input regions shorter than this [0]'
        print '   -n, --min-markers <num>               Filter input regions with fewer marker than this [0]'
        print '   -o, --outfile <file>                  Output file name [plot.png]'
        print '   -q, --min-qual <num>                  Filter input regions with quality smaller than this [0]'
        print '   -r, --region [^]<chr|chr:beg-end>     Plot this chromosome/region only'
        print '   -s, --samples <file>                  List of samples to show, rename or group: "name[\\tnew_name[\\tgroup]]"'
        print '   -h, --help                            This usage text'
        print 'Matplotlib options:'
        print '   +adj, --adjust <str>     Set plot adjust [bottom=0.18,left=0.07,right=0.98]'
        print '   +dpi, --dpi <num>        Set bitmap DPI [150]'
        print '   +sxt, --show-xticks      Show x-ticks (genomic coordinate)'
        print '   +xlb, --xlabel <str>     Set x-label'
        print '   +xli, --xlimit <num>     Extend x-range by this fraction [0.05]'
    else:
        print msg
    sys.exit(1)

dir         = None
regs        = None
min_length  = 0
min_markers = 0
min_qual    = 0
interactive = False
sample_file = None
highlight   = None
outfile     = None
adjust      = 'bottom=0.18,left=0.07,right=0.98'
dpi         = 150
xlim        = 0.05
show_xticks = False
xlabel      = None

if len(sys.argv) < 2: usage()
args = sys.argv[1:]
while len(args):
    if args[0]=='-r' or args[0]=='--region': 
        args = args[1:]
        regs = args[0]
    elif args[0]=='-i' or args[0]=='--interactive': 
        interactive = True
    elif args[0]=='-l' or args[0]=='--min-length': 
        args = args[1:]
        min_length = float(args[0])
    elif args[0]=='-n' or args[0]=='--min-markers': 
        args = args[1:]
        min_markers = float(args[0])
    elif args[0]=='-o' or args[0]=='--outfile': 
        args = args[1:]
        outfile = args[0]
    elif args[0]=='-q' or args[0]=='--min-qual': 
        args = args[1:]
        min_qual = float(args[0])
    elif args[0]=='-H' or args[0]=='--highlight': 
        args = args[1:]
        highlight = args[0]
    elif args[0]=='-s' or args[0]=='--samples': 
        args = args[1:]
        sample_file = args[0]
    elif args[0]=='-?' or args[0]=='-h' or args[0]=='--help':
        usage()
    elif args[0]=='+adj' or args[0]=='--adjust': 
        args = args[1:]
        adjust = args[0]
    elif args[0]=='+dpi' or args[0]=='--dpi': 
        args = args[1:]
        dpi = float(args[0])
    elif args[0]=='+xlb' or args[0]=='--xlabel': 
        args = args[1:]
        xlabel = args[0]
    elif args[0]=='+sxt' or args[0]=='--show-xticks': 
        show_xticks = True
    elif args[0]=='+xli' or args[0]=='--xlimit': 
        args = args[1:]
        xlim = float(args[0])
    else:
        dir = args[0]
    args = args[1:]

if interactive and outfile!=None: usage("Use -i, --interactive or -o, --outfile, but not both")
if not interactive and outfile==None: outfile = 'plot.png'

def wrap_hash(**args): return args
adjust = eval("wrap_hash("+adjust+")")


import matplotlib as mpl
for gui in ['TKAgg','GTKAgg','Qt4Agg','WXAgg','MacOSX']:
    try:
        mpl.use(gui,warn=False, force=True)
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        break
    except:
        continue

cols = [ '#337ab7', '#5cb85c', '#5bc0de', '#f0ad4e', '#d9534f', 'grey', 'black' ]
mpl.rcParams['axes.color_cycle'] = cols

globstr = os.path.join(dir, '*.txt.gz')
fnames = glob.glob(globstr)
if len(fnames)==0: usage("No data files found in \""+dir+"\"")

def parse_regions(str):
    if str==None: return None
    regs = { 'inc':[], 'exc':[] }
    list = str.split(',')
    key = 'inc'
    if list[0][0]=='^': 
        key = 'exc'
        list[0] = list[0][1:]
    for reg in list:
        x = reg.split(':')
        chr = x[0]
        beg = 0
        end = (1<<32)-1
        if len(x)>1:
            (beg,end) = x[1].split('-')
            beg = float(beg)
            end = float(end)
        regs[key].append({'chr':chr,'beg':beg,'end':end})
    return regs

def region_overlap(regs,chr,beg,end):
    if regs==None: return (beg,end)
    if len(regs['exc'])>0:
        for reg in regs['exc']:
            if chr==reg['chr']: return None
        return (beg,end)
    if len(regs['inc'])==0: return (beg,end)
    for reg in regs['inc']:
        if chr!=reg['chr']: continue
        if beg>reg['end']: continue
        if end<reg['beg']: continue
        if beg<reg['beg']: beg = reg['beg']
        if end>reg['end']: end = reg['end']
        return (beg,end)
    return None

def parse_outfile(fname):
    files = re.split(r',',fname)
    bname = re.search(r'^(.+)\.[^.]+$', files[0]).group(1)
    for i in range(len(files)-1):
        files[i+1] = bname+"."+files[i+1]
    return files

def next_region(rgs):
    min = None
    for smpl in rgs:
        if len(rgs[smpl])==0: continue
        reg = rgs[smpl][0]
        if min==None:
            min = [0,0]
            min[0] = reg[0]
            min[1] = reg[1]
        if min[0] > reg[0]: min[0] = reg[0]
    if min==None: return None
    for smpl in rgs:
        if len(rgs[smpl])==0: continue
        reg = rgs[smpl][0]
        if min[1] > reg[1]: min[1] = reg[1]
        if min[1] > reg[0] - 1 and min[0] != reg[0]: min[1] = reg[0] - 1
    return min;

def merge_regions(rg):
    rgs = copy.deepcopy(rg)
    out = {}
    while True:
        min = next_region(rgs)
        if min==None: break
        beg = min[0]
        end = min[1]
        smpls = []
        for smpl in rgs:
            if len(rgs[smpl])==0: continue
            reg = rgs[smpl][0]
            if reg[0] > end: continue
            if reg[1] > end: 
                rgs[smpl][0][0] = end + 1
            else:
                rgs[smpl] = rgs[smpl][1:]
            if smpl not in out: out[smpl] = []
            smpls.append(smpl)
        if len(smpls)>1:
            for smpl in smpls: out[smpl].append([beg,end])
    return out

def prune_regions(groups,regions):
    regs = {'+':{},'-':{}}
    for smpl in regions:
        grp = groups[smpl]
        for reg in regions[smpl]:
            key = str(reg[0])+"-"+str(reg[1])   # reg=[beg,end] -> "beg-end"
            if key not in regs[grp]: regs[grp][key] = 0
            regs[grp][key] += 1
    nexp = 0
    for smpl in groups:
        if groups[smpl]=='+': nexp += 1
    for smpl in regions:
        rm = []
        for reg in regions[smpl]:
            key = str(reg[0])+"-"+str(reg[1]) 
            if key in regs['-']: rm.append(reg)
            elif key not in regs['+'] or regs['+'][key]!=nexp: rm.append(reg)
        for reg in rm:
            if reg in regions[smpl]:
                regions[smpl].remove(reg)
    return regions

def parse_samples(fname,highlight):
    if fname==None: return (None,None,{})
    samples = {}
    groups  = {}
    grp2sgn = {}
    smpl2y  = {}
    # parse "+name" to create a map "name":"+"
    if highlight!=None:
        for grp in re.split(r',', highlight):
            if grp[0]!='+' and grp[0]!='-': usage("Expected + or - before the group name: "+grp)
            grp2sgn[grp[1:]] = grp[0]
    # read samples, renaming them
    with open(fname) as f:
        for line in f:
            row  = re.split(r'\s+', line.rstrip('\n'))
            smpl = row[0]
            if len(row)==1: samples[smpl] = smpl
            else:
                samples[smpl] = row[1]
            if len(row)==3:
                grp = row[2]
                if grp in grp2sgn:
                    grp = grp2sgn[grp]
                else:
                    grp = '+'
                groups[smpl] = grp
            y = len(smpl2y)
            smpl2y[smpl] = y
    if highlight==None: groups = None
    return (samples,groups,smpl2y)

regs = parse_regions(regs)
(samples,groups,smpl2y) = parse_samples(sample_file,highlight)

dat_gt = {}
dat_rg = {}
chrs   = []
for fname in fnames:
    f = gzip.open(fname, 'rb')
    reader = csv.reader(f, 'tab')
    for row in reader:
        if row[0]=='GT':
            chr  = row[1]
            pos  = int(row[2])
            reg  = region_overlap(regs,chr,pos,pos)
            if reg==None: continue
            for i in range(3,len(row),2):
                smpl = row[i]
                if samples!=None and smpl not in samples: continue
                gt   = row[i+1]
                x = gt.split('/')
                if x[0]=='.': continue          # missing genotype ./.
                dsg = 2
                if x[0]!=x[1]: dsg = 1
                elif x[0]=='0': continue        # skip HomRef 0/0 genotypes
                if chr not in dat_gt: 
                    dat_gt[chr] = {}
                    chrs.append(chr)
                if smpl not in dat_gt[chr]: 
                    dat_gt[chr][smpl] = []
                if smpl not in smpl2y:
                    y = len(smpl2y)
                    smpl2y[smpl] = y
                dat_gt[chr][smpl].append([pos,dsg])
        elif row[0]=='RG':
            smpl  = row[1]
            if samples!=None and smpl not in samples: continue
            chr   = row[2]
            beg   = int(row[3])
            end   = int(row[4])
            length= int(row[5])
            nmark = int(row[6])
            qual  = float(row[7])
            if length < min_length: continue
            if nmark < min_markers : continue
            if qual < min_qual : continue
            reg = region_overlap(regs,chr,beg,end)
            if chr not in dat_rg: dat_rg[chr] = {}
            if smpl not in dat_rg[chr]: dat_rg[chr][smpl] = []
            if reg!=None:
                if beg<reg[0]: beg = reg[0]
                if end>reg[1]: end = reg[1]
            dat_rg[chr][smpl].append([beg,end])

if samples==None: 
    samples = {}
    for smpl in smpl2y: samples[smpl] = smpl

# list the samples in the same order as encountered in the file, from top to bottom
for smpl in smpl2y:
    smpl2y[smpl] = len(smpl2y) - smpl2y[smpl] - 1

off_list = []
off_hash = {}
off = 0
off_sep = 0
dat_rg1 = {}
for chr in chrs:
    if chr in dat_rg:
        rg1 = merge_regions(dat_rg[chr])
        if groups!=None: 
            rg1 = prune_regions(groups,rg1)
        if len(rg1)!=0: dat_rg1[chr] = rg1
    off_hash[chr] = off
    max_pos = 0
    for smpl in dat_gt[chr]:
        if max_pos < dat_gt[chr][smpl][-1][0]: max_pos = dat_gt[chr][smpl][-1][0]
    if off_sep==0: off_sep = max_pos*0.1
    off += max_pos + off_sep
    off_list.append(off)

height = len(smpl2y)
if len(smpl2y)>5: heigth = 5
wh = 20,height

def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

def format_coord(x, y):
    chr = None
    off = 0
    for i in range(len(off_list)):
        chr = chrs[i]
        if off_list[i] > x: break
        off = off_list[i]
    return 'chr%s:%s'%(chr,bignum(int(x - off)))


fig, ax1 = plt.subplots(1, 1, figsize=wh, num=dir)
ax1.yaxis.set_ticks_position('none')
ax1.format_coord = format_coord
xtick_lbl = []
xtick_pos = []
max_x = 0
for chr in dat_gt:
    off  = off_hash[chr]
    icol = 0
    max  = 0
    for smpl in dat_gt[chr]:
        y = smpl2y[smpl]
        if chr in dat_rg and smpl in dat_rg[chr]:
            for rg in dat_rg[chr][smpl]:
                rect = patches.Rectangle((rg[0]+off,3*y+0.5), rg[1]-rg[0]+1, 2, color='#dddddd')
                ax1.add_patch(rect)
        if chr in dat_rg1 and smpl in dat_rg1[chr]:
            for rg in dat_rg1[chr][smpl]:
                rect = patches.Rectangle((rg[0]+off,3*y+0.5), rg[1]-rg[0]+1, 2, color='#d9534f')
                ax1.add_patch(rect)
        ax1.plot([x[0]+off for x in dat_gt[chr][smpl]],[x[1]+3*y for x in dat_gt[chr][smpl]],'.',color=cols[icol])
        if max_x < dat_gt[chr][smpl][-1][0]+off: max_x = dat_gt[chr][smpl][-1][0]+off
        if max < dat_gt[chr][smpl][-1][0]: max = dat_gt[chr][smpl][-1][0]
        icol += 1
        if icol >= len(cols): 0
    xtick_lbl.append(chr)
    xtick_pos.append(off)
ytick_lbl = []
ytick_pos = []
for chr in dat_gt:
    for smpl in dat_gt[chr]:
        ytick_lbl.append(samples[smpl])
        ytick_pos.append(3*smpl2y[smpl]+1)
    break
if xlim!=0:
    ax1.set_xlim(0,max_x+xlim*max_x)
lbl_pos = 3*(len(smpl2y)-1)
ax1.annotate('   HomAlt ',xy=(max_x,lbl_pos-1),xycoords='data',va='center')
ax1.annotate('   Het',xy=(max_x,lbl_pos-2),xycoords='data',va='center')
if not show_xticks:
    ax1.set_xticks(xtick_pos)
    ax1.set_xticklabels(xtick_lbl)
if xlabel!=None:
    ax1.set_xlabel(xlabel)
ax1.set_yticks(ytick_pos)
ax1.set_yticklabels(ytick_lbl)
ax1.set_ylim(0,3*len(smpl2y)+0.5)
plt.subplots_adjust(**adjust)
if interactive:
    plt.show()
else:
    files = parse_outfile(outfile)
    for file in (parse_outfile(outfile)):
        plt.savefig(file,dpi=dpi)
    plt.close()


