#!/usr/bin/python

import glob, gzip, csv, sys, os, copy
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

dir         = None
reg         = {'chr':None,'beg':0,'end':(1<<32)-1}
min_length  = 0
min_markers = 0
min_qual    = 0
interactive = False
if len(sys.argv) < 2: 
    print 'Usage: plot.py [OPTIONS] <dir>'
    print 'Options:'
    print '   -i, --interactive                 Run interactively'
    print '   -l, --min-length <num>            Filter input regions shorter than this [0]'
    print '   -n, --min-markers <num>           Filter input regions with fewer marker than this [0]'
    print '   -q, --min-qual <num>              Filter input regions with quality smaller than this [0]'
    print '   -r, --region <chr|chr:beg-end>    Plot this chromosome/region only'
    sys.exit(1)
args = sys.argv[1:]
while len(args):
    if args[0]=='-r' or args[0]=='--region': 
        args = args[1:]
        x = args[0].split(':')
        reg['chr'] = x[0]
        if len(x)>1:
            (reg['beg'],reg['end']) = x[1].split('-')
            reg['beg'] = float(reg['beg'])
            reg['end'] = float(reg['end'])
    elif args[0]=='-i' or args[0]=='--interactive': 
        interactive = True
    elif args[0]=='-l' or args[0]=='--min-length': 
        args = args[1:]
        min_length = float(args[0])
    elif args[0]=='-n' or args[0]=='--min-markers': 
        args = args[1:]
        min_markers = float(args[0])
    elif args[0]=='-q' or args[0]=='--min-qual': 
        args = args[1:]
        min_qual = float(args[0])
    else:
        dir = args[0]
    args = args[1:]

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


smpl2y = {}
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
            if reg['chr']!=None and (chr!=reg['chr'] or pos<reg['beg'] or pos>reg['end']): continue
            smpl = row[3]
            gt   = row[4]
            x = gt.split('/')
            dsg = 2
            if x[0]!=x[1]: dsg = 1
            elif x[0]=='0': dsg = 0
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
            chr   = row[2]
            beg   = int(row[3])
            end   = int(row[4])
            length= int(row[5])
            nmark = int(row[6])
            qual  = float(row[7])
            if length < min_length: continue
            if nmark < min_markers : continue
            if qual < min_qual : continue
            if reg['chr']!=None and (chr!=reg['chr'] or end<reg['beg'] or beg>reg['end']): continue
            if chr not in dat_rg: dat_rg[chr] = {}
            if smpl not in dat_rg[chr]: dat_rg[chr][smpl] = []
            if reg['chr']!=None:
                if beg<reg['beg']: beg = reg['beg']
                if end>reg['end']: end = reg['end']
            dat_rg[chr][smpl].append([beg,end])

off_list = []
off_hash = {}
off = 0
off_sep = 0
dat_rg1 = {}
for chr in chrs:
    if chr in dat_rg:
        rg1 = merge_regions(dat_rg[chr])
        if len(rg1)!=0: dat_rg1[chr] = rg1
    off_hash[chr] = off
    max_pos = 0
    for smpl in dat_gt[chr]:
        if max_pos < dat_gt[chr][smpl][-1][0]: max_pos = dat_gt[chr][smpl][-1][0]
    if off_sep==0: off_sep = max_pos*0.1
    off += max_pos + off_sep
    off_list.append(off)

height = len(fnames)
if len(fnames)>5: heigth = 5
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

if interactive:
    fig, ax1 = plt.subplots(1, 1, figsize=wh, num=dir)
    ax1.yaxis.set_ticks_position('none')
    ax1.format_coord = format_coord
    xtick_lbl = []
    xtick_pos = []
    for chr in dat_gt:
        off  = off_hash[chr]
        xtick_lbl.append(chr)
        xtick_pos.append(off)
        icol = 0
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
            icol += 1
            if icol >= len(cols): 0
    ytick_lbl = []
    ytick_pos = []
    for chr in dat_gt:
        for smpl in dat_gt[chr]:
            ytick_lbl.append(smpl)
            ytick_pos.append(3*smpl2y[smpl]+1)
        break
    ax1.set_xticks(xtick_pos)
    ax1.set_xticklabels(xtick_lbl)
    ax1.set_yticks(ytick_pos)
    ax1.set_yticklabels(ytick_lbl)
    ax1.set_ylim(0,3*len(smpl2y)+0.5)
    plt.subplots_adjust(bottom=0.18,left=0.05,right=0.98)
    plt.show()
else:
    for chr in dat_gt:
        fig, ax1 = plt.subplots(1, 1, figsize=wh)
        ax1.yaxis.set_ticks_position('none')
        ax1.format_coord = format_coord
        tick_lbl = []
        tick_pos = []
        for smpl in dat_gt[chr]:
            y = smpl2y[smpl]
            tick_lbl.append(smpl)
            tick_pos.append(3*y+1)
            if chr in dat_rg and smpl in dat_rg[chr]:
                for rg in dat_rg[chr][smpl]:
                    rect = patches.Rectangle((rg[0],3*y+0.5), rg[1]-rg[0]+1, 2, color='#dddddd')
                    ax1.add_patch(rect)
            if chr in dat_rg and smpl in dat_rg1[chr]:
                for rg in dat_rg1[chr][smpl]:
                    rect = patches.Rectangle((rg[0],3*y+0.5), rg[1]-rg[0]+1, 2, color='#d9534f')
                    ax1.add_patch(rect)
            ax1.plot([x[0] for x in dat_gt[chr][smpl]],[x[1]+3*y for x in dat_gt[chr][smpl]],'.')
        ax1.set_yticks(tick_pos)
        ax1.set_yticklabels(tick_lbl)
        ax1.set_xlabel('chr'+chr)
        ax1.set_ylim(0,3*len(smpl2y)+0.5)
        plt.subplots_adjust(bottom=0.18,left=0.1,right=0.95)
        plt.savefig('rmme-chr'+chr+'.png',dpi=150)
        plt.close()


