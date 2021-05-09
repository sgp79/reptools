import os
import reptools

def count_dirs(csvfile,paireddirs=(),unpaireddirs=(),pairsuffixes=('_1','_2'),filetype=None,overwrite=False,basename=False):
	
	if not overwrite and os.path.exists(csvfile):
		raise IOError('Target file already exists. To overwrite, set overwrite=True.')
	
	filetypes = reptools.select_filetypes(filetype)
	counts={}
	for dir in paireddirs:
		filelist = [os.path.join(dir,fn) for fn in os.listdir(dir) if os.path.splitext(fn)[1] in filetypes]
		counts[dir] = {}
		for fn in filelist:
			root = os.path.splitext(os.path.split(fn)[1])[0][:-len(pairsuffixes[0])]
			if root not in list(counts.keys()):
				counts[dir][root]=0
			if os.path.splitext(fn)[1] in ['.fastq','.fq']:
				counts[dir][root] += reptools.fastqcounter(fn)
			else:
				counts[dir][root] += reptools.fascounter(fn)
	for dir in unpaireddirs:
		counts[dir] = {}
		filelist = [os.path.join(dir,fn) for fn in os.listdir(dir) if os.path.splitext(fn)[1] in filetypes]
		for fn in filelist:
			root = os.path.splitext(os.path.split(fn)[1])[0]
			if os.path.splitext(fn)[1] in ['.fastq','.fq']:
				counts[dir][root] = reptools.fastqcounter(fn)
			else:
				counts[dir][root] = reptools.fascounter(fn)
	#
	allroots = sorted(set([k for dict in counts for k in counts[dict] ]))
	
	if basename:
		titles = [os.path.split(path)[1] for path in paireddirs] + [os.path.split(path)[1] for path in unpaireddirs]
	else:
		titles = paireddirs+unpaireddirs
	
	with open(csvfile,'w') as out_handle:
		out_handle.write('root,%s\n' % (','.join(titles) ))
		for root in allroots:
			out_handle.write('%s,%s\n' % (root,','.join([str(counts[dir][root]) if root in counts[dir] else '0' for dir in paireddirs+unpaireddirs] ) ) )

