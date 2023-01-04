import reptools
import subprocess

def remove_adapter_bbduk(
                   infile,outfile,adapters,pairfile=False,outpair=False,adapter_rc=False,
                   discard_trimmed=False,bbduk_path='bbduk.sh'
                   ):
    if outpair and not pairfile: raise ValueError
    
    call_list = [bbduk_path]
    call_list.append('-Xmx1g')
    
    if pairfile:
        if not outpair: raise ValueError
        call_list.append('in1={}'.format(infile))
        call_list.append('out1={}'.format(outfile))
        call_list.append('in2={}'.format(pairfile))
        call_list.append('out2={}'.format(outpair))
    else:
        call_list.append('in='.format(infile))
        call_list.append('out='.format(outfile))
    
    call_list.append('literal={}'.format(','.join(adapters)))
    
    if not adapter_rc:
        call_list.append('rcomp=f')
    
    call_list.append('hdist=2')
    call_list.append('hdist2=1')
    adapterlength = min([len(a) for a in adapters])
    call_list.append('k={}'.format(adapterlength-2 if adapterlength>=20 else adapterlength))
    call_list.append('mink={}'.format(max([adapterlength-5,12])))
    
    if not discard_trimmed:
        call_list.append('ktrim=l')
    
    subprocess.call(call_list)


def filter_CDR3(infile,outfile,frame=True,productive=True,start='C',end='FW'):
    def Writer(outhandle,seq_tuple):
        outhandle.write('>%s\n%s\n' % (seq_tuple[0],seq_tuple[1]))
    
    with open(infile) as in_handle, open(outfile,'w') as out_handle:
        for title,seq in reptools.FASTAparser(in_handle):
            if frame:
                if len(seq.strip())%3!=0:
                    continue #next sequence (don't output this one)
            aa = str(reptools.trans(seq.strip()))
            if productive:
                if '*' in aa:
                    continue #next sequence (don't output this one)
            if start:
                correctstart = False
                for codon in start:
                    if aa[0]==codon:
                        correctstart = True
                if not correctstart:
                    continue #next sequence (don't output this one)
            if end:
                correctend = False
                for codon in end:
                    if aa[-1]==codon:
                        correctend = True
                if not correctend:
                    continue #next sequence (don't output this one)
            
            Writer(out_handle,(title,seq)) #sequence is good, write it

