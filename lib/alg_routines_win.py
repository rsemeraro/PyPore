
#! /usr/bin/env python
import argparse
import multiprocessing
import os, warnings
import re
import shutil
import subprocess
import sys
import time
import ctypes
from time import time
from logging_module import log

if __name__ == '__main__':
    if not len(sys.argv) > 7:
        sys.exit(1)     
    work_dir = sys.argv[1]
    fast_Q_file = sys.argv[2]
    ref = sys.argv[3]
    stats_trigger = sys.argv[4]
    prefix = sys.argv[5]
    out_dir = sys.argv[6]
    als = sys.argv[7]
    th = int(sys.argv[8])
    warnings.filterwarnings('ignore', category=UserWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    import plotly
    import plotly.graph_objs as go
    from numpy import sort, digitize, arange, diff, argmax, mean, asarray, int64
    from plotly import tools
    tmpdir = 'alg.tmp'

    def verbosity(ver):
        if ver == 0:
            lib.logging_module.log.setLevel(lib.logging_module.logging.WARN)
        elif ver == 1:
            lib.logging_module.log.setLevel(lib.logging_module.logging.INFO)
        elif ver == 2:
            lib.logging_module.log.setLevel(lib.logging_module.logging.DEBUG)
        elif ver == 3:
            sys.stdout = open(os.devnull, 'w')
            lib.logging_module.log.setLevel(lib.logging_module.logging.DEBUG)
            lib.logging_module.ch.setLevel(lib.logging_module.logging.DEBUG)
    
    def plot_stats(out_dict, s_unmap, s_map, c_c_dict, odir):
        fig = tools.make_subplots(rows=3, cols=2, specs=[[{}, {}], [{'colspan': 2}, None], [{'colspan': 2}, None]],
                                shared_xaxes=False,
                                shared_yaxes=False, vertical_spacing=0.1, print_grid=False)
    
        trace1 = go.Bar(
            x=map(lambda x: x[1], s_map),
            y=map(lambda x: '_' + str(x[0]), s_map),
            name='Mapped',
            orientation='h',
            showlegend=True,
            visible=True,
            text=map(lambda x: str(x[1]), s_map),
            marker=dict(color='rgba(50, 171, 96, 0.6)',
                        line=dict(
                            color='rgba(50, 171, 96, 1.0)',
                            width=0.3)),
            hoverinfo="text+name"
        )
        fig.append_trace(trace1, 1, 1)
        trace2 = go.Bar(
            x=map(lambda x: x[2], s_map),
            y=map(lambda x: '_' + str(x[0]), s_map),
            name='Unmapped',
            orientation='h',
            showlegend=True,
            visible=True,
            text=map(lambda x: str(x[2]), s_map),
            marker=dict(color='rgba(133, 37, 25, 0.6)',
                        line=dict(
                            color='rgba(133, 37, 25, 1.0)',
                            width=0.3)),
            hoverinfo="text+name"
        )
        fig.append_trace(trace2, 1, 1)
    
        trace3 = go.Scatter(
            x=[x if x is not 0    else None for x in s_unmap],
            y=map(lambda x: '_' + str(x[0]), s_map),
            showlegend=False,
            mode='lines+text',
            name='Mapped Fraction',
            fill='toself',
            hoverinfo='text+name',
            text=[x if x is not 0    else None for x in s_unmap],
            line=dict(width=3,
                    color='rgb(188,189,220)'),
            xaxis='x2',
            yaxis='y'
        )
        fig.append_trace(trace3, 1, 2)
    
        cols = ['rgb(166,206,227)', 'rgb(31,120,180)', 'rgb(178,223,138)', 'rgb(51,160,44)', 'rgb(251,154,153)',
                'rgb(227,26,28)', 'rgb(253,191,111)', 'rgb(255,127,0)', 'rgb(202,178,214)', 'rgb(106,61,154)',
                'rgb(204,204,0)']
    
        chromosmes = map(lambda x: 'chr' + str(x), range(1, 23))
        last_pos = 0
        c_counter = 0
        difference = 0
        mean_val_for_axes = 0
        for idx, c in enumerate(chromosmes):
            pos_vec = [list(x) for x in zip(*sorted(c_c_dict[str(c)], key=lambda pair: pair[0]))][0]
            if idx > 0:
                difference = (pos_vec[0] + last_pos) - (last_pos + 1)
            cov_vec = [list(x) for x in zip(*sorted(c_c_dict[str(c)], key=lambda pair: pair[0]))][1]
            if mean([mean_val_for_axes, max(cov_vec)]) > mean_val_for_axes:
                mean_val_for_axes = mean(cov_vec)
                mean_val_for_axes = mean([mean_val_for_axes, max(cov_vec)])
            new_pos_vec = [i + last_pos - difference for i in asarray(pos_vec, dtype=int64)]
            last_pos = new_pos_vec[-1]
            GapIdx = argmax(diff(new_pos_vec))
            rep_idxs = [GapIdx, GapIdx + 1]
            repl = [None, None]
            for repidx, rl in zip(rep_idxs, repl):
                new_pos_vec[repidx] = rl
                cov_vec[repidx] = rl
            if idx == 11:
                c_counter = 0
            trace5 = go.Scatter(
                x=new_pos_vec,
                y=[round(i, 3) for i in cov_vec if i is not None],
                name=str(c),
                fill='tozeroy',
                showlegend=False,
                text=map(lambda x: str(x) + 'x', [round(i, 3) for i in cov_vec if i is not None]),
                hoverinfo="x+text+name",
                line=dict(width=0.5,
                        color=cols[c_counter])
            )
            c_counter += 1
            fig.append_trace(trace5, 3, 1)
    
        vars_cols = ['rgba(255, 144, 14, 1)', 'rgb(255, 65, 54)', 'rgb(93, 164, 214)']
        var_dict = {'I': 'Insertions',
                    'M': 'SNVs',
                    'D': 'Deletions'}
        for n, v in enumerate(['M', 'D', 'I']):
            trace6 = go.Scatter(
                x=sorted(out_dict[v].keys()),
                y=map(lambda x: float(out_dict[v].get(x)), sorted(out_dict[v].keys())),
                name=var_dict.get(v),
                fill='tozeroy',
                line=dict(width=1,
                        color=vars_cols[n],
                        ),
                text=[str(i) + 'bp' for i in sorted(out_dict[v].keys()) if i is not None],
                hoverinfo="text+y+name",
            )
            fig.append_trace(trace6, 2, 1)
    
        fig['layout'].update(
            height=1920,
            width=1920,
            title='Alignment Summary',
            spikedistance=-1,
            yaxis1=dict(
                title='Sequence size (bp)',
                showgrid=False,
                showline=True,
                fixedrange=True,
                hoverformat='.2s',
                showspikes=True,
                spikethickness=1.5,
                spikemode='across',
                spikecolor='grey',
                spikedash='solid',
                tickvals=map(lambda x: '_' + str(x[0]), s_map),
                ticktext=list(zip(*s_map)[0]),
            ),
            yaxis2=dict(
                showgrid=False,
                showline=False,
                linecolor='rgba(102, 102, 102, 0.8)',
                fixedrange=True,
                hoverformat='.2s'
            ),
            xaxis1=dict(
                title='Sequence Number',
                zeroline=False,
                showline=False,
                showgrid=True,
                fixedrange=True,
            ),
            xaxis2=dict(
                title='Mapped Sequences Fraction',
                zeroline=False,
                showline=False,
                showgrid=True,
                fixedrange=True,
                side='top',
            ),
            xaxis4=dict(showspikes=False, title='Position (bp)', range=[0, last_pos + 500], anchor='y4',
                        rangemode='nonnegative', constrain='range', exponentformat='SI'),
            yaxis4=dict(showspikes=False, title='Coverage (X)', range=[0, round(mean_val_for_axes, 2)],
                        rangemode='nonnegative', constrain='range'),
            xaxis3=dict(title='Position (bp)', anchor='y3', exponentformat='SI', domain=[0, 1], type='log',
                        hoverformat='.2s', showspikes=True, spikethickness=1.5, spikemode='across', spikecolor='grey',
                        spikedash='solid'),
            yaxis3=dict(title='Error rate', showline=True, showspikes=False, hoverformat='.0%', tickformat='.0%'),
            legend=dict(
                x=-0.009000000000000001,
                y=1.0222754491017962,
                orientation='h',
                traceorder='normal',
                font=dict(
                    size=10,
                ),
            ),
            margin=dict(
                l=150,
                r=100,
                t=120,
                b=150,
            ),
            barmode='stack',
            hovermode='closest'
        )
        config = {'displaylogo': False,
                'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud', 'select2d', 'lasso2d',
                                             'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleHover', 'zoomIn2d',
                                             'zoomOut2d', 'toggleSpikelines']}
        for i in fig.data:
            for k, v in i.iteritems():
                phrase = str(k) + ':' + str(v)
                if phrase == 'yaxis:y2':
                    i[k] = 'y'
        out_stats = os.path.join(odir, (prefix + '_alignment_stats.html'))
        plotly.offline.plot(fig, filename=out_stats, auto_open=False, config=config, show_link=False)
    
    
    def sam_parser(bwastdout, out_dir):
        log.info('Error estimation...')
        bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900,
                1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 35000, 40000,
                45000, 50000, 100000, 200000, 300000, 400000, 500000, 750000, 1000000000]
        mapped_frac_size = {k: [0, 0] for k in bins}
        pos_dict = {k: {} for k in arange(0, 250000000, 1000000)}
        OutDict = {k: {} for k in ['M', 'I', 'D']}
        header_list = []
        proc_lists = {k: [] for k in range(th)}
        chr_cov_dict = {}
        file_list = []
        fill_val = 0
        Flag = True
        while Flag:
            line = bwastdout.stdout.readline()
            Threshold = 50000
            while len(proc_lists[th - 1]) <= Threshold:
                proc_check_first = map(lambda x: len(proc_lists[x]), range(th - 1))
                line_counter = 0
                if line.strip() == '':
                    Flag = False
                    break
                if line[0] == '@':
                    header_list.append(line)
                    line = bwastdout.stdout.readline()
                else:
                    while line_counter < len(proc_lists):
                        proc_lists[line_counter].append(line)
                        line_counter += 1
                        line = bwastdout.stdout.readline()
                    line_counter = 0
                    line = bwastdout.stdout.readline()
                proc_check_second = map(lambda x: len(proc_lists[x]), range(th - 1))
                if all(v == 0 for v in proc_check_second) == False:
                    if proc_check_second == proc_check_first:
                        time.sleep(5)
                        proc_check_second = map(lambda x: len(proc_lists[x]), range(th - 1))
                        if proc_check_second == proc_check_first:
                            break
            fill_list = (fill_val, fill_val + (th - 1) - map(lambda x: len(proc_lists[x]), range(th - 1)).count(0))
            fill_val = fill_val + (th - 1) - map(lambda x: len(proc_lists[x]), range(th - 1)).count(0)
            res_obj = error_wrap(proc_lists, header_list, fill_list)
            for ro in res_obj:
                p_d, OD, m_f_s, fname = ro
                file_list.append(fname)
                for k, v in p_d.iteritems():
                    if v is not {}:
                        for ch, l in v.iteritems():
                            pos_dict[k][ch] = pos_dict[k].get(ch, []) + l
                for var_k, cnt in OD.iteritems():
                    if cnt is not {}:
                        for k, v in cnt.iteritems():
                            outdictget = OutDict.get(var_k)
                            if outdictget.get(int(k)) == None:
                                outdictget[int(k)] = v
                            else:
                                outdictget[int(k)] = [outdictget.get(int(k))[0], int(outdictget.get(int(k))[1]) + int(v[1]),
                                                    int(outdictget.get(int(k))[2]) + int(v[2])]
                for mfs in m_f_s:
                    mapped_frac_size[mfs[0]] = map(sum, zip(mapped_frac_size[mfs[0]], mfs[1:]))
                proc_lists = {k: [] for k in range(th)}
        for ev in ['I', 'M', 'D']:
            OutDict[ev] = dict(zip(sorted(OutDict[ev].keys()),
                                 zip([sorted(OutDict[ev].keys())[0]] + list(diff(sorted(OutDict[ev].keys()))), [
                                     max(zip(*map(lambda x: OutDict[ev].get(x), sorted(OutDict[ev].keys())))[1]) if x[0] <
                                                                                                                    zip(
                                                                                                                        *map(
                                                                                                                            lambda
                                                                                                                                x:
                                                                                                                            OutDict[
                                                                                                                                ev].get(
                                                                                                                                x),
                                                                                                                            sorted(
                                                                                                                                OutDict[
                                                                                                                                    ev].keys())))[
                                                                                                                        1].index(
                                                                                                                        max(
                                                                                                                            zip(
                                                                                                                                *map(
                                                                                                                                    lambda
                                                                                                                                        x:
                                                                                                                                    OutDict[
                                                                                                                                        ev].get(
                                                                                                                                        x),
                                                                                                                                    sorted(
                                                                                                                                        OutDict[
                                                                                                                                            ev].keys())))[
                                                                                                                                1])) else
                                     x[1] for x in
                                     enumerate(zip(*map(lambda x: OutDict[ev].get(x), sorted(OutDict[ev].keys())))[1])],
                                     zip(*map(lambda x: OutDict[ev].get(x), sorted(OutDict[ev].keys())))[2])))
        for k, v in pos_dict.iteritems():
            for k1, v1 in v.iteritems():
                chr_cov_dict[str(k1)] = chr_cov_dict.get(str(k1), []) + [(k, sum(v1) / 1000000.)]
        mapped_frac_size = [[k[0], k[1][0], k[1][1], round(1 - float(k[1][1]) / (float(k[1][0]) + float(k[1][1])), 3)] for k
                            in sorted(mapped_frac_size.iteritems()) if k[1][0] != 0 or k[1][1] != 0]
        sorted_unmapfraqseq = zip(*mapped_frac_size)[3]
        for ev in OutDict.keys():
            for k, v in sorted(OutDict[ev].iteritems()):
                OutDict[ev][k] = round(float(v[2]) / (int(v[0]) * float(v[1])), 4)
        plot_stats(OutDict, sorted_unmapfraqseq, mapped_frac_size, chr_cov_dict, out_dir)
        finalfile = os.path.join(tmpdir, (prefix + '.bam'))
        bamsfile = os.path.join(out_dir, 'to_merge.txt')
        file = open(bamsfile, 'w')
        for line in file_list:
            file.write(os.path.join(work_dir, line) + '\n')
        file.close()
        bammergeline = subprocess.Popen(['lib\SAMtools\samtools.exe', 'merge', '-cp', '-@', str(th), '-b', bamsfile, finalfile], stdout=subprocess.PIPE).communicate()
        for b in file_list:
            os.remove(b)
        os.remove(bamsfile)
        log.debug('[Alignment][minimap2al] - samtools index %s' % (finalfile))
        samindexline = subprocess.Popen(['lib\SAMtools\samtools.exe', 'index', finalfile], stdout=subprocess.PIPE).communicate()
        for f in os.listdir(tmpdir):
        try:
            shutil.move(os.path.join(tmpdir,f), os.path.join(work_dir, out_dir))
        except shutil.Error:
            log.error("Unable to move %s" % tmpdir,f)            
        shutil.rmtree(tmpdir, ignore_errors=True)
    
    
    def error_wrap(pl, hl, fl):
        from MPC import Consumer, error_calc
        os.mkdir(tmpdir)
        ctypes.windll.kernel32.SetFileAttributesW(tmpdir, 2)
        Q = multiprocessing.JoinableQueue()
        R = multiprocessing.Queue()
        consumers = [Consumer(Q, R)
                     for ix in xrange(th)]
        for w in consumers:
            w.start()
        joblists = range(fl[0], fl[1] + 1)
        for p in range(len(joblists)):
            if len(pl[p]) > 0:
                Q.put(error_calc(pl[p], hl, joblists[p], out_dir, prefix))
        for ix in xrange(th):
            Q.put(None)
        Q.join()
        while not Q.empty():
            Q.get()
        for r in range(R.qsize()):
            res = R.get()
            if res is not None:
                p_d, OD, m_f_s, fn = res
                out_data = [p_d, OD, m_f_s, fn]
                yield out_data
    
    
    def minimap2al():
        mm_ext_dir = os.path.join(work_dir, out_dir, 'minimap2_alignments')
        if not os.path.exists(mm_ext_dir):
            os.makedirs(mm_ext_dir)
        mmi_ref = os.path.splitext(ref)[0] + '.mmi'
        if not os.path.exists(mmi_ref):
            log.debug('[Alignment][minimap2al] - minimap2 -x map-ont -t %s -d %s %s' % (th, mmi_ref, ref))
            miniline = ('lib\minimap2\minimap2.exe -x map-ont -t', str(th), '-d', str(mmi_ref), str(ref))
            minrun = ' '.join(miniline)
            subprocess.Popen(minrun, stdout=subprocess.PIPE, shell=True).wait()
        bam_mm2_file = os.path.join(mm_ext_dir, (prefix + '.bam'))
        if not os.path.exists(bam_mm2_file):
            if stats_trigger in ['y', 'yes']:
                log.debug('[Alignment][minimap2al] - minimap2 -ax map-ont --MD -L -t %s -R %s %s %s' % (
                th, str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)), str(mmi_ref),
                str(fast_Q_file)))
                minimap2line = subprocess.Popen(['lib\minimap2\minimap2.exe', '-ax', 'map-ont', '--MD', '-L', '-t', str(th), '-R', str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)), str(mmi_ref), str(fast_Q_file)], stdout=subprocess.PIPE)
                PairDict = sam_parser(minimap2line, mm_ext_dir)
            else:
                sam_mm2_file = os.path.join(mm_ext_dir, (prefix + '.sam'))
                log.debug('[Alignment][minimap2al] - minimap2 -ax map-ont --MD -L -t %s -R %s %s %s > %s' % (
                th, str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)), str(mmi_ref),
                str(fast_Q_file), str(sam_mm2_file)))
                with open(sam_mm2_file,"w") as out:
                minimap2line = subprocess.Popen(['lib\minimap2\minimap2.exe', '-ax', 'map-ont', '--MD', '-L', '-t', str(th), '-R',
                                                 str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)),
                                                 str(mmi_ref), str(fast_Q_file)], stdout=out).communicate()                            
                outputfilebam = os.path.join(mm_ext_dir, (prefix + '.tmp.bam'))
                log.debug('[Alignment][minimap2al] - samtools view -Sb -@ %s %s -o %s' % (th, sam_mm2_file, outputfilebam))
                samviewline = subprocess.Popen(['lib\SAMtools\samtools.exe', 'view', '-Sb', '-@', str(th), sam_mm2_file, '-o', outputfilebam], stdout=subprocess.PIPE).communicate()
                os.remove(sam_mm2_file)
                samsortline = subprocess.Popen(['lib\SAMtools\samtools.exe', 'sort', '-@', str(th), outputfilebam, '-o', bam_mm2_file], stdout=subprocess.PIPE).communicate()
                log.debug('[Alignment][minimap2al] - samtools index %s -@%s' % (bam_mm2_file, str(th)))
                samindexline = subprocess.Popen(['lib\SAMtools\samtools.exe', 'index', bam_mm2_file], stdout=subprocess.PIPE).communicate()
                os.remove(outputfilebam)
        else:
            log.warning('[Alignment][minimap2al] - file %s already exists!' % bam_mm2_file)
        #try:
        #        shutil.move(mm_ext_dir, os.path.join(work_dir, out_dir))
        #except shutil.Error:
        #        log.error("Unable to move %s" % mm_ext_dir)
    
    
    def bwaal():
        log.warning('[Alignment][bwaal] - Currently not available for Windows, sorry.')
        sys.exit(1)
        bwa_ext_dir = os.path.join(work_dir, out_dir,'bwa_alignments')
        if not os.path.exists(bwa_ext_dir):
            os.makedirs(bwa_ext_dir)
        bam_bwa_file = os.path.join(bwa_ext_dir, (prefix + '.bam'))
        if not os.path.exists(bam_bwa_file):
            if stats_trigger in ['y', 'yes']:
                log.debug('[Alignment][bwaal] - bwa mem -x ont2d -t %s %s %s' % (th, ref, fast_Q_file))
                bwaline = subprocess.Popen(['lib\BWA\bwa.exe', 'mem', '-x', 'ont2d', '-t', str(th), str(ref), str(fast_Q_file)],
                                         stdout=subprocess.PIPE)
                PairDict = sam_parser(bwaline, bwa_ext_dir)
            else:
                sam_bwa_file = os.path.join(bwa_ext_dir, (prefix + '.sam'))
                log.debug('[Alignment][bwaal] - bwa mem -x ont2d -t %s %s %s > %s' % (th, ref, fast_Q_file, sam_bwa_file))
                bwaline = subprocess.Popen(
                    ['bwa', 'mem', '-x', 'ont2d', '-t', str(th), str(ref), str(fast_Q_file), '>', str(sam_bwa_file)],
                    stdout=subprocess.PIPE)
                outputfilebam = os.path.join(bwa_ext_dir, (prefix + '.tmp.bam'))
                log.debug('[Alignment][bwaal] - samtools view -Sb -@ %s %s -o %s' % (th, sam_bwa_file, outputfilebam))
                pysam.view("-Sb", "-@%s" % str(th), sam_bwa_file, "-o%s" % outputfilebam, catch_stdout=False)
                os.remove(sam_bwa_file)
                pysam.sort(outputfilebam, "-o%s" % bam_bwa_file, catch_stdout=False)
                log.debug('[Alignment][bwaal] - samtools index %s -@%s' % (bam_bwa_file, str(th)))
                pysam.index(bam_bwa_file, "-@%s" % str(th), catch_stdout=False)
        else:
            log.warning('[Alignment][bwaal] - file %s already exists!' % bam_bwa_file)
        try:
            shutil.move(bwa_ext_dir, os.path.join(work_dir, out_dir))
        except shutil.Error:
            log.error("Unable to move %s" % bwa_ext_dir)
    
     
    def ngmlral():
        log.warning('[Alignment][ngmlral] - Currently not available for Windows, sorry.')
        sys.exit(1)        
        ngmlr_ext_dir = os.path.join(work_dir, out_dir, 'ngmlr_alignments')
        if not os.path.exists(ngmlr_ext_dir):
            os.makedirs(ngmlr_ext_dir)
        bam_ngmlr_file = os.path.join(ngmlr_ext_dir, (prefix + '.bam'))
        if not os.path.exists(bam_ngmlr_file):
            if stats_trigger in ['y', 'yes']:
                log.debug('[Alignment][ngmlral] - ngmlr -t %s -r %s -q %s -x ont' % (th, ref, fast_Q_file))
                ngmlrline = subprocess.Popen(['ngmlr', '-t', str(th), '-r', str(ref), '-q', str(fast_Q_file), '-x ont'],
                                             stdout=subprocess.PIPE)
                PairDict = sam_parser(ngmlrline, ngmlr_ext_dir)
            else:
                sam_nglmr_file = os.path.join(ngmlr_ext_dir, (prefix + '.sam'))
                log.debug(
                    '[Alignment][ngmlral] - ngmlr -t %s -r %s -q %s -o %s -x ont' % (th, ref, fast_Q_file, sam_ngmlr_file))
                ngmlrline = subprocess.Popen(
                    ['ngmlr', '-t', str(th), '-r', str(ref), '-q', str(fast_Q_file), '-o', str(sam_ngmlr_file), '-x ont'],
                    stdout=subprocess.PIPE)
                outputfilebam = os.path.join(ngmlr_ext_dir, (prefix + '.tmp.bam'))
                log.debug('[Alignment][ngmlral] - samtools view -Sb -@ %s %s -o %s' % (th, sam_nglmr_file, outputfilebam))
                pysam.view("-Sb", "-@%s" % str(th), sam_nglmr_file, "-o%s" % outputfilebam, catch_stdout=False)
                os.remove(sam_nglmr_file)
                pysam.sort(outputfilebam, "-o%s" % bam_ngmlr_file, catch_stdout=False)
                log.debug('[Alignment][ngmlral] - samtools index %s -@%s' % (bam_ngmlr_file, str(th)))
                pysam.index(bam_ngmlr_file, "-@%s" % str(th), catch_stdout=False)
        else:
            log.warning('[Alignment][ngmlral] - file %s already exists!' % bam_ngmlr_file)
        try:
            shutil.move(ngmlr_ext_dir, os.path.join(work_dir, out_dir))
        except shutil.Error:
            log.error("Unable to move %s" % ngmlr_ext_dir)
    
    
    def als_parser(als):
        if not os.path.exists(os.path.join(work_dir, out_dir)):
            os.makedirs(os.path.join(work_dir, out_dir))
        tooldict = {'b': bwaal,
                    'm': minimap2al,
                    'n': ngmlral}
    
        if isinstance(als, list):
            for a in als:
                log.info('[Alignment] - Start alignment with %s...' % str(tooldict.get(a)).split(' ')[1][:-2])
                tooldict[a]()
                log.info('[Alignment] - Finish alignment with %s...' % str(tooldict.get(a)).split(' ')[1][:-2])
        else:
            log.info('[Alignment] - Start alignment with %s...' % str(tooldict.get(als)).split(' ')[1][:-2])
            tooldict[als]()
            log.info('[Alignment] - Finish alignment with %s...' % str(tooldict.get(als)).split(' ')[1][:-2])
    
    als_parser(als)