#! /usr/bin/env python
import argparse
import multiprocessing
import multiprocessing.queues
import os, re, shutil, sys, time, warnings
import subprocess
from time import time
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import plotly
import plotly.graph_objs as go
from plotly import tools
from logging_module import log
from numpy import sort, digitize, arange, diff, argmax, mean
from numpy import diff as diff__
#import lib
if os.name != 'nt':
    import pysam
from IPython import embed

class SharedCounter(object):
    """ A synchronized shared counter.
    The locking done by multiprocessing.Value ensures that only a single
    process or thread may read or write the in-memory ctypes object. However,
    in order to do n += 1, Python performs a read followed by a write, so a
    second process may read the old value before the new one is written by the
    first process. The solution is to use a multiprocessing.Lock to guarantee
    the atomicity of the modifications to Value.
    This class comes almost entirely from Eli Bendersky's blog:
    http://eli.thegreenplace.net/2012/01/04/shared-counter-with-pythons-multiprocessing/
    """

    def __init__(self, n = 0):
        self.count = multiprocessing.Value('i', n)

    def increment(self, n = 1):
        """ Increment the counter by n (default = 1) """
        with self.count.get_lock():
            self.count.value += n

    @property
    def value(self):
        """ Return the value of the counter """
        return self.count.value


class Queue(multiprocessing.queues.Queue):
    """ A portable implementation of multiprocessing.Queue.
    Because of multithreading / multiprocessing semantics, Queue.qsize() may
    raise the NotImplementedError exception on Unix platforms like Mac OS X
    where sem_getvalue() is not implemented. This subclass addresses this
    problem by using a synchronized shared counter (initialized to zero) and
    increasing / decreasing its value every time the put() and get() methods
    are called, respectively. This not only prevents NotImplementedError from
    being raised, but also allows us to implement a reliable version of both
    qsize() and empty().
    """

    def __init__(self, *args, **kwargs):
        super(Queue, self).__init__(*args, **kwargs)
        self.size = SharedCounter(0)

    def put(self, *args, **kwargs):
        self.size.increment(1)
        super(Queue, self).put(*args, **kwargs)

    def get(self, *args, **kwargs):
        self.size.increment(-1)
        return super(Queue, self).get(*args, **kwargs)

    def qsize(self):
        """ Reliable implementation of multiprocessing.Queue.qsize() """
        return self.size.value

    def empty(self):
        """ Reliable implementation of multiprocessing.Queue.empty() """
        return not self.qsize()

    def clear(self):
        """ Remove all elements from the Queue. """
        while not self.empty():
            self.get()

class Consumer(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, main='.'):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.main = main

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                self.task_queue.task_done()
                break
            next_task.main = self.main
            next_task.sub = proc_name
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return
      
class error_calc_w32(object):
    def __init__(self, sampled_data, header_lines, file_number, mpc_out_dir, samp_pref):
        self.sampled_data = sampled_data
        self.header_lines = header_lines
        self.file_number = file_number
        self.mpc_out_dir = mpc_out_dir
        self.samp_pref = samp_pref
        self.tmpdir = 'alg.tmp'

    def __call__(self):
        variants_dict = {k: [] for k in ['M', 'I', 'D']}
        RCounter = []
        SeenNames = set()
        pos_dict = {k: {} for k in arange(0, 250000000, 1000000)}
        unmapped = []
        file_out = os.path.join(str(self.tmpdir), str(self.samp_pref) + '.tmp.' + str(self.file_number) + '.sam')
        file = open(file_out, mode='wb')
        for hl in self.header_lines:
            if hl != '':
                file.write(hl)
        for ll in self.sampled_data:
            if ll.strip() is not '':
                file.write(ll)
                pair = ll.strip().split()
                flag = pair[1]
                if int(flag) not in [272, 256]:
                    a_cigar = pair[5]
                    map_quality = int(pair[4])
                    chrome = str(pair[2])
                    read_name = str(pair[0])
                    read_len = len(pair[9])
                    most_left = int(pair[3])
                    if read_name not in SeenNames:
                        SeenNames.add(read_name)
                        if chrome in map(str, range(1, 23)) + map(lambda x: 'chr' + str(x), range(1, 23)):
                            kkk = pos_dict.get(most_left,
                                               pos_dict[min(pos_dict.keys(), key=lambda k: abs(k - most_left))])
                            kkk[chrome] = kkk.get(chrome, []) + [read_len]
                        SClip = 0
                        if a_cigar != '*':
                            try:
                                SClip = int(re.findall("[0-9]*S", a_cigar)[0][:-1])
                            except:
                                log.debug('[Alignment][error_calc] - No soft clipping for this read...')
                                continue
                            a_tags = dict([(x[:2], x.split(':')[-1]) for x in pair[11:]])
                            mdgsub = re.sub("([\\^]*[ACGT]+)[0]*", " \\1 ", a_tags['MD']).split()
                            MissMatchesDelsPos = miss_match_founder(mdgsub, SClip)
                            Insertions = parseCIGAR(a_cigar, SClip)
                            for p, e, s in MissMatchesDelsPos:
                                variants_dict[str(e)].append((p, s))
                            variants_dict['I'].extend(list(Insertions))
                            RCounter.append([read_name, read_len, map_quality])
                        else:
                            unmapped.append([read_name, read_len, map_quality])
        file.close()
        outputfilebam = os.path.join(str(self.tmpdir), str(self.samp_pref) + '.tmp.' + str(self.file_number) + '.bam')
        log.debug('[Alignment][error_calc] - samtools view -Sb %s -o %s' % (file_out, outputfilebam))
        samviewline = subprocess.Popen([os.path.abspath(os.path.dirname(__file__)) + '\SAMtools\samtools.exe', 'view', '-Sb', file_out, '-o', outputfilebam], stdout=subprocess.PIPE).communicate()
        os.remove(file_out)
        outputfilebamsorted = os.path.join(str(self.tmpdir), str(self.samp_pref) + '.' + str(self.file_number) + '.bam')
        log.debug('[Alignment][error_calc] - samtools sort %s -o %s' % (outputfilebam, outputfilebamsorted))
        samsortline = subprocess.Popen([os.path.abspath(os.path.dirname(__file__)) + '\SAMtools\samtools.exe', 'sort', outputfilebam, '-o', outputfilebamsorted], stdout=subprocess.PIPE).communicate()
        os.remove(outputfilebam)
        bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800,
                900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 35000,
                40000, 45000, 50000, 100000, 200000, 300000, 400000, 500000, 750000, 1000000000]
        d_lengths = {}
        d_unm = {}
        r_lengths = map(lambda x: int(x[1]), RCounter)
        if len(unmapped) > 0:
            unm_lengths = map(int, zip(*unmapped)[1])
            binned_unm_data = digitize(unm_lengths, bins)
            for l in binned_unm_data:
                d_unm[bins[l]] = d_unm.get(bins[l], 0) + 1
        else:
            d_unm = {k: 0 for k in bins}
        binned_lengths_data = digitize(r_lengths, bins)
        mapped_frac_size = []
        for l in binned_lengths_data:
            d_lengths[bins[l]] = d_lengths.get(bins[l],
                                               0) + 1  ## Conto le read di lunghezza bins[l] e le storo in un dizionario con la lunghezza come chiave
        for le in sorted(set(d_lengths.keys() + d_unm.keys())):
            mapped_frac_size.append(
                [le, d_lengths.get(le, 0), d_unm.get(le, 0)])  # creo righe con lunghezza e conte per mapped e unmapped
        d_lengths = dict(zip(sorted(d_lengths.keys()),
                             [sum(zip(*sorted(d_lengths.iteritems(), reverse=True))[1][:i]) for i in
                              range(1, len(zip(*sorted(d_lengths.iteritems(), reverse=True))[1]) + 1)][::-1]))
        OutDict = {k: {} for k in ['M', 'I', 'D']}
        from numpy import diff as diff__
        for ev in ['I', 'M', 'D']:
            d_somme = {}
            binned_data = digitize(map(int, (zip(*variants_dict[str(ev)])[0])), bins)
            for index, vi in enumerate(binned_data):
                d_somme[bins[vi]] = d_somme.get(bins[vi], 0) + int(variants_dict[str(ev)][index][1])
            bin_sizes = [sort(d_somme.keys())[0]] + list(diff__(sort(d_somme.keys())))
            for i, kv in enumerate(sorted(d_somme.iteritems())):
                kkk = OutDict.get(ev)
                kkk[str(kv[0])] = kkk.get(str(kv[0]), []) + [int(bin_sizes[i]), int(
                    d_lengths.get(kv[0], d_lengths[min(d_lengths.keys(), key=lambda k: abs(k - int(kv[0])))])),
                                                             int(kv[1])]
        out_data = [pos_dict, OutDict, mapped_frac_size, outputfilebamsorted]
        return (out_data)

class error_calc(object):
    def __init__(self, sampled_data, header_lines, file_number):
        self.sampled_data = sampled_data
        self.header_lines = header_lines
        self.file_number = file_number

    def __call__(self):
        variants_dict = {k: [] for k in ['M', 'I', 'D']}
        RCounter = []
        SeenNames = set()
        pos_dict = {k: {} for k in arange(0, 250000000, 1000000)}
        unmapped = []
        file_out = 'tmp.' + str(self.file_number) + '.sam'
        file = open(file_out, mode='wb')
        for hl in self.header_lines:
            if hl != '':
                file.write(hl)
        for ll in self.sampled_data:
            if ll.strip() is not '':
                file.write(ll)
                pair = ll.strip().split()
                flag = pair[1]
                if int(flag) not in [272, 256]:
                    a_cigar = pair[5]
                    map_quality = int(pair[4])
                    chrome = str(pair[2])
                    read_name = str(pair[0])
                    read_len = len(pair[9])
                    most_left = int(pair[3])
                    if read_name not in SeenNames:
                        SeenNames.add(read_name)
                        if chrome in map(str, range(1, 23)) + map(lambda x: 'chr' + str(x), range(1, 23)):
                            kkk = pos_dict.get(most_left,
                                               pos_dict[min(pos_dict.keys(), key=lambda k: abs(k - most_left))])
                            kkk[chrome] = kkk.get(chrome, []) + [read_len]
                        SClip = 0
                        if a_cigar != '*':
                            try:
                                SClip = int(re.findall("[0-9]*S", a_cigar)[0][:-1])
                            except:
                                log.debug('[Alignment][error_calc] - No soft clipping for this read...')
                                continue
                            a_tags = dict([(x[:2], x.split(':')[-1]) for x in pair[11:]])
                            mdgsub = re.sub("([\\^]*[ACGT]+)[0]*", " \\1 ", a_tags['MD']).split()
                            MissMatchesDelsPos = miss_match_founder(mdgsub, SClip)
                            Insertions = parseCIGAR(a_cigar, SClip)
                            for p, e, s in MissMatchesDelsPos:
                                variants_dict[str(e)].append((p, s))
                            variants_dict['I'].extend(list(Insertions))
                            RCounter.append([read_name, read_len, map_quality])
                        else:
                            unmapped.append([read_name, read_len, map_quality])
        file.close()
        outputfilebam = 'tmp.' + str(self.file_number) + '.bam'
        log.debug('[Alignment][error_calc] - samtools view -Sb %s -o %s' % (file_out, outputfilebam))
        pysam.view("-Sb", file_out, "-o%s" % outputfilebam, catch_stdout=False)
        os.remove(file_out)
        outputfilebamsorted = os.path.join(out_dir, str(self.file_number) + '.bam')
        log.debug('[Alignment][error_calc] - samtools sort %s -o %s' % (outputfilebam, outputfilebamsorted))
        pysam.sort(outputfilebam, "-o%s" % outputfilebamsorted, catch_stdout=False)
        os.remove(outputfilebam)
        bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800,
                900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 35000,
                40000, 45000, 50000, 100000, 200000, 300000, 400000, 500000, 750000, 1000000000]
        d_lengths = {}
        d_unm = {}
        r_lengths = map(lambda x: int(x[1]), RCounter)
        if len(unmapped) > 0:
            unm_lengths = map(int, zip(*unmapped)[1])
            binned_unm_data = digitize(unm_lengths, bins)
            for l in binned_unm_data:
                d_unm[bins[l]] = d_unm.get(bins[l], 0) + 1
        else:
            d_unm = {k: 0 for k in bins}
        binned_lengths_data = digitize(r_lengths, bins)
        mapped_frac_size = []
        for l in binned_lengths_data:
            d_lengths[bins[l]] = d_lengths.get(bins[l],
                                               0) + 1  ## Conto le read di lunghezza bins[l] e le storo in un dizionario con la lunghezza come chiave
        for le in sorted(set(d_lengths.keys() + d_unm.keys())):
            mapped_frac_size.append(
                [le, d_lengths.get(le, 0), d_unm.get(le, 0)])  # creo righe con lunghezza e conte per mapped e unmapped
        d_lengths = dict(zip(sorted(d_lengths.keys()),
                             [sum(zip(*sorted(d_lengths.iteritems(), reverse=True))[1][:i]) for i in
                              range(1, len(zip(*sorted(d_lengths.iteritems(), reverse=True))[1]) + 1)][::-1]))
        OutDict = {k: {} for k in ['M', 'I', 'D']}
        for ev in ['I', 'M', 'D']:
            d_somme = {}
            binned_data = digitize(map(int, (zip(*variants_dict[str(ev)])[0])), bins)
            for index, vi in enumerate(binned_data):
                d_somme[bins[vi]] = d_somme.get(bins[vi], 0) + int(variants_dict[str(ev)][index][1])
            bin_sizes = [sort(d_somme.keys())[0]] + list(diff__(sort(d_somme.keys())))
            for i, kv in enumerate(sorted(d_somme.iteritems())):
                kkk = OutDict.get(ev)
                kkk[str(kv[0])] = kkk.get(str(kv[0]), []) + [int(bin_sizes[i]), int(
                    d_lengths.get(kv[0], d_lengths[min(d_lengths.keys(), key=lambda k: abs(k - int(kv[0])))])),
                                                             int(kv[1])]
        out_data = [pos_dict, OutDict, mapped_frac_size, outputfilebamsorted]
        return (out_data)


def parseCIGAR(s, clipping):
    match = re.findall(r'(\d+)(\w)', s)
    MovingIndex = 0
    for n, e in match:
        if e == 'I':
            yield ((MovingIndex + clipping), n)
            pass
        elif e == 'D':
            MovingIndex += int(n)
        elif e in ['M']:
            MovingIndex += int(n)


def miss_match_founder(readMd, clipping):
    MovingIndex = 0
    for i in readMd:
        if i.isdigit():
            MovingIndex += int(i)
        else:
            if i.startswith('^'):
                yield ((MovingIndex + clipping), 'D', len(i) - 1)
                MovingIndex += len(i) - 1
            else:
                MovingIndex += 1
                yield ((MovingIndex + clipping), 'M', 1)


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
        x=[x if x is not 0  else None for x in s_unmap],
        y=map(lambda x: '_' + str(x[0]), s_map),
        showlegend=False,
        mode='lines+text',
        name='Mapped Fraction',
        fill='toself',
        hoverinfo='text+name',
        text=[x if x is not 0  else None for x in s_unmap],
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
    if 'chr' not in c_c_dict.keys()[0]:
        chromosmes = map(str, range(1, 23))    
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
        new_pos_vec = [i + last_pos - difference for i in pos_vec]
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
    finalfile = os.path.join(out_dir, (prefix + '.bam'))
    bamsfile = os.path.join(out_dir, 'to_merge.txt')
    file = open(bamsfile, 'w')
    for line in file_list:
        file.write(os.path.join(work_dir, line) + '\n')
    file.close()
    pysam.merge("-cp", "-@%s" % str(th), "-b%s" % bamsfile, finalfile, catch_stdout=False)
    for b in file_list:
        os.remove(b)
    os.remove(bamsfile)


def error_wrap(pl, hl, fl):
    Q = multiprocessing.JoinableQueue()
    R = Queue()
    consumers = [Consumer(Q, R)
                 for ix in xrange(th)]
    for w in consumers:
        w.start()
    joblists = range(fl[0], fl[1] + 1)
    for p in range(len(joblists)):
        if len(pl[p]) > 0:
            Q.put(error_calc(pl[p], hl, joblists[p]))
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
    mm_ext_dir = os.path.join(work_dir, 'minimap2_alignments')
    if not os.path.exists(mm_ext_dir):
        os.makedirs(mm_ext_dir)
    mmi_ref = os.path.splitext(ref)[0] + '.mmi'
    if not os.path.exists(mmi_ref):
        log.debug('[Alignment][minimap2al] - minimap2 -x map-ont -t %s -d %s %s' % (th, mmi_ref, ref))
        miniline = ('minimap2 -x map-ont -t', str(th), '-d', str(mmi_ref), str(ref))
        minrun = ' '.join(miniline)
        subprocess.Popen(minrun, stdout=subprocess.PIPE, shell=True).wait()
    bam_mm2_file = os.path.join(mm_ext_dir, (prefix + '.bam'))
    if not os.path.exists(bam_mm2_file):
        if stats_trigger in ['y', 'yes']:
            log.debug('[Alignment][minimap2al] - minimap2 -ax map-ont --MD -L -t %s -R %s %s %s' % (
            th, str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)), str(mmi_ref),
            str(fast_Q_file)))
            minimap2line = subprocess.Popen(['minimap2', '-ax', 'map-ont', '--MD', '-L', '-t', str(th), '-R',
                                             str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)),
                                             str(mmi_ref), str(fast_Q_file)], stdout=subprocess.PIPE)
            PairDict = sam_parser(minimap2line, mm_ext_dir)
        else:
            sam_mm2_file = os.path.join(mm_ext_dir, (prefix + '.sam'))
            log.debug('[Alignment][minimap2al] - minimap2 -ax map-ont --MD -L -t %s -R %s %s %s > %s' % (
            th, str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)), str(mmi_ref),
            str(fast_Q_file), str(sam_mm2_file)))
            minimap2line = subprocess.Popen(['minimap2', '-ax', 'map-ont', '--MD', '-L', '-t', str(th), '-R',
                                             str('@RG\\tID:minimap2\\tLB:NxEr\\tPL:MinION\\tPU:NA\\tSM:' + str(prefix)),
                                             str(mmi_ref), str(fast_Q_file), '>', str(sam_mm2_file)],
                                            stdout=subprocess.PIPE)
            outputfilebam = os.path.join(mm_ext_dir, (prefix + '.tmp.bam'))
            log.debug('[Alignment][minimap2al] - samtools view -Sb -@ %s %s -o %s' % (th, sam_mm2_file, outputfilebam))
            pysam.view("-Sb", "-@%s" % str(th), sam_mm2_file, "-o%s" % outputfilebam, catch_stdout=False)
            os.remove(sam_mm2_file)
            pysam.sort(outputfilebam, "-o%s" % bam_mm2_file, catch_stdout=False)
            log.debug('[Alignment][minimap2al] - samtools index %s -@%s' % (bam_mm2_file, str(th)))
            pysam.index(bam_mm2_file, "-@%s" % str(th), catch_stdout=False)
    else:
        log.warning('[Alignment][minimap2al] - file %s already exists!' % bam_mm2_file)
    try:
        shutil.move(mm_ext_dir, os.path.join(work_dir, out_dir))
    except shutil.Error:
        log.error("Unable to move %s" % mm_ext_dir)


def bwaal():
    bwa_ext_dir = os.path.join(work_dir, 'bwa_alignments')
    if not os.path.exists(bwa_ext_dir):
        os.makedirs(bwa_ext_dir)
    bam_bwa_file = os.path.join(bwa_ext_dir, (prefix + '.bam'))
    if not os.path.exists(bam_bwa_file):
        if stats_trigger in ['y', 'yes']:
            log.debug('[Alignment][bwaal] - bwa mem -x ont2d -t %s %s %s' % (th, ref, fast_Q_file))
            bwaline = subprocess.Popen(['bwa', 'mem', '-x', 'ont2d', '-t', str(th), str(ref), str(fast_Q_file)],
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
    ngmlr_ext_dir = os.path.join(work_dir, 'ngmlr_alignments')
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


def run(arguments):
    if not len(arguments) > 7:
        sys.exit(1)
    global work_dir, fast_Q_file, out_dir, ref, th, prefix, stats_trigger        
    work_dir = arguments[0]
    fast_Q_file = arguments[1]
    ref = arguments[2]
    stats_trigger = arguments[3]
    prefix = arguments[4]
    out_dir = arguments[5]
    als = arguments[6]
    th = int(arguments[7])
    als_parser(als)
