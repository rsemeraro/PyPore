import multiprocessing, threading
import sys, os, shutil, traceback, re, ntpath, random, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)
import time
from numpy import sort, mean
from Bio.SeqUtils import GC
from Bio import SeqIO
import plotly
import plotly.graph_objs as go
from logging_module import log

####################
## Define classes ##
####################


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


class fastq_reader(object):
    def __init__(self, datas, fqfolder):
        self.datas = datas
        self.fqfolder = fqfolder
        
    def __call__(self):
        datas = self.datas
        GC_list = []
        for fqs in datas:
            with open(os.path.join(self.fqfolder, fqs), "rU") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    idx = int(record.description.split()[3][3:])-1
                    GC_list.append(GC(record.seq))           
        return GC_list


############################
## Define main functions ###
############################ 


def HeatTrigger2(html_file_2):
    log.debug("Triggering stats report...")
    try:
        fh = file(html_file_2, 'r')
    except IOError:
        e = sys.exc_info()
        tb = traceback.format_exception(*e)
        log.error(
            "No such file or directory!\n\033[36m##################################\n\033[93m%s\033[32m%s\033[93m%s\033[36m##################################\n\033[0m" % (
                tb[0], tb[1], tb[2]))
        sys.exit(1)
    subject = fh.read()
    fh.close()
    Idpattern = '(?<=\<div id=")(.*?)"'
    try:
        FileId = re.search(Idpattern, subject).group(1)
    except:
        e = sys.exc_info()
        tb = traceback.format_exception(*e)
        log.error(
            "Can't find ID in html file!\n\033[36m##################################\n\033[93m%s\033[32m%s\033[93m%s\033[36m##################################\n\033[0m" % (
                tb[0], tb[1], tb[2]))
        sys.exit(1)
    FixLine = "<script type = 'text/javascript'>var graph = document.getElementById('%(FileId)s');graph.on('plotly_hover', function(data){var divname = data.points.map(function(d){return (d.data.name)});if(divname == 'ReadsPerHour' || divname == 'BasesPerHour' || divname == 'MeanReadLengthPerHour'){if(data.xvals){Plotly.Fx.hover(graph, {xval:data.xvals[0] }, ['x5y', 'x5y3','x5y5']);}}});</script>" % locals()
    SubLine = FixLine + '</body></html>'
    pattern = re.compile('</body></html>')
    result = pattern.sub(SubLine, subject)
    f_out = file(html_file_2, 'w')
    f_out.write(result)
    f_out.close()
    try:
        shutil.move(html_file_2, os.path.join(work_dir, out_dir))
    except shutil.Error:
        log.error("Unable to move %s" % html_file_2)


def FastStats(ChannelOut, GCTable = None):
    log.debug('Preparing data for summary...')
    HoursDict = {str(il): [0, 0, 0] for il in range(0, 48)}
    SmallDict = dict([(k, v) for k, v in ChannelOut.items() if len(v) > 0])
    FailSum = sum(map(lambda x: int(SmallDict[str(x)][2]), SmallDict.keys()))
    if GCTable == None:
        for c in SmallDict:
            Sd = dict([(k, v) for k, v in ChannelOut[str(c)][4].items() if len(v) > 0])
            for k, v in Sd.items(): #### itera sui pori
                for o, r in v[0].items(): #### itera sulle ore
                    HoursDict[o][0] += r[0]
                    HoursDict[o][1] += r[1]
    else:            
        for c in SmallDict:
            Sd = dict([(k, v) for k, v in ChannelOut[str(c)][4].items() if len(v) > 0])
            for o, r in Sd.items(): #### itera sulle ore
                HoursDict[o][0] += r[0]
                HoursDict[o][1] += r[1]
    PassSum = sum(map(lambda x: int(HoursDict[x][0]), HoursDict.keys()))
    for x in sort([int(k) for k, v in HoursDict.items() if v != [0, 0, 0]]):
        HoursDict[str(x)][2] = str(round(float(HoursDict[str(x)][1]) / float(HoursDict[str(x)][0])))
    MaxMinvec = []
    LengthsVec = [0]
    QualsVec = [0]
    if GCTable == None:
        GCVec = [0]
    else:
        GCVec = GCTable
    for c in ChannelOut:
        try:
            Bait = map(int, zip(*ChannelOut[c][3])[1])
            LengthsVec.extend(Bait)
            QualsVec.extend(map(float, zip(*ChannelOut[c][3])[2]))
            if GCTable == None:
                GCVec.extend(map(float, zip(*ChannelOut[c][3])[4]))
            MaxMinvec.append(max(Bait))
            MaxMinvec.append(min(Bait))            
        except: pass    
    SmallestRead = min(filter(lambda a: a > 0, MaxMinvec))
    GCVec = filter(lambda a: a > 0, GCVec)  
    BiggestRead = max(MaxMinvec)
    MeanReadLength = round(mean(map(float, map(lambda x: HoursDict[str(x)][2],
                                               sort([int(k) for k, v in HoursDict.items() if v != [0, 0]])))), 1)
    TotalThroughput = sum(map(float, map(lambda x: HoursDict[str(x)][1],
                                                   sort([int(k) for k, v in HoursDict.items() if v != [0, 0]]))))                                               
    sortingindex = [y[0] for y in sorted(enumerate(map(int, HoursDict.keys())), key=lambda z: z[1])]
    SortKeys = [HoursDict.keys()[y] for y in sortingindex]
    HoursDict = [HoursDict[y] for y in SortKeys]
    log.debug('Plotting stats summary...')
    fig = plotly.tools.make_subplots(rows=16, cols=2, shared_xaxes=False,
                                     specs=[[{'rowspan': 3}, {}], [None, None], [None, None], [{'rowspan': 3}, {}],
                                            [None, None], [None, None], [{'rowspan': 3}, {}],
                                            [None, None], [None, None], [None, None], [None, None],
                                            [{'colspan': 2}, {}], [None, None], [None, None], [{'colspan': 2}, {}],
                                            [None, None]],
                                     vertical_spacing=0.001, print_grid=False)
    log.debug('Plotting pie chart...')
    PieChartReads = go.Pie(values=[round(((float(PassSum) / (FailSum + PassSum))) * 100, 2),
                                   round(((float(FailSum) / (FailSum + PassSum))) * 100, 2)],
                           labels=[
                               'Pass',
                               'Fail',
                           ], marker={'colors': ['rgba(25, 133, 37, 1.0)',
                                                 'rgba(133, 37, 25, 1.0)']},
                           domain={'x': [.70, .82], 'y': [.82, 1]},
                           name='Reads',
                           hoverinfo='labels+percent+name',
                           hole=.4,
                           textfont={'color': 'rgb(255,255,255)', 'size': 12},
                           type='pie',
                           visible=True,
                           showlegend=True)

    trace = go.Scatter(
        x=map(int, SortKeys),
        y=zip(*HoursDict)[0],
        name='ReadsPerHour',
        line=dict(
            color=('rgb(50,136,189)')),
        showlegend=False,
        hoverinfo='y',
        visible=True,
        xaxis='x5',
        yaxis='y1'
    )
    trace1 = go.Scatter(
        x=map(int, SortKeys),
        y=zip(*HoursDict)[1],
        name='BasesPerHour',
        line=dict(
            color=('rgb(253,174,97)')),
        showlegend=False,
        hoverinfo='y',
        visible=True,
        xaxis='x5',
        yaxis='y3'
    )
    trace2 = go.Scatter(
        x=map(int, SortKeys),
        y=map(float, zip(*HoursDict)[2]),
        name='MeanReadLengthPerHour',
        line=dict(
            color=('rgb(21,112,76)')),
        showlegend=False,
        visible=True,
        hoverinfo='y',
        xaxis='x5',
        yaxis='y5'
    )

    val = [[FailSum + PassSum], [PassSum], [FailSum], [TotalThroughput], [MeanReadLength], [SmallestRead], [BiggestRead]]
    log.debug('Plotting table...')
    table = go.Table(
        type='table',
        columnorder=[1, 2, 3, 4, 5, 6, 7],
        columnwidth=[5] * 7,
        visible=True,
        domain={'x': [0, 1], 'y': [0, .12]},
        header=dict(values=[['<b>READS TOTAL <br> NUMBER</b>'], ['<b>PASS</b>'], ['<b>FAIL</b>'], ['<b>TOTAL BP</b>'],
                            ['<b>MEAN READ <br>LENGTH (bp)</b>'], ['<b>SHORTEST READ (bp)</b>'],
                            ['<b>LONGEST READ (bp)</b>']],
                    line=dict(color='#7D7F80'),
                    align=['center'] * 7,
                    fill=dict(color='#00083E'),  # 00083e
                    font=dict(color='white', size=12), height=30),
        cells=dict(values=val, line=dict(color='#7D7F80'),
                   align=['center'] * 7,
                   fill=dict(color='#EDEDEE'),  # EDFAFF
                   font=dict(color='#506784', size=12), height=25))

    log.debug('Plotting length histogram...')
    LengthHist = go.Histogram(x=LengthsVec, nbinsx=100, nbinsy=10, opacity=0.75, showlegend=False, histnorm='percent',
                              name='% Reads', marker=dict(color='rgb(0, 200, 200)',
                                                          line=dict(width=0.5, color='rgba(255, 255, 255, 0.3)')),
                              xaxis='x4', yaxis='y4', hoverinfo='y+x+name')

    log.debug('Plotting quality histogram...')
    QualHist = go.Histogram(x=QualsVec, nbinsx=9, xbins=dict(start=6, end=14, size=1), nbinsy=10, opacity=0.75,
                            showlegend=False, histnorm='percent', name='% Reads',
                            error_y={'visible': True, 'thickness': 2, 'color': 'rgb(133, 37, 25)'},
                            marker=dict(color='rgb(0, 0, 255)', line=dict(width=0.5, color='rgba(255, 255, 255, 0.3)')),
                            xaxis='x6', yaxis='y6', hoverinfo='y+x+name')

    log.debug('Plotting GC histogram...')
    GCHist = go.Histogram(x=GCVec, nbinsx=14, xbins=dict(start=0, end=100, size=5), nbinsy=10, opacity=0.75,
                          showlegend=False, histnorm='percent', name='% Reads',
                          marker=dict(color='rgba(154, 0, 0, 1.0)',
                                      line=dict(width=0.5, color='rgba(255, 255, 255, 0.3)')),
                          xaxis='x7', yaxis='y7', hoverinfo='y+x+name')

    data = [trace, trace1, trace2, PieChartReads, LengthHist, QualHist, GCHist, table]

    layout = fig.layout
    fig = go.Figure(data=data, layout=layout)
    fig['layout'].update(title='Sperimental Summary', titlefont=dict(size=20), margin=dict(r=80, t=80, b=80, l=80),
                         hovermode='compare',
                         spikedistance=-1,
                         xaxis1=dict(visible=False, showspikes=False),
                         xaxis2=dict(visible=False, showspikes=False),
                         xaxis3=dict(visible=False, showspikes=False),
                         xaxis4=dict(title='Length', hoverformat='.2s', fixedrange=True, domain=[.55, 1]),
                         xaxis5=dict(visible=True, autotick=False, title='hours', showspikes=True, spikethickness=1.5,
                                     spikemode='across', spikecolor='grey', spikedash='solid'),
                         xaxis6=dict(visible=True, title='Quality', hoverformat='.2s', range=[5, 15], fixedrange=True,
                                     domain=[.55, 1]),
                         xaxis7=dict(visible=True, title='GC-Content', hoverformat='.2s', fixedrange=True,
                                     domain=[.3, .7]),
                         yaxis1=dict(fixedrange=True, title='Reads', anchor='x5', hoverformat='.2s'),
                         yaxis2=dict(visible=False),
                         yaxis3=dict(title='Bases', fixedrange=True, anchor='x5', hoverformat='.2s'),
                         yaxis4=dict(title='%', hoverformat='.2f', anchor='x4', fixedrange=False, domain=[.632, .732]),
                         yaxis5=dict(title='Mean Length', fixedrange=True, hoverformat='.2s'),
                         yaxis6=dict(title='%', hoverformat='.2f', anchor='x6', fixedrange=True, domain=[.445, .545]),
                         yaxis7=dict(title='%', hoverformat='.2f', anchor='x7', fixedrange=True),
                         legend=dict(x=1, y=1, xanchor='right'),
                         )
    config = {'displaylogo': False,
              'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud', 'zoom2d', 'pan2d', 'select2d', 'lasso2d',
                                         'hoverClosestCartesian', 'toggleHover', 'autoScale2d']}
    out_html_2 = os.path.join(work_dir, prefix + '_sequencing_summary.html')
    plotly.offline.plot(fig, filename=out_html_2, auto_open=False, show_link=False, config=config)
    HeatTrigger2(out_html_2)
    log.info('Finish')


def HeatTrigger(html_file):
    log.debug("Triggering heatmaps...")
    try:
        fh = file(html_file, 'r')
    except IOError:
        e = sys.exc_info()
        tb = traceback.format_exception(*e)
        log.error(
            "No such file or directory!\n\033[36m##################################\n\033[93m%s\033[32m%s\033[93m%s\033[36m##################################\n\033[0m" % (
                tb[0], tb[1], tb[2]))
        sys.exit(1)
    subject = fh.read()
    fh.close()
    Idpattern = '(?<=\<div id=")(.*?)"'
    try:
        FileId = re.search(Idpattern, subject).group(1)
    except:
        e = sys.exc_info()
        tb = traceback.format_exception(*e)
        log.error(
            "Can't find ID in html file!\n\033[36m##################################\n\033[93m%s\033[32m%s\033[93m%s\033[36m##################################\n\033[0m" % (
                tb[0], tb[1], tb[2]))
        sys.exit(1)
    if summary_flag == False:
        FixLine = "<script type='text/javascript'>var gd = document.getElementById('%(FileId)s');gd.on('plotly_click', function(data){var infotext = data.points.map(function(d){return (d.text)});var divname = data.points.map(function(d){return (d.data.name)}); var chan = data.points.map(function(d){return (d.text.split(' ')[1].replace('<br', ''))}); if(divname == 'ReadNHeat' || divname == 'BaseNHeat'){function range(start, stop, step){var a=[start], b=start;while(b<stop){b+=step;a.push(b)}return a;};var editindeces=range(10*(chan-1), 10*chan, 1); editindeces.pop(); var numeroPerguntas = 5122; var anyBoxesChecked = new Array(numeroPerguntas).fill(false); for (i in editindeces) { idx = parseInt(editindeces[i]); anyBoxesChecked[idx]=true}; var Newtitle = {title: infotext.toString()}; if (gd.data[0].visible === true) {var hv=[true, false]} else {var hv=[false, true]}; var exit = hv.concat(anyBoxesChecked); var update = {'visible': exit}; Plotly.update(gd ,update, Newtitle);};});gd.on('plotly_hover', function(data){var divname = data.points.map(function(d){return (d.data.name)});if(divname == 'Pore1' || divname == 'Pore2' || divname == 'Pore3' || divname == 'Pore4'){if(data.xvals){Plotly.Fx.hover(gd, {xval:data.xvals[0] }, ['x3y2', 'x3y3']);}}});</script>" % locals()  
    else:
        FixLine = "<script type='text/javascript'>var gd = document.getElementById('%(FileId)s');gd.on('plotly_click', function(data){var infotext = data.points.map(function(d){return (d.text)});var divname = data.points.map(function(d){return (d.data.name)}); var chan = data.points.map(function(d){return (d.text.split(' ')[1].replace('<br', ''))}); if(divname == 'ReadNHeat' || divname == 'BaseNHeat'){function range(start, stop, step){var a=[start], b=start;while(b<stop){b+=step;a.push(b)}return a;};var editindeces=range(4*(chan-1), 4*chan, 1); editindeces.pop(); var numeroPerguntas = 2050; var anyBoxesChecked = new Array(numeroPerguntas).fill(false); for (i in editindeces) { idx = parseInt(editindeces[i]); anyBoxesChecked[idx]=true}; var Newtitle = {title: infotext.toString()}; if (gd.data[0].visible === true) {var hv=[true, false]} else {var hv=[false, true]}; var exit = hv.concat(anyBoxesChecked); var update = {'visible': exit}; Plotly.update(gd ,update, Newtitle);};});gd.on('plotly_hover', function(data){var divname = data.points.map(function(d){return (d.data.name)});if(divname == 'Channel'){if(data.xvals){Plotly.Fx.hover(gd, {xval:data.xvals[0] }, ['x3y2', 'x3y3']);}}});</script>" % locals()
    SubLine = FixLine + '</body></html>'
    pattern = re.compile('</body></html>')
    result = pattern.sub(SubLine, subject)
    f_out = file(html_file, 'w')
    f_out.write(result)
    f_out.close()            
    try:
        shutil.move(html_file, os.path.join(work_dir, out_dir))
    except shutil.Error:
        log.error("Unable to move %s" % html_file)


def Bargen(PassReads, FailReads):
    trace1 = go.Bar(
        x=[FailReads],
        y=['Fail/Pass'],
        width=[0.4],
        name='Fail',
        text=str(FailReads) + "% Fail",
        orientation='h',
        showlegend=False,
        visible=False,
        hoverinfo="text",
        marker=dict(color='rgba(133, 37, 25, 1.0)'),
        xaxis='x4',
        yaxis='y4'
    )
    trace2 = go.Bar(
        x=[PassReads],
        y=['Fail/Pass'],
        width=[0.4],
        name='Pass',
        text=str(PassReads) + '% Pass',
        orientation='h',
        showlegend=False,
        visible=False,
        hoverinfo="text",
        marker=dict(color='rgba(25, 133, 37, 1.0)'),
        xaxis='x4',
        yaxis='y4'
    )
    bardata = [trace1, trace2]
    return bardata


def mux_fixer(lista):
    new_list = [None]*len(lista)
    if lista[0] != 0:
        new_list[0] = lista[0] 
    if lista[-1] != 0:   
        new_list[-1] = lista[-1]
    for x,u in enumerate(lista):
        if 0 < x < (len(lista)/2):  #-1
            if lista[x] != 0:
                new_list[x] = lista[x]    
            elif lista[x] == 0 and lista[x+1] == 0 and lista[x-1] != 0:
                new_list[x] = lista[x]
            elif lista[x] == 0 and lista[x+1] != 0 and lista[x-1] == 0:
                new_list[x] = None
        elif (len(lista)/2)-1 < x < len(lista)-1:
            if lista[x] != 0:
                new_list[x] = lista[x]    
            elif lista[x] == 0 and lista[x+1] == 0 and lista[x-1] != 0:
                new_list[x] = lista[x]
            elif lista[x] == 0 and lista[x+1] != 0 and lista[x-1] == 0:
                new_list[x] = lista[x]                                       
    return new_list                               


def ChannelPlotter_f(c, c_plo):
    if summary_flag == False:
        MuxList = [[0] * 48] * 8
        try:
            log.debug("Plotting channel %s stats..." % c)
            Bait = c_plo[str(c)][4]
            for kk, vv in Bait.items():
                if vv != []:
                    MuxLine1 = [0] * 48
                    MuxLine2 = [0] * 48
                    for u, z in vv[0].items():
                        MuxLine1[int(u)] = int(z[0])
                        MuxLine2[int(u)] = int(z[1])
                    MuxList[int(kk) - 1] = MuxLine1
                    MuxList[(int(kk) - 1) + 4] = MuxLine2
        except:
            pass 
        trace2Mux1 = go.Scatter(
            y=mux_fixer(MuxList[0]),
            name='Pore1',
            line=dict(
                color=('rgb(50,136,189)')),
            legendgroup='Pore1Chr' + str(c),
            showlegend=False,
            text=map(lambda x: str(x[0]) + '<br>' + str(x[1]) + 'h', zip(mux_fixer(MuxList[0]), hour_vec)),
            hoverinfo="text+name",
            fill='tozeroy',
            visible=False,
            xaxis='x3',
            yaxis='y2'
        )
        trace2Mux2 = go.Scatter(
            y=mux_fixer(MuxList[1]),
            name='Pore2',
            line=dict(
                color=('rgb(253,174,97)')),
            legendgroup='Pore2Chr' + str(c),
            showlegend=False,
            text=map(lambda x: str(x[0]) + '<br>' + str(x[1]) + 'h', zip(mux_fixer(MuxList[1]), hour_vec)),
            hoverinfo="text+name",
            fill='tozeroy',
            visible=False,
            xaxis='x3',
            yaxis='y2'
        )
        trace2Mux3 = go.Scatter(
            y=mux_fixer(MuxList[2]),
            name='Pore3',
            line=dict(
                color=('rgb(21,112,76)')),
            legendgroup='Pore3Chr' + str(c),
            showlegend=False,
            text=map(lambda x: str(x[0]) + '<br>' + str(x[1]) + 'h', zip(mux_fixer(MuxList[2]), hour_vec)),
            hoverinfo="text+name",
            fill='tozeroy',
            visible=False,
            xaxis='x3',
            yaxis='y2'
        )
        trace2Mux4 = go.Scatter(
            y=mux_fixer(MuxList[3]),
            name='Pore4',
            line=dict(
                color=('rgb(165,15,21)')),
            legendgroup='Pore4Chr' + str(c),
            showlegend=False,
            text=map(lambda x: str(x[0]) + '<br>' + str(x[1]) + 'h', zip(mux_fixer(MuxList[3]), hour_vec)),
            hoverinfo="text+name",
            fill='tozeroy',
            visible=False,
            xaxis='x3',
            yaxis='y2'
        )
        trace3Mux1 = go.Scatter(
            x=hour_vec,
            y=mux_fixer(MuxList[4]),
            name = 'Pore1',
            line=dict(
                color=('rgb(50,136,189)')),
            legendgroup='Pore1Chr' + str(c),
            visible=False,
            text=map(lambda x: str(x) + 'h', hour_vec),
            hoverinfo="y+text+name",
            fill='tozeroy',
            xaxis='x3',
            yaxis='y3'
        )
        trace3Mux2 = go.Scatter(
            x=hour_vec,
            y=mux_fixer(MuxList[5]),
            name='Pore2',
            line=dict(
                color=('rgb(253,174,97)')),
            legendgroup='Pore2Chr' + str(c),
            visible=False,
            text=map(lambda x: str(x) + 'h', hour_vec),
            hoverinfo="y+text+name",
            fill='tozeroy',
            xaxis='x3',
            yaxis='y3'
        )
        trace3Mux3 = go.Scatter(
            x=hour_vec,
            y=mux_fixer(MuxList[6]),
            name='Pore3',
            line=dict(
                color=('rgb(21,112,76)')),
            legendgroup='Pore3Chr' + str(c),
            visible=False,
            text=map(lambda x: str(x) + 'h', hour_vec),
            hoverinfo="y+text+name",
            fill='tozeroy',
            xaxis='x3',
            yaxis='y3'
        )
        trace3Mux4 = go.Scatter(
            x=hour_vec,
            y=mux_fixer(MuxList[7]),
            name='Pore4',
            line=dict(
                color=('rgb(165,15,21)')),
            legendgroup='Pore4Chr' + str(c),
            visible=False,
            text=map(lambda x: str(x) + 'h', hour_vec),
            hoverinfo="y+text+name",
            fill='tozeroy',
            xaxis='x3',
            yaxis='y3'
        )
        VisibleIndexes = range(10 * (c - 1), 10 * c)
        BarIdx = map(lambda x: x + 2, VisibleIndexes[-2:])
        fixIdx = map(lambda x: x + 2, VisibleIndexes[0:4])
        FailReads = 100
        if len(c_plo[str(c)]) > 0:
            FailReads = (int(c_plo[str(c)][2]) * 100) / int(c_plo[str(c)][0])
        PassReads = 100 - FailReads
        bardata = Bargen(PassReads, FailReads)
        tracedata = [trace2Mux1, trace2Mux2, trace2Mux3, trace2Mux4, trace3Mux1, trace3Mux2, trace3Mux3, trace3Mux4]
        outdata = tracedata + bardata
        Indexeddata = (outdata, fixIdx, BarIdx)
    else:
        try:
            MuxList1 = [0] * 48
            MuxList2 = [0] * 48
            Bait = c_plo[str(c)][4]
            for kk, vv in Bait.items():
                if vv != []:
                    nr, nb = vv
                    MuxList1[int(kk)] = int(nr)
                    MuxList2[int(kk)] = int(nb)
        except:
            pass
        trace2Mux1 = go.Scatter(
            y=mux_fixer(MuxList1),
            name='Channel',
            line=dict(
                color=('rgb(50,136,189)')),
            legendgroup='Pore1Chr' + str(c),
            showlegend=False,
            text=map(lambda x: str(x[0]) + '<br>' + str(x[1]) + 'h', zip(mux_fixer(MuxList1), hour_vec)),
            hoverinfo="text+name",
            fill='tozeroy',
            visible=False,
            xaxis='x3',
            yaxis='y2'
        )
        trace3Mux1 = go.Scatter(
            x=hour_vec,
            y=mux_fixer(MuxList2),
            name = 'Channel',
            line=dict(
                color=('rgb(50,136,189)')),
            legendgroup='Pore1Chr' + str(c),
            visible=False,
            text=map(lambda x: str(x) + 'h', hour_vec),
            hoverinfo="y+text+name",
            fill='tozeroy',
            xaxis='x3',
            yaxis='y3'
        )
        VisibleIndexes = range(4 * (c - 1), 4 * c)
        BarIdx = map(lambda x: x + 2, VisibleIndexes[-2:])
        fixIdx = map(lambda x: x + 2, VisibleIndexes[0:4])
        FailReads = 100
        if len(c_plo[str(c)]) > 0:
            FailReads = (int(c_plo[str(c)][2]) * 100) / int(c_plo[str(c)][0])
        PassReads = 100 - FailReads
        bardata = Bargen(PassReads, FailReads)
        tracedata = [trace2Mux1, trace3Mux1]
        outdata = tracedata + bardata
        Indexeddata = (outdata, fixIdx, BarIdx)          
    return Indexeddata


def result_plotting(ch_par):
    def rngrabber(t):
        try:
            Value = int(Nanoline[t][0])
        except:
            Value = 0
        return Value

    def bgrabber(t):
        try:
            Value = int(Nanoline[t][1])
        except:
            Value = 0
        return Value

    log.info('Plotting results...')
    log.debug("Preparing heatmaps...")
    NanoLayout = [['1', '2', '3', '4', '5', '6', '7', '8', '40', '39', '38', '37', '36', '35', '34', '33'],
                  ['9', '10', '11', '12', '13', '14', '15', '16', '48', '47', '46', '45', '44', '43', '42', '41'],
                  ['17', '18', '19', '20', '21', '22', '23', '24', '56', '55', '54', '53', '52', '51', '50', '49'],
                  ['25', '26', '27', '28', '29', '30', '31', '32', '64', '63', '62', '61', '60', '59', '58', '57'],
                  ['449', '450', '451', '452', '453', '454', '455', '456', '488', '487', '486', '485', '484', '483',
                   '482',
                   '481'],
                  ['457', '458', '459', '460', '461', '462', '463', '464', '496', '495', '494', '493', '492', '491',
                   '490',
                   '489'],
                  ['465', '466', '467', '468', '469', '470', '471', '472', '504', '503', '502', '501', '500', '499',
                   '498',
                   '497'],
                  ['473', '474', '475', '476', '477', '478', '479', '480', '512', '511', '510', '509', '508', '507',
                   '506',
                   '505'],
                  ['385', '386', '387', '388', '389', '390', '391', '392', '424', '423', '422', '421', '420', '419',
                   '418',
                   '417'],
                  ['393', '394', '395', '396', '397', '398', '399', '400', '432', '431', '430', '429', '428', '427',
                   '426',
                   '425'],
                  ['401', '402', '403', '404', '405', '406', '407', '408', '440', '439', '438', '437', '436', '435',
                   '434',
                   '433'],
                  ['409', '410', '411', '412', '413', '414', '415', '416', '448', '447', '446', '445', '444', '443',
                   '442',
                   '441'],
                  ['321', '322', '323', '324', '325', '326', '327', '328', '360', '359', '358', '357', '356', '355',
                   '354',
                   '353'],
                  ['329', '330', '331', '332', '333', '334', '335', '336', '368', '367', '366', '365', '364', '363',
                   '362',
                   '361'],
                  ['337', '338', '339', '340', '341', '342', '343', '344', '376', '375', '374', '373', '372', '371',
                   '370',
                   '369'],
                  ['345', '346', '347', '348', '349', '350', '351', '352', '384', '383', '382', '381', '380', '379',
                   '378',
                   '377'],
                  ['257', '258', '259', '260', '261', '262', '263', '264', '296', '295', '294', '293', '292', '291',
                   '290',
                   '289'],
                  ['265', '266', '267', '268', '269', '270', '271', '272', '304', '303', '302', '301', '300', '299',
                   '298',
                   '297'],
                  ['273', '274', '275', '276', '277', '278', '279', '280', '312', '311', '310', '309', '308', '307',
                   '306',
                   '305'],
                  ['281', '282', '283', '284', '285', '286', '287', '288', '320', '319', '318', '317', '316', '315',
                   '314',
                   '313'],
                  ['193', '194', '195', '196', '197', '198', '199', '200', '232', '231', '230', '229', '228', '227',
                   '226',
                   '225'],
                  ['201', '202', '203', '204', '205', '206', '207', '208', '240', '239', '238', '237', '236', '235',
                   '234',
                   '233'],
                  ['209', '210', '211', '212', '213', '214', '215', '216', '248', '247', '246', '245', '244', '243',
                   '242',
                   '241'],
                  ['217', '218', '219', '220', '221', '222', '223', '224', '256', '255', '254', '253', '252', '251',
                   '250',
                   '249'],
                  ['129', '130', '131', '132', '133', '134', '135', '136', '168', '167', '166', '165', '164', '163',
                   '162',
                   '161'],
                  ['137', '138', '139', '140', '141', '142', '143', '144', '176', '175', '174', '173', '172', '171',
                   '170',
                   '169'],
                  ['145', '146', '147', '148', '149', '150', '151', '152', '184', '183', '182', '181', '180', '179',
                   '178',
                   '177'],
                  ['153', '154', '155', '156', '157', '158', '159', '160', '192', '191', '190', '189', '188', '187',
                   '186',
                   '185'],
                  ['65', '66', '67', '68', '69', '70', '71', '72', '104', '103', '102', '101', '100', '99', '98', '97'],
                  ['73', '74', '75', '76', '77', '78', '79', '80', '112', '111', '110', '109', '108', '107', '106',
                   '105'],
                  ['81', '82', '83', '84', '85', '86', '87', '88', '120', '119', '118', '117', '116', '115', '114',
                   '113'],
                  ['89', '90', '91', '92', '93', '94', '95', '96', '128', '127', '126', '125', '124', '123', '122',
                   '121']]

    NanoBasePlate = []
    NanoPlate = []
    for i in NanoLayout:
        Nanoline = [ch_par[_] for _ in i]
        Exlin = map(lambda x: rngrabber(x), range(0, len(Nanoline)))
        NanoPlate.append(Exlin)
        Exlin = map(lambda x: bgrabber(x), range(0, len(Nanoline)))
        NanoBasePlate.append(Exlin)

    owncolorscale = [[0.0, 'rgb(0,0,0)'], [0.001, 'rgb(140,81,10)'], [0.1111111111111111, 'rgb(140,81,10)'],
                     [0.2222222222222222, 'rgb(191,129,45)'], [0.33333333333, 'rgb(223,194,125)'],
                     [0.44444444444, 'rgb(254,224,139)'], [0.5555555555555555, 'rgb(246,232,195)'],
                     [0.6666666666666666, 'rgb(199,234,229)'], [0.7777777777777778, 'rgb(128,205,193)'],
                     [0.8888888888888888, 'rgb(53,151,143)'], [1.0, 'rgb(1,102,94)']]

    xh = map(list, zip(*NanoLayout))
    yh = map(list, zip(*NanoPlate))
    zh = map(list, zip(*NanoBasePlate))
    hovertext=[[] for i in range(len(xh))]

    for x1, y1 in enumerate(xh):
        for x2, y2 in enumerate(y1):
            hovertext[x1].append('<b>CHANNEL:</b> {}<br /><b>RN:</b> {}<br /><b>BN:</b> {}'.format(y2, yh[x1][x2], zh[x1][x2]))

    trace = go.Heatmap(name="ReadNHeat", z=map(list, zip(*NanoPlate)), x=range(1, 33), y=range(1, 17),
                       text=hovertext, hoverinfo="text", xgap=23, ygap=8, colorscale=owncolorscale,
                       colorbar=dict(y=0.77, len=0.470, exponentformat="SI"))  # y=0.77,
    trace1 = go.Heatmap(name="BaseNHeat", z=map(list, zip(*NanoBasePlate)), x=range(1, 33), y=range(1, 17),
                        text=hovertext, hoverinfo="text", xgap=23, ygap=8, colorscale=owncolorscale,
                        colorbar=dict(y=0.77, len=0.470, exponentformat="SI"), visible=False)
    fig = plotly.tools.make_subplots(rows=13, cols=1, shared_xaxes=False,
                                     specs=[[{'rowspan': 6}], [None], [None], [None], [None], [None], [None],
                                            [{'rowspan': 2}], [None], [{'rowspan': 2}], [None], [None],
                                            [{'rowspan': 1}]],
                                     vertical_spacing=0.001, print_grid=False)

    log.debug("Heatmaps ready!")

    VisibleData = [False] * (513 * 10)

    updatemenus = list([
        dict(type="buttons",
             active=-1,
             buttons=list([
                 dict(label='Read Number',
                      method='update',
                      args=[{'visible': [True, False] + VisibleData},
                            {'annotations': trace}]),
                 dict(label='Base Pair Number',
                      method='update',
                      args=[{'visible': [False, True] + VisibleData},
                            {'annotations': trace1}]),
             ]),
             x=-0.03,
             y=0.99
             )
    ])

    fig.append_trace(trace, 1, 1)
    fig.append_trace(trace1, 1, 1)

    log.debug("Computing pore performance...")

    ChannelsInfoList = []
    X3Index = []
    X4index = []
    global hour_vec
    hour_vec = range(48)
    for c in range(1, 513):
        Res, x3i, x4i = ChannelPlotter_f(c, ch_par)
        ChannelsInfoList += Res
        X3Index += x3i
        X4index += x4i

    fig.data.extend(ChannelsInfoList)

    for fi in X3Index:
        fig.data[fi].update({'xaxis': 'x3'})

    for fi in X4index:
        fig.data[fi].update({'xaxis': 'x4'})

    log.debug("Plotting pore performance...")

    fig['layout'].update(title='Channel Productivity', spikedistance=-1, yaxis1=dict(range=[16.5, 0.5], visible=False, fixedrange=True),
                         xaxis1=dict(visible=False, fixedrange=True, range=[0.5, 32.5], hoverformat='.2s'),
                         yaxis2=dict(fixedrange=False, title='Reads', hoverformat='.2s', anchor='x3'),
                         xaxis2=dict(visible=False, fixedrange=True),
                         xaxis3=dict(autotick=False, title='hours', range=[0, 48], showspikes=True, spikethickness=1.5,
                                     spikemode='across', spikecolor='grey', spikedash='solid'),
                         yaxis3=dict(fixedrange=False, title='Base Pairs', hoverformat='.2s'),      
                         xaxis4=dict(visible=False, fixedrange=True, range=[0, 100]),
                         yaxis4=dict(fixedrange=True, visible=False),
                         legend=dict(x=-0.062, y=0.439, yanchor='yaxis2', xanchor='right'), updatemenus=updatemenus,
                         barmode='stack')  # y=0.450   range=[0, 1500000] range=[0, 600]
    fig.layout.update({'hovermode': 'compare'})
    config = {'displaylogo': False,
              'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud', 'zoom2d', 'pan2d', 'select2d', 'lasso2d',
                                         'hoverClosestCartesian', 'toggleHover', 'autoScale2d']}
    try:
        out_html = os.path.join(work_dir, prefix + '_pore_activity_map.html')
        plotly.offline.plot(fig, filename=out_html, show_link=False, auto_open=False, config=config)
    except:
        e = sys.exc_info()
        tb = traceback.format_exception(*e)
        log.error(
            "Something went wrong during plotting...\n\033[36m##################################\n\033[93m%s\033[32m%s\033[93m%s\033[36m##################################\n\033[0m" % (
                tb[0], tb[1], tb[2]))
        sys.exit(1)
    HeatTrigger(out_html)


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
        

def fastq_writer(tmp_dir):
    finalfile = os.path.join(work_dir, (prefix + '.fastq.gz'))
    file_list = os.listdir(tmp_dir)
    file_list.sort(key=lambda x: int(ntpath.basename(x)[4:-9]))
    with open(finalfile, 'wb') as outfile:
        for infile in map(lambda x: os.path.join(tmp_dir, x), file_list):
            shutil.copyfileobj(open(infile, 'rb'), outfile, 1024*1024*10)
    for b in map(lambda x: os.path.join(tmp_dir, x), file_list):
        try:
            os.remove(b)
        except Exception:
            log.warning("File %s doesn't exist!" % b)
    if not os.path.exists(os.path.join(work_dir, out_dir)):
        os.makedirs(os.path.join(work_dir, out_dir))            
    try:
        shutil.move(finalfile, os.path.join(work_dir, out_dir))
    except shutil.Error:
        log.error("Unable to move %s" % finalfile)
    shutil.rmtree(tmp_dir, ignore_errors=True)


def min_time_catcher(f):
    import h5py
    from dateutil import parser as par
    with h5py.File(f, 'r') as f5:        
        StartTime = f5['UniqueGlobalKey/tracking_id'].attrs['exp_start_time'].decode('UTF-8')
    dt = par.parse(StartTime)
    StartTConv = time.mktime(dt.timetuple())
    global RefTime 
    if RefTime is None or float(StartTConv) < float(RefTime):
        RefTime = str(StartTConv)


def summary_reader(filein):
    global ChannelDict
    ChannelDict = {str(il): () for il in range(1, 513)}
    ref_start = 1000000000
    lis=[[] for _ in range(512)]
    with open(filein) as fp:
        first_line = fp.readline()
        for line in fp:
            idx = int(line.split()[3])-1
            if float(line.split()[4]) < ref_start:
                ref_start = float(line.split()[4])   
            lis[idx].append(line.split()[3:-5])
    for l in lis:
        if l != []:
            l[:] = [l[z] for z in (y[0] for y in sorted(enumerate(zip(*l)[1]), key=lambda z: float(z[1])))]
            TimeVec = map(lambda x: int(time.strftime('%H', time.localtime(float(x) - float(ref_start)))), zip(*l)[1])
            for e in range(len(TimeVec)):
                hour = TimeVec[e]
                if e > 1:
                    if hour < TimeVec[e - 1]:
                        TimeVec[e] = hour + 24
            chr_list = zip(TimeVec, list(zip(*l)[3]), list(zip(*l)[-1]))
            ReadPerChannel = len(chr_list)
            BasesPerChannel = sum(map(int, zip(*chr_list)[1]))
            FailedReads = sum(float(i[2]) < 7.00 for i in chr_list)
            hour_table_productivity = {}
            for h in set(zip(*chr_list)[0]):
                hour_tab = [i for i in chr_list if i[0]==int(h)]
                n_reads = len(hour_tab)
                n_bases = sum(map(int,zip(*hour_tab)[1]))
                hour_table_productivity[str(h)] = (n_reads, n_bases)
                ChannelDict[str(l[0][0])] = (ReadPerChannel, BasesPerChannel, FailedReads, chr_list, hour_table_productivity)      
    return ChannelDict        


def deamon_table_reader(file_to_read):
    global deamon_thread
    deamon_thread = threading.Thread(target=summary_reader, args=file_to_read)
    deamon_thread.start()


def dir_surfer(file_folder):
    if summary_flag == False:
        for subdir, dirs, files in os.walk(file_folder):
            if dirs != []:
                my_subdir = map(lambda y: os.path.join(subdir, y), files)
                my_subdir.sort(key=lambda x: os.path.getmtime(x))
                files = map( lambda x: os.path.basename(x), my_subdir)
                for en, fi in enumerate(files):
                    if not fi.startswith('.'):
                        f = os.path.join(subdir, fi)
                        if en == 0 and int(os.path.basename(subdir)) == 0 and multiread_flag == False:
                            min_time_catcher(f)
                        yield f
            else:
                my_subdir = map(lambda y: os.path.join(subdir, y), files)
                my_subdir.sort(key=lambda x: os.path.getmtime(x))
                files = map( lambda x: os.path.basename(x), my_subdir)
                for en, fi in enumerate(files):
                    if not fi.startswith('.'):
                        f = os.path.join(subdir, fi)
                        if en == 0 and multiread_flag == False:
                            min_time_catcher(f)
                        yield f
    else:
        for subdir, dirs, files in os.walk(file_folder):
            for fi in files:
                if not fi.startswith('.') and fi.endswith('.fastq'):
                    f = os.path.join(subdir, fi)
                    if 'pass/' in f or 'fail/' in f:
                        yield f

def fast5module(work_dir, file_folder, prefix, tmp_dir, Fast_flag, size):    
    import fast5_reader as f5r 
    ChannelDict = {str(il): () for il in range(1, 513)}
    lis=[[] for _ in range(512)]

    for e,x in enumerate(dir_surfer(file_folder)):
        i = x.split('_')[-2]
        idx = int(i)-1
        lis[idx].append(x)

    if not os.path.exists(os.path.join(tmp_dir)):
        os.makedirs(os.path.join(tmp_dir))   

    Q = multiprocessing.JoinableQueue()
    R = multiprocessing.Queue()

    process_consumers_number = [Consumer(Q, R)
                                for ix in range(size)]

    for w in process_consumers_number:
        w.start()

    c_count = 0
    
    for e,chunk in enumerate(lis):
        if chunk != []:
            c_count += 1
            Q.put(f5r.parsing_func(chunk, e, tmp_dir, Fast_flag, RefTime))
        else: log.warn('Channel %s inactive!' % (e))    

    for i in range(size):
        Q.put(None)
    Q.join()

    while not Q.empty():
        Q.get()
   
    for r in range(c_count):
        result = R.get()
        ChannelDict[str(int(result[0])+1)] = result[1:6]
        
    if Fast_flag == True:
        fastq_writer(tmp_dir)

    result_plotting(ChannelDict)
    FastStats(ChannelDict)



def multifast5reader(work_dir, file_folder, prefix, tmp_dir, Fast_flag, size):    
    import multi_read_fast5_reader as mf5r
    ChannelDict = {str(il): () for il in range(1, 513)}
    lis=[[] for _ in range(512)]
    support_list=[]
    RefStart = None

    if not os.path.exists(os.path.join(tmp_dir)):
        os.makedirs(os.path.join(tmp_dir))   

    Q = multiprocessing.JoinableQueue()
    R = multiprocessing.Queue()

    process_consumers_number = [Consumer(Q, R)
                                for ix in range(size)]

    for w in process_consumers_number:
        w.start()

    c_count = 0

    for en,mff in enumerate(dir_surfer(file_folder)):
        c_count += 1
        False_flag = False if 'fail' not in mff else True
        Q.put(mf5r.mf5_reader(mff, en, tmp_dir, Fast_flag, False_flag))   

    for i in range(size):
        Q.put(None)
    Q.join()

    while not Q.empty():
        Q.get()
   
    for r in range(en):
        result = R.get()
        support_list.extend(result)
   
    for x in support_list:
        if type(x) is not float:
            idx = int(x[0][0])-1
            lis[idx].append(x[0][1:])
        else:
            RefStart=x     
    
    Q = multiprocessing.JoinableQueue()
    R = multiprocessing.Queue()

    process_consumers_number = [Consumer(Q, R)
                                for ix in range(size)]

    for w in process_consumers_number:
        w.start()

    c_count = 0

    for e,chunk in enumerate(lis):
        if chunk != []:
            c_count += 1
            Q.put(mf5r.channel_parser(chunk, e, RefStart))
        else: log.warn('Channel %s inactive!' % (e))    

    for i in range(size):
        Q.put(None)
    Q.join()

    while not Q.empty():
        Q.get()
   
    for r in range(c_count):
        result = R.get()
        ChannelDict[str(int(result[0])+1)] = result[1:]

    if Fast_flag == True:
        fastq_writer(tmp_dir)

    result_plotting(ChannelDict)
    FastStats(ChannelDict)
    os._exit(0)


def summary_module(fqfolder, summary_table_file, size):
    if size == 1:
        size = 2

    myList = list(dir_surfer(fqfolder))

    main_GC_list = []

    chunk_number = size - 1

    deamon_table_reader([summary_table_file])

    my_chunks = [myList[(i*len(myList))//chunk_number:((i+1)*len(myList))//chunk_number] for i in range(chunk_number)]
    Q = multiprocessing.JoinableQueue()
    R = multiprocessing.Queue()

    process_consumers_number = [Consumer(Q, R)
                                for ix in range(chunk_number)]

    for w in process_consumers_number:
        w.start()

    for chunk in my_chunks:
        Q.put(fastq_reader(chunk, fqfolder))

    for i in range(chunk_number):
        Q.put(None)
    Q.join()

    while not Q.empty():
        Q.get()

    for r in range(len(my_chunks)):
        result = R.get()
        main_GC_list.extend(result)

    while deamon_thread.isAlive() == True:
        time.sleep(1)

    result_plotting(ChannelDict)
    FastStats(ChannelDict, main_GC_list)         


def run(arguments):
    global summary_flag
    global work_dir
    global prefix
    global out_dir
    global multiread_flag
    if not len(arguments) > 4:
        sys.exit(1)   
    work_dir = arguments[0]
    file_folder = arguments[1]
    prefix = arguments[2]
    out_dir = arguments[3]
    ver = arguments[4]
    tmp_dir = os.path.join(work_dir, '.sts_temp')
    Fast_flag = eval(arguments[5])
    size = int(arguments[6])
    summary_flag = eval(arguments[7])
    multiread_flag = eval(arguments[9])
    if summary_flag == True:
        summary_table = eval(arguments[8])[0]

    log.info('Start Main Parser')
    verbosity(ver)

    if not os.path.exists(os.path.join(work_dir, out_dir)):
        os.makedirs(os.path.join(work_dir, out_dir))

    if summary_flag == False:
        global RefTime
        RefTime = None
        if multiread_flag == False: 
            fast5module(work_dir, file_folder, prefix, tmp_dir, Fast_flag, size)
        else:
            multifast5reader(work_dir, file_folder, prefix, tmp_dir, Fast_flag, size)    
    else:
        deamon_thread = None
        summary_module(file_folder, summary_table, size)

