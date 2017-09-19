
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 14:49:09 2017

@author: Casper van Elteren

The zebrafish viewer takes as input:
    data:   a dict containing containing a subdict
            Minimum input fields:
                datatype        : line or scatter (make name unique)
                data            : containing the data entries
            Data types:
                scatter:
                    Minimal input:
                            coordinates     : where in 3d to plot them
                    Optional inputs:
                            x      : time axis for figures
                            xlabel : corresponding figure label
                line:
                    Mininal input:
                        data
                    Optional data:
                        ...TBC

            
    
"""
# TODO: remove the coordinate dependency and put it into the scatter class, 
# make the scatter control group dependent on the scatter
 # TODO: fix this import mess 
from matplotlib import use
use('Qt5Agg')
from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt, QCoreApplication
from vispy import scene, visuals, app, gloo, color
from matplotlib.pyplot import *
from numpy import *
from matplotlib.pyplot import cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.figure import Figure
from tqdm import tqdm
#from seaborn import color_palette
import scipy

class EditLine(scene.visuals.LinePlot):
    '''
    Used for drawing the ROI
    To do:
        try except is uggly work around, linePlot cannot be initialized with empty data
    '''
    def __init__(self, *args, **kwargs):
        scene.visuals.LinePlot.__init__(self, *args, **kwargs)
        self.unfreeze()
    def getData(self):
        return self.data

    def setData(self, newData):
        # hacky way to allow init of line without initial data
        try:
            self.data = vstack((self.data, newData))
        except:
            self.data = newData
        self.visible = True

    def delData(self, kind):
        try:
            if kind == 'single':
                self.data = self.data[:-1, :]
            elif kind == 'all':
                self.__delattr__('data')
                self.visible = False
        except:
            self.visible = False

    def drawData(self):
        try:
            self.set_data(data = self.data)
        except:
            self.visible = False

    def getPolygon(self):
        '''
        Draws the polygon after confirming
        '''
        self.setData(self.data[0, :])
        self.drawData()

class Line(scene.SceneCanvas):
    def __init__(self, lineData, label = None):
        scene.SceneCanvas.__init__(self, keys = 'interactive', size = (300,300), show = False)
        self.unfreeze()
        grid          = self.central_widget.add_grid(margin = 10)
        # camera settings
        # view          = self.central_widget.add_view()
        view          = grid.add_view(1,1, col_span = 3)
        view.bgcolor  = 'black'
        view.camera   = 'panzoom'                                               # 2d view problem does  not fill entire screen
        view.padding  = 0
        xaxis = scene.AxisWidget(orientation='bottom',
                                axis_label='Time[step]',
                                axis_font_size=12,
                                axis_label_margin=45,
                                tick_label_margin=5)


        yaxis = scene.AxisWidget(orientation='left',
                                 axis_label='',
                                 axis_font_size=40,
                                 axis_label_margin=60,
                                 tick_label_margin=5)
        xaxis.height_max = 40
        yaxis.width_max = 40

        right_padding = grid.add_widget(row=1, col=2, row_span=1)
        right_padding.width_max = 5

        grid.spacing = 1
        [grid.add_widget(ax, *c) for ax, c in zip([xaxis, yaxis], ((2,1,1,3), (1,0)))]
        if label != None:
            title = scene.Label(label, color = 'white')
            title.height_max = 20
            grid.add_widget(title, 0, 1, col_span = 3)

        xaxis.link_view(view)
        yaxis.link_view(view)

        # add grid
        scene.visuals.GridLines(color = (1,1,1,.5), parent = view.scene)
        # tmp data
        if 1 in lineData.shape or len(lineData.shape) == 1:
            x        = list(range(len(lineData)))                                   #place holder for time
            lineData = vstack((x, lineData.flatten())).T
        line = scene.visuals.Line(pos = lineData,\
                                  width = 2,\
                                  color = '#0033cc',\
                                  parent = view.scene,\
                                  method = 'gl',\
                                  antialias = True)

        # update viewbox
        ylim = array([lineData.min(0)[-1], lineData.max(0)[-1]])
#        ylim += .7 *  ylim
        view.camera.viewbox.camera.set_range(x = [0, len(lineData)], y = ylim)

        tmp = vstack(([0,0], ylim))
        red = (1,0,0,.5)

        # set properties
        scrubCfg  = {'width': 3, 'face_color': (0,0,0,0), 'color':red,
                          'symbol': 'disc'}
        self.scrub     = scene.visuals.LinePlot(data = tmp.T, **scrubCfg, parent = view.scene)
        self.scrub.unfreeze()
        self.scrub.cfg = scrubCfg

        self.view      = view
        self.data      = lineData
        self.view      = view
        self.ylim      = ylim
        print('Constructed line plot')

class Scatter(scene.SceneCanvas):
    '''
    This object contains the 3d scatter data points
    '''
    def __init__(self, coordinates, data, cmap, parent, cluster = False, x = None, xlabel = ''):
        scene.SceneCanvas.__init__(self, keys = 'interactive',
                                   show = False,
                                   vsync = True,
                                   parent = parent)
#        cmap = 'summer' #TODO : change to default
        self.unfreeze() # add properties to class
        self.parent = parent
        if type(x)  == type(None):
            x = arange(len(data.T))
        self.x      = x
        self.xlabel = xlabel
        # camera, background setup
        view          = self.central_widget.add_view()
        view.bgcolor  = '#45454f' # 'black' #'#EFEBE7'#'#16161D'# 
        view.camera   = 'turntable'
        view.padding  = 0

        # init scatters
        self.coordinates         = coordinates.copy()
        self.t           = 0
        self.nC, self.nT = data.shape[:2]
        self.plotIdx     = arange(0, self.nC)

        #create colormap
        nColors       = 11
        colormap = color.get_colormap(cmap) # clip larger smaller 0,1
        self.colormap = colormap
        
        # plot contour TODO: clean this
#        from scipy import spatial
#        tmp = [unique(i) for i in coordinates.T] # obtain uniques
#        shape = (int(tmp[0].max()) + 1, int(tmp[1].max()) + 1)
#        m = zeros(shape) # create xy map
#        # make density of map
#        for x,y,z in coordinates:
##            if z > 10 : # leave out initial layers?
#             m[int(x),int(y)] += 1
#
#        m /= len(tmp[-1])     # normalize according to number of layers
#        l = sort(m[m>0].flatten()) # prevent spurious cells by picking top 90 percentile
##        print(len(l), len(m[m>0].flatten()))
#        n  = int(.7* len(l))  # percentile
#        i,j = where(m > l[n]) # get coordinates
#        tmp = vstack((i,j)).T # stack x to y
#        v = spatial.ConvexHull(tmp).vertices # create the convex hull
#        tmp = tmp[v,:]                       # border points
##        fig, ax = subplots(); ax.hist(m.flatten())
##        fig, ax = subplots(); ax.plot(*tmp.T)
#        scene.visuals.Line(pos = tmp, parent = view.scene, color = (1,1,1,.2)) # plot the outline
        print('Generating the color map')
        # this can be heavy on the memory
        # self.cdata  = zeros((*data.shape, 5), dtype = float16)
        # this seems to run the smoothest, other direct methods will yield memory glitches[not user friendly]
        # this however runs really slow
        # this implementation creates a linspace without outliers generating offshoots
#        self.cdata = computeColors(data, colormap) #added multithreading
        sigmoid = lambda x: 1 / ( 1 + exp( -x ))
        nColors = 41
        colorScale = linspace(0,1, nColors) # color scalings
        colors     = colormap.map(colorScale[:,None])
        
        self.cdata = zeros((*data.shape, 5))
        mindata = data.min(); maxdata = data.max()
        if  not abs(data.min()) >= 0 and not abs(data.max()) <= 1:
            print(allclose(data.min(), -1, 1e-2),  allclose(data.max(), 1, 1e-2))
            print('normalizing data')
#            data = sigmoid(data)
#            data = (data - data.min(0)) / (data.max(0) - data.min(0))
            for i, ti in enumerate(tqdm(data)):
                tmp = (ti - ti.min()) / (ti.max() - ti.min())
#                tmp = abs(tanh(ti))
                c   = colormap.map(tmp[:,None])
                if not cluster:
                    c[...,-1] = abs(tmp)     # update the alpha values (default 1)
                c = hstack((c, ti[:, None]))
                self.cdata[i,:,:] = c
        else:
            c         = array([colormap.map(abs(d)) for d in tqdm(data)]).reshape((*data.shape, -1))
            
            if not cluster:
                c[...,-1] = clip(abs(data), 0, 1)
            self.cdata[..., :-1] = c
            self.cdata[...,-1] = data
#        tmp = digitize(data, colorScale, right = True)

        y = np.max(coordinates[:,1])
        cbar = scene.visuals.ColorBar(cmap = cmap,
                               orientation = 'left',
                               size = (y, 10),
                               pos = (-20 , .5 * y),
                               clim = [0, 1],
                               parent = view.scene)
#        cbar.label.color = color.Color('white')
        self.cbar = cbar
        # add orientation aid
        ax = scene.visuals.XYZAxis(parent = view.scene)
        tr = visuals.transforms.MatrixTransform()
        ax.transform = tr
        ax.transform.scale(self.coordinates.max(0))
        # set markers
        # CONFIG FOR SCATTERS
        self.cfg = {'pos'       : self.coordinates,
                    'symbol'    :'disc',
                    'face_color':self.cdata[:, self.getT(), :4],
                    'edge_color': (1,1,1,.01),
                    'size'      : 5}
        self.mark = scene.visuals.Markers(**self.cfg, parent = view.scene)
        
        self.text = scene.visuals.Text(xlabel + str(self.x[self.t]), 'white', parent = view.scene, pos = (0, -20))

        #creating outline
        # create roi line
        self.roiLine = EditLine(data = array([[0], [0], [0]]).T,
                                symbol = 'disc',
                                color = 'white',
                                parent = view.scene)
        self.roiLine.visible = False
        # connect events
        self.events.mouse_press.connect(self.on_mouse_press)
        self.events.key_press.connect(self.on_key_press)

        # center the camera in the center of the grid
        minmax = array([(i.min(), i.max()) for i in self.coordinates.T])
        view.camera.viewbox.camera.set_range(x = minmax[0,:],\
                                             y = minmax[1,:],\
                                             z = minmax[2,:])

        # storage for cell time series
        self.cellPlots  = {}
        self.space = 0
        # update view
        self.view     = view
        print('Constructed scatter plot')

    def getT(self):
        'returns current time step of the data'
        return self.t
    def setT(self, value):
        'sets the current tiem step'
        if value < self.nT:
#            print('Value error\nRange {0}-{1}'.format(0, self.nT))
#        else:
            self.t = value
  
    def on_mouse_press(self, event):
        '''
        Ctrl + left mouse will add point
        Alt  + left mouse will destroy point
        Alt  + right mouse will destroy all points
        '''
        modifier = event.modifiers
        print(event.buttons)
        if len(modifier) == 1:
            modifier = modifier[0] # may be a tuple
#            print(event.last_event)
            if 'Control' in modifier.name and event.button == 1:
               position = self.view.scene.transform.imap(event.pos)[:3] # map from screen to world
               self.roiLine.setData(position[None,:])
               self.roiLine.drawData()
            elif 'Control' in modifier.name and event.button == 2:
                self.roiLine.delData(kind = 'single')
                self.roiLine.drawData()
            elif 'Control' in modifier.name and event.button == 3:
                self.roiLine.delData(kind = 'all')

        # plot time trace when clicking on a coordinate
        elif modifier != ():
            names = [i.name for i in modifier]
            print('keys active', names)
            # create new window
            if 'Control' in names and 'Shift' in names and event.button == 1:
                idx = self.findClosestToMouseClick(event)
                # separate plot
                fig, ax = subplots()
                fig.canvas.mpl_connect('close_event', self.figClose) # connect with close function
                fig.canvas.mpl_connect('button_press_event', self.determineActive) # connect with active determination
                
                # allow for multidimensional x 
                try:
                    ax.plot(self.x[idx, :], self.cdata[idx,:,-1])
                except:
                    ax.plot(self.x, self.cdata[idx,:,-1])
                ax.legend(['Cell {0}'.format(idx)])
                ax.set_xlabel(self.xlabel)
                
                self.plotActive   = len(self.cellPlots) # set active to the latest window
                print(self.plotActive)
                color             = cm.Accent(self.plotActive) # TODO: cannot have more than 255 plots open
                fig.suptitle('Selection {0}'.format(self.plotActive), color = color)
                cfg = {'pos': self.coordinates[[idx], :],
                       'face_color': color,
                       'size':       self.cfg['size'] * 3, # make marker slightly bigger to show what is in the selection
                       'edge_color': None}
                markers = scene.visuals.Markers(**cfg, parent = self.view.scene)
                show(block = False) # prevent from blocking the output

                self.cellPlots[self.plotActive] = {'fig': fig, 'markers' : markers, 'cfg' : cfg, 'indices': array(idx)}
                print('new selection added')
            # plot in active window
            elif 'Control' in names and 'Shift' in names and event.button == 2:
                # find point closest to mouse
                idx      = self.findClosestToMouseClick(event)
                activePlot = self.cellPlots[self.plotActive] # current active figure # TODO: make this according to what plot is currently viewed?
                # check if cell is already in the subset
                if idx in activePlot['indices']:
                    print('Index already in selection')
                # add it
                else:
                    ax       = self.cellPlots[self.plotActive]['fig'].axes[0] #current active plot
                    ax.plot(self.x, self.cdata[idx,:,-1])
                    leg      = ax.get_legend()
                    legTexts = leg.get_texts()
                    labels   = [i.get_text() for i in legTexts]
                    labels.append('Cell {0}'.format(idx))
                    ax.legend(labels)
                    # update indices
                    activePlot['cfg']['pos'] =  vstack((activePlot['cfg']['pos'], self.coordinates[[idx],:]))
                    activePlot['markers'].set_data(**activePlot['cfg'])
                    activePlot['indices'] = hstack((activePlot['indices'], idx))
                    # update the plot
                    draw()
        if event.button == 3:
            # TODO: remove this
            print('event ', event.pos)
            print('scene ', self.view.scene.transform.imap(event.pos))
            print('camera ', self.view.camera.transform.imap(event.pos))
            print(self.visual_at(event.pos))
            
    def on_key_press(self, event):
        '''
        Ctr  + s will autocomplete to first point
        Alt  + s will confirm and store selection
        '''
        # to do : key press to confirm selection
        # note : here no check for length modifier is made
        move = self.move
        if 'Control' in event.modifiers:
            if event.key.name == 'S':
                print('Saving')
                self.roiLine.getPolygon()
        elif 'Alt' in event.modifiers:
            if event.key.name == 'S':
                self.plotRoi(self.roiLine.getData())
        elif event.key.name == 'Left':
            move(self.getT() - 1)
        elif event.key.name == 'Right':
            move(self.getT() + 1)
        elif event.key.name == 'Space':
            self.parent.runAnimation(self.parent.autoTimeScroll)
        elif event.key.name == 'R':
            self.parent.runAnimation(self.parent.autoRotate)
            
    def move(self, value):
        if value >= 0 and value <= self.nT:
            self.parent.updateAll(value)
            
    def determineActive(self, event):
        '''Clicking on a window will make it active. Please note, dragging the window doesnt'''
        num = int(event.canvas.parent().windowTitle().split(' ')[-1]) # grab figure number
        for plot, properties in self.cellPlots.items():
            if properties['fig'].number == num:
                self.plotActive = plot
        
    def figClose(self, event):
        num = int(event.canvas.parent().windowTitle().split(' ')[-1]) # current active plot
        for plot, properties in self.cellPlots.items():
            if properties['fig'].number == num:
                properties['markers'].visible = False # maybe this needs a destroy
        self.view.scene.update()


    def findClosestToMouseClick(self, event):
        '''
        Find scene coordinate that is closest to the mouse click
        It indexes into what is currently plotted on the screen
        '''
#        state = self.view.camera.get_state()
#        azimuth = state['azimuth']/360
#        elevation = state['elevation']/360
#        viewBoxPos = self.view.scene.transform.imap(event.pos)[:3]
        viewBoxPos = self.view.scene.transform.imap(event.pos)[:2]
#        viewBoxPos = self.view.scene.transform.imap(viewBoxPos)[:3]

#        viewBoxPos[1:] *= array([azimuth, elevation])
        distance   = sqrt(sum((self.coordinates[self.plotIdx,:2] - viewBoxPos)**2,1))
        idx        = argmin(distance)
        idx        = self.plotIdx[idx]
        print('Scene coordinate {0}\nMarker coordinates {1}'.format(viewBoxPos, self.coordinates[idx, ...]))
        return idx

    def setPlotIdx(self, threshold):
        '''
        adjusts the plotIdx to specific alpha value threshold
        '''
        # for -1 or None plot entire range else plot value
        if  logical_and(threshold.min() >= 0, threshold.max() <= 1):
            
            alpha = self.cdata[:, self.getT(), -2]          # current alpha values in plot
            idx   = where(logical_and(alpha >= threshold[0], alpha <= threshold[1]))[0]            # disregard the rest
            # if none are bigger than threshold, don't show markers
            if len(idx) == 0:
                self.mark.visible = False
            else:
                self.plotIdx      = arange(0, len(alpha))[idx]             # update the plotIdx
                self.mark.visible = True
#        elif threshold == None or threshold < 0:
#            self.plotIdx = arange(self.nC)
        else:
            print('Threshold not understood\n Value should be in 0<=alpha<=1 or -1 for all')

    def updateCfg(self, cfg, update):
        for key, value in update.items():
            if key in cfg.keys():
                cfg[key] = value
            else:
                print('key not found in cfg, check the settings')

    def plotRoi(self, roi):
        from scipy.spatial import Delaunay
        print(self.plotIdx)
        try:
#            from visualize import plotLocations
#            plotLocations(self.locs[self.tmp,:])
            self.plotIdx = Delaunay(roi[:,:2]).find_simplex(self.coordinates[:, :2]) >= 0 #figure this out
        except scipy.spatial.qhull.QhullError: #tmp hack to assume you want slice in y direction
            tmp = []
            self.roi = roi
            for c in roi:
                tmp.append(logical_and(self.coordinates[:,0] < c[0], self.coordinates[:,1] < c[1]))

            idx = where(all(tmp))
            print(idx, tmp)
            self.plotIdx = idx
        #update config file
        updateCfg = {'face_color': self.cdata[self.plotIdx, self.getT(), :4],
                     'pos':  self.coordinates[self.plotIdx,:]}
        self.updateCfg(self.cfg, updateCfg)
        # show markers
        self.mark.set_data(**self.cfg)

# TODO : add menu with movie recording possiblities
class Menu(QMenu):
    def __init__(self, parent, *args, **kwargs):
        print(parent)
        QMenu.__init__(parent = parent, *args, **kwargs)
        
        
class ZebraFishViewer (QMainWindow):
    def __init__(self, data, cmap = 'summer'):
        # init window
        QMainWindow.__init__(self)
        screenResolution = QDesktopWidget().screenGeometry(-1)
        
        # screen resolution fix
        self.resize(screenResolution.width() * .25,screenResolution.height() * .5) # TODO: changes this
        self.setWindowTitle('zebrafish viewer')
        
        # add slider
        horSlider = QSlider(Qt.Horizontal)
        self.lab = QLabel(parent = horSlider)
        self.lab.setAlignment(Qt.AlignRight) # TODO: find a way to center this, ALignHCenter does not work
        self.timer       = app.Timer() # internal timer used for animation
        # add subplot lines
        layout = QGridLayout()
        widgets     = [] # init widgets
        # plots = {key : [] for key in data.keys()} # new storage vector
        plots  = {}
        for plotType, cfg in data.items():
            data = cfg['data']
            if 'scatter' in plotType:
                # create list entry if it does not exist
                if 'scatter' not in plots.keys():
                    plots['scatter'] = []
                widget = Scatter(**cfg, cmap = cmap, parent = self)
                plotType = 'scatter' #TODO: maybe clean this up
            elif 'line' in plotType:
                # make personal entry in plot indices
                if 'line' not in plots.keys():
                    plots['line'] = []
                # check if title of plot is given
                print(cfg.keys())
                if 'label' in cfg.keys():
                    plotTitle = cfg['label']
                else:
                    plotTitle = ''
                widget = Line(data, label = plotTitle)
                plotType = 'line'
            plots[plotType].append(widget) # new storage vector
            widgets.append(widget.native)
        self.plots = plots
        horSlider.setMinimum(0)
        
        # TODO: need a way to either check the data dimensions
        # and then estimate the length
        if len(data.shape) > 1:
            if data.shape[-1] != 1:    
                maxT = len(data.T) - 1
            else:
                maxT = len(data) - 1
        else:
            maxT = len(data) - 1
            
        horSlider.setMaximum(maxT)
#        self.lab.setText(str(horSlider.minimum()) * len(str(horSlider.maximum()))) # TODO: if you dont set the size the larger number will be clipped
        horSlider.setTickPosition(QSlider.TicksBothSides)
        horSlider.setTickInterval(maxT //10)
        self.nT = maxT
        # horSlider.setGeometry(10,10, 0,0)
        widgets.append(horSlider)
        # init control group
        if len(self.plots['scatter']) != 0:
            controlGroup = self.createControlGroup('Controls', self.plots['scatter'][0].coordinates)
            widgets.append(controlGroup)
        horSlider.valueChanged.connect(self.updateAll)
        # add widgets - generate and put it in layout
        step = 1
        for row, widget in enumerate(widgets):
            # %if row != 0: row += step
            # groupbox tends to be too big, fix it with conditional
            if type(widget) != QGroupBox:
                coordinate = (row, 0, step, 3)
            else:
                coordinate = (row, 0, step, 1)
            layout.addWidget(widget, *coordinate)
            
        # increase single plot size
        if 'line' not in plots.keys():
            if len(plots['scatter']) == 1:
                plots['scatter'][0].size = [i * 1.2 for i in plots['scatter'][0].size]
        self.horSlider   = horSlider
#        self.createMenu()
        # initialize layout
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)
#        self.setLayout(layout) # needed
        self.movie = []
        self.layout     = layout
        self.on_flatten()      # reset camera
        self.show()            # show the app
        print('Showing viewer')
        
    def createMenu(self):
        menu  =     QMenu('test', parent = self)
#        menu.addAction('test', self.about)
#        print('test')
        self.menuBar().addMenu(menu)
    def createControlGroup(self, title, coordinates):
        """
        Adds control group
        """
        # TODO: this needs a cleaning up, i.e. maybe controller function t hat makes the button 
        # create group add title
        zs = unique(coordinates[...,-1])
        controlGroup = QGroupBox(title)
        ZSliceControl = QLineEdit()
        ZSliceControl.objectNameChanged.connect(self.selectSlice)
        ZSliceControl.returnPressed.connect(self.selectSlice) # respond to enter commands
        
        sliceLabel   = QLabel('Slice\nmin {0} max {1}'.format(0, len(zs)-1))
        sliceButton  = QPushButton('Slice\nmin {0} max {1}'.format(0, len(zs)-1))
        sliceButton.clicked.connect(self.selectSlice)

        # flatten option
        flatButton = QPushButton('Flatten')
        flatButton.clicked.connect(self.on_flatten)
        flatButton.setAutoDefault(False)
        # reset all cells
        resetButton = QPushButton('Reset')
        resetButton.clicked.connect(self.on_reset)
        resetButton.setAutoDefault(False)
        # alpha control 
        alphaButton  = QPushButton('Alpha')
        alphaLabel   = QLabel('Alpha\nmin 0 max 1')
        alphaControl = QLineEdit()
#        tmp =  QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
#        alphaControl.setSizePolicy(tmp)
#        alphaControl.setMaximumHeight(40)
        alphaControl.objectNameChanged.connect(self.setAlpha)
        alphaControl.returnPressed.connect(self.setAlpha)  # respond to enter commands
        alphaButton.clicked.connect(self.setAlpha)
        
        frameControlLabel = QLabel('Time window slice\nmin {0} max {1}'.format(0, self.nT))
        frameControlButton = QPushButton('Time window slice\nmin {0} max {1}'.format(0, self.nT))
        frameControlButton.clicked.connect(self.playback)
        frameControl      = QLineEdit() # playback loop
        frameControl.objectNameChanged.connect(self.playback)
        frameControl.returnPressed.connect(self.playback)
        self.frameControl = frameControl
        
#        self.recordButton = BinaryButton('Record', parent = self)
#        self.recordButton.setStyleSheet("QPushButton:checked { background-color: red; }")
        
#        whatever = QLabel('test')
        
#        coordinates = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,2), (2,0)]
        widgets     = [(sliceLabel, ZSliceControl, alphaLabel, alphaControl),
                       (frameControlLabel, frameControl, 
                        resetButton, flatButton), 
                       ]
        
        controlLayout = QGridLayout()
        controlLayout.geometry().setWidth(10)
        controlLayout.geometry().setHeight(10)
        for row, widgetSet in enumerate(widgets):
            for column, widget in enumerate(widgetSet):
                controlLayout.addWidget(widget, row, column)
#        [controlLayout.addWidget(widget, *coordinate) for widget, coordinate in zip(widgets, coordinates)]
        controlGroup.setLayout(controlLayout)
        self.ZSliceControl = ZSliceControl
        self.alphaControl  = alphaControl
        return controlGroup
    

    def setAlpha(self):
        '''Controls the alpha parameter of the scatter plots (if present)'''
        import re
        text = self.alphaControl.text()
        l = re.findall("[-+]?\d*\.\d+|\d+", text)
        l = [float(i) for i in l] # convert to floats
        if len(l) < 2 : # append 1 if only one value is given alpha in [0,1]
            l.append(1)
        print('Setting alpha :', l)
        for scatter in self.plots['scatter']:
            t = scatter.getT()
            scatter.setPlotIdx(array(l))
        self.updateAll(t)
        
    def playback(self):
        text = self.getDigits(self.frameControl)
        # if empty reset the values
        if len(text) == 0:
            start = stop = None
        else:
            text = [int(i) if str.isdecimal(i) else None for i in text]
            start, stop = text[::2]
        print('Start / stop ', start, stop)
        # in case only one number is given
        if start == stop:
            stop = None
        from functools import partial
        self.autoTimeScroll = partial(self.autoTimeScroll, start = start,
                                      stop = stop)
#        print(func)
#        self.runAnimation(self.autoTimeScroll)
        
    def getDigits(self, control):
        # TODO use regex for this
        import itertools
        text =  ["".join(x) for _, x in itertools.groupby(control.text(), key = str.isdigit)]
        return text

    def selectSlice(self):
        '''
        Controls slice selection
        Checks for boundary conditions
        Slices are assumed to be integer values(!)
        '''
        import itertools
        # split text on numbers and non-numbers
        text = self.getDigits(self.ZSliceControl)

        # TODO: update this for independent scatter (multiple) (done?)
        for scatter in self.plots['scatter']:
            # get the correct dimensions from the scatter plot
            zs = unique(scatter.coordinates[...,-1])
            zdims = [zs.min(), 
                     zs.max()]
            # test boundary conditions
            try:
                # note this has the added benefit that if len(text) === 1, single slice is returned
                minZ, maxZ  = int(text[0]), int(text[-1])
                minZ = zs[minZ]
                maxZ = zs[maxZ]
                if minZ > maxZ:
                    print('Error, minimal slice should be smaller than max')
                #update markers
                elif zdims[0] <= minZ <= zdims[-1] and zdims[0] <= maxZ <= zdims[-1] :
                    # find intersection of locations and layer
                     idx = where(logical_and(scatter.coordinates[:,-1] >= minZ, scatter.coordinates[:,-1] <= maxZ))[0]
                     print('Slicing {0} - {1}'.format(minZ, maxZ))
                     if len(idx) == 0: # slice is empty
                         scatter.mark.visible = False
                     else:
                         scatter.mark.visible = True
                         # find overlap with entire set
                         allCell = tuple(arange(0, scatter.nC))
                         idx     = tuple(idx)
                         scatter.plotIdx = list(set(allCell) & set(idx))
                         coordinates = scatter.coordinates[scatter.plotIdx, :]
                         minmax      = array([coordinates.min(0), coordinates.max(0)]).T
                         scatter.view.camera.set_range(x = minmax[0,:], y = minmax[1,:], z = minmax[2,:])

                     self.updateAll(scatter.getT())

                else:
                    print('Error, z-range is {0}-{1}'.format(zdims[0], zdims[1]))

            except ValueError:
                print('Input not understood\nFormat should be min Z - max Z')
    def on_reset(self):
        # TODO: multiple scatter plots add functionality (done)
        # scatter = self.plots['scatter'][0]
        for scatter in self.plots['scatter']:
            scatter.plotIdx = arange(0, scatter.nC)
            self.updateAll(scatter.getT())

    def on_flatten(self):
        # TODO: multiple scatter plots add functionality (done)
        # scatter = self.plots['scatter'][0]
        for scatter in self.plots['scatter']:
            scatter.view.camera.elevation = 180
            scatter.view.camera.azimuth   = 180
        
    def autoRotate(self, event, increment = 1, incrementSign = 1):
        '''Starts app timer to begin rotation along the elvation direction'''
        for scatter in self.plots['scatter']:
            # this logic will go wrong when different scatter plots have a phase shift of 90 degrees
            # TODO: due to partial making a new function, it is not removed from the timer function
            # so you currently im manually removing it
            if scatter.view.camera.elevation >= 90:
                from functools import partial
                self.runAnimation(self.autoRotate)
                self.autoRotate = partial(self.autoRotate, incrementSign = -1)
                self.runAnimation(self.autoRotate)
                incrementSign = -1 # This needs to change
            elif scatter.view.camera.elevation <= -90:
                from functools import partial
                self.runAnimation(self.autoRotate)
                self.autoRotate = partial(self.autoRotate, incrementSign = 1)
                self.runAnimation(self.autoRotate)
            scatter.view.camera.elevation += increment * incrementSign
         
    def autoTimeScroll(self, event, start = None, stop = None, record = False):
        '''Automatically scrolls through the 2nd dimension of the data [assumed time]'''
        #TODO: currently the reset works only through reinvoking run animation
        # this needs to be fixed either by separating the checks of the function or coming up with another mechanism
        #TODO: only works when the scatter has focus
        
        # update very scatter plot
        for scatter in self.plots['scatter']:
            end   = scatter.nT      # default stopping values
            begin = scatter.getT()  # get current value
            # if start or stop is given set overwrite the default values
            if start != None:
                begin = start if begin < start else begin # only jump to start if t is smaller, else loops overwrites
            if stop != None:
                end = stop
            t = begin + 1 # step forward
            # reset the animation if out of bounds
            if t >= end:
                t = start if start != None else 0
            scatter.setT(t)
        single = self.updateAll(t)
        if record:
            self.movie.append(single)
        
    def runAnimation(self, func):
#        print('Starting animation')
        # check whether the timer is running or not
        if not self.timer.running:
            self.timer.functions = []
            self.timer.start()
        # hack to start the animation or not
        if func not in self.timer.functions:
            self.timer.connect(func)
            self.timer.functions.append(func)
        else:
            self.timer.disconnect(func)
            self.timer.functions.remove(func)
#        print(self.timer.functions)

    def updateAll(self, value):
        '''
        Binds together the slider value with the corresponding graphs
        '''
        #TODO: change this to be dynamic, i.e. check for object type? (done sorta)
        # update markers
        renders = []
        if 'scatter' in self.plots.keys():
            for scatter in self.plots['scatter']:
                scatter.setT(value)
                updateCfg = {'pos': scatter.coordinates[scatter.plotIdx,:], 'face_color' :  scatter.cdata[scatter.plotIdx, scatter.getT(), :4]}
                scatter.updateCfg(scatter.cfg, updateCfg)
                scatter.mark.set_data(**scatter.cfg)
                scatter.view.scene.update()
                scatter.text.text = scatter.xlabel + str(scatter.x[value])
                renders.append(scatter.view.scene.canvas.render())
        # update slider text
        self.lab.setText(str(value))
        self.horSlider.setValue(value) # hack to let it work with vispy canvas, i.e. left right when canvas is active will move this

        # update scrub lines
        if 'line' in self.plots.keys():
            for line in self.plots['line']:
                tmpData = vstack(([value, value], line.ylim)).T
                line.scrub.set_data(data = tmpData, **line.scrub.cfg)
        return renders
    
    def closeEvent(self, event):
        if 'scatter' in self.plots.keys():
            for scatter in self.plots['scatter']:
                try:
                    for cfg in scatter.cellPlots.values():
                        cfg['fig'].close()
                except:
                    pass
        self.timer.stop()
        self.deleteLater()

class NewQlineEdit(QLineEdit):
    def __init__(self, defaultText, *args, **kwargs):
        QLineEdit.__init__(self, *args, **kwargs)
        self.setText(defaultText)
        self.defaultText = defaultText
    def focusInEvent(self, event):
        self.clear()
    def focusOutEvent(self, event):
        if len(self.text()) == 0:
            self.setText(self.defaultText)
class BinaryButton(QPushButton):
    def __init__(self, *args, **kwargs):
        QPushButton.__init__(self, *args, **kwargs)
        self.state = 0
        self.clicked.connect(self.changeState)
    def changeState(self):
        self.state = 1 if self.state == 0 else 0
        [self.setStyleSheet("background-color: red") if self.state == 1 else 
         self.setStyleSheet("background-color: ")   ]
'''
Multiprocessing helping functions for computing the computing the colors
'''

def computeColorsMultiProc(ti, colormap):
    tmp = tanh(ti*1.5) # scal alpha value 0,1 for neg and pos
    s = tmp
    c   = colormap(s[:,None])
    c[...,-1] = abs(tmp)     # update the alpha values (default 1)
    c = hstack((c, ti[:, None]))
    return c
def computeColors(data, colormap):
    # sigmoid = lambda x: 1 / ( 1 + exp( -x ))
    from multiprocess import Pool, cpu_count
    from functools import partial
    func = partial(computeColorsMultiProc, colormap = colormap.map)
    pool = Pool(cpu_count())
    cdata = pool.map(func, data)
    return array(cdata)


def createQt(cfg, object = ZebraFishViewer):
    '''
    I don't fully understand how to correctly create QT object
    This seems to be a work around
    '''
    import sys
    from PyQt5.QtCore import QCoreApplication
    app = QCoreApplication.instance()
#    app = QCoreApplication(sys.argv)
    if app is None:
        app = QApplication(sys.argv)
        w = object(**cfg)
    else:
        print('QApplication instance already exists: %s' % str(app))
        w = object(**cfg)
    app.exec_()
    return w
    
#%% ONLY DEBUGGING STUFF: messy due to autoreload issue in python
if __name__ =='__main__':
    # example load
    nC = 500
    nT = 1000
    # load locations
#    print(coordinates.shape)
    cc = random.randint(0, 100, size = (nC,3))
    colormap = color.get_colormap('viridis')
#    data = sin(array([linspace(0, nT, nT) for i in range(nC)]) * 2 * pi) + random.randn(nC, nT) * 5 - 2
    ddd = random.rand(nC, nT) * 10 - 1
    behaviorData = random.randn(nT,1)
    motorData    = random.rand(*behaviorData.shape) * 15
    print(motorData.shape)
#    from multiprocessing import Pool, cpu_count
#    data = array(Pool(cpu_count()).map(calcCdata, data))
    # note weird names to prevent overwriting
    appp = QCoreApplication.instance()
    tmpScatterData = {'scatter main': {'data' : ddd, 'coordinates' : cc, 'xlabel': 'distance micrometer '}}
    tmpPlotData    = {'line motor' :
                                    {'data' : motorData,
                                     'label' : 'motor dummy data'}}
    cfg1 = {'data' :{**tmpScatterData, **tmpPlotData}} 
#    cfg1 = {'coordinates': coordinates, 'data' : {'scatter main' : {'data' : tmp, 'x': z, 'xlabel': 'distance $\mu$m'}}}
#    cfg1 = {'coordinates' : coordinates[:,:], 'data' : {'scatter main': {'data': abs(rr[:,None]), 'cluster' : False}}}
    w    = createQt(cfg1)
#    s = w.plots['scatter'][0]

##%%
#z = unique(coordinates[:,-1])
#from mpl_toolkits.mplot3d import Axes3D
#fig, ax = subplots(subplot_kw = {'projection':'3d'})
#x  = []
#y  = []
#for zi in z:
#    i = where(coordinates[:,-1] == zi)[0] #find indices
#    tmp = coordinates[i,:2]
#    hull = spatial.ConvexHull(tmp)
#    v = hull.vertices
#    verts = tmp[v,:]
#    k = spatial.distance.squareform(spatial.distance.pdist(verts))
#    k[eye(*k.shape) == 1] = nan
#    n = nanmin(k,0)
#    mu = n.mean()
#    s = n.std()
#    j = where(logical_and(n <= mu + 2 * s, n>= mu - 2 * s))[0]
#    v = v[j]
#    c = tmp[v,:]
#    l = ones((c.shape[0], 1))*zi
#    c = hstack((c, l))
#    x.append(c[:,0]); y.append(c[:,1])
#    ax.plot(*c.T, '.')
#x = array(x); y= array(y)
##%%
#tmp = [unique(i) for i in coordinates.T]
#m = zeros((tmp[0].max()+1, tmp[1].max()+1))
#for x,y,z in coordinates:
#    m[x,y] += 1
#m /= len(tmp[-1])
