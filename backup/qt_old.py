# -*- coding: utf-8 -*-
"""
Created on Wed May 31 14:49:09 2017

@author: caspe
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May 19 23:09:30 2017
to do:
    Create slider object v
    Create canvas for 3d plot v
    connect the update value of the slider to the data v

    -vanuit terminal krijg je error dat er een Qapplicatie gemaakt moet worden
    -de subplots zijn nog niet goed geframed
    -goede color scheme / de visualisatie zelfs ziet er nog niet goed uit
    - add flat button in control group (xy view)

@author: caspe
"""

from PyQt5.QtWidgets import QGroupBox, QWidget, QSlider, QInputDialog, QHBoxLayout, QDialog, QGridLayout, QLabel,QDesktopWidget, QApplication, QSpinBox, QLineEdit, QPushButton, QMainWindow
from PyQt5.QtCore import Qt, QCoreApplication
from vispy import scene, visuals, app, gloo
#app.use_app('pyqt5').
#import pylab as plt
#import matplotlib.pyplot as plt
from pylab import * 
from numpy import *
from matplotlib.pyplot import cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.figure import Figure
from tqdm import tqdm
from seaborn import color_palette


import scipy
#plt.ioff()
#print(plt.matplotlib.is_interactive())
# to do :
#    Qapplication moet nog worden opgezet
#    De plots in vispy zijn nu niet goed gezoomed en de resolutie vind ik nog niks


        
        

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
        scene.SceneCanvas.__init__(self, keys = 'interactive', size = (300,300), show = True)
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
        self.scrub     = scene.visuals.LinePlot(data = tmp.T,
                                                width = 3,
                                                face_color = (0,0,0,0),
                                                color = red,
                                                symbol = 'disc',
                                                parent = view.scene)
        self.view      =  view
        self.ylim      =  ylim



class Scatter(scene.SceneCanvas):
    '''
    This object contains the 3d scatter data points
    '''
    def __init__(self, coordinates, data):
        scene.SceneCanvas.__init__(self, keys = 'interactive', show = False, size = (400,600), vsync = True)
        self.unfreeze() # laat toe dat je andere modules kan toevoegen

        # camera, background setup
        view          = self.central_widget.add_view()
        view.bgcolor  = 'black'
        view.camera   = 'turntable'
        view.padding  = 0

        # init scatters
        self.locs         = coordinates.copy()
#        self.locs[...,-1]*= 4
        self.t           = 0
        self.nC, self.nT = data.shape[:2]
        self.plotIdx     = arange(0, self.nC)

        #create colormap
        nColors       = 11
#        self.colormap = array([cm.RdBu(int(i)) for i in linspace(0,255,nColors)])
        colormap = color_palette('RdBu_r', nColors)
#        colormap = [c + (1,) for c in colormap]
        self.colormap = array(colormap)

        print('Generating the color map')
        # rescale to 0,1 parse to colormap, and hack in the alpha as activity
        numStd = 1
        self.cdata  = zeros((*data.shape, 5), dtype = float16)

        # this seems to run the smoothest, other direct methods will yield memory glitches[not user friendly]
        # this however runs really slow
        # this implementation createas a linspace without outliers generating offshoots
        for i, ti in enumerate(tqdm(data.T)):
            m        = ti.mean()
            s        = ti.std()
            r        = linspace(m - numStd*s, m + numStd*s, nColors - 1)
            colorBin = digitize(ti, r)# bin data to match to colors
            c = array([self.colormap[i,:] for i in colorBin])
#            ti       = (ti - ti.min()) / (ti.max() - ti.min())
#            ti       = abs(ti[:, None]*2-1)
#            print(ti.min(), ti.max())
#            for i in ti: print(i)
#            ti =  abs(ti*2-1)
            c = hstack((c, abs(tanh(ti*1.5))[:, None], ti[:, None]))
            self.cdata[:,i,:] = c

        # add orientation aid
        ax = scene.visuals.XYZAxis(parent = view.scene)
        tr = visuals.transforms.MatrixTransform()
        ax.transform = tr
        ax.transform.scale(self.locs.max(0))
        # set markers
        # CONFIG FOR SCATTERS
        self.cfg = {'pos'       : self.locs,
                    'symbol'    :'disc',
                    'face_color':self.cdata[:, self.getT(), :4],
                    'edge_color': None,
                    'size'      : 10}
        self.mark = scene.visuals.Markers(**self.cfg, parent = view.scene)
        # create roi line
        self.roiLine = EditLine(data = array([[0], [0], [0]]).T, symbol = 'disc', color = 'white', parent = view.scene)
        self.roiLine.visible = False
        # connect events
        self.events.mouse_press.connect(self.on_mouse_press)
        self.events.key_press.connect(self.on_key_press)

        # center the camera in the center of the grid
        minmax = array([(i.min(), i.max()) for i in self.locs.T])
        view.camera.viewbox.camera.set_range(x = minmax[0,:],\
                                             y = minmax[1,:],\
                                             z = minmax[2,:])
        
        # storage for cell time series
        self.cellPlots  = {}
         
        # update view
        self.view     = view


    def getT(self):
        return self.t
    
    def setT(self, value):
        if value > self.nT:
            print('Value error\nRange {0}-{1}'.format(0, self.nT))
        else:
            self.t = value

    def on_mouse_press(self, event):
        '''
        Ctrl + left mouse will add point
        Alt  + left mouse will destroy point
        Alt  + right mouse will destroy all points
        '''
        modifier = event.modifiers
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
            
            # create new window
            if 'Control' in names and 'Alt' in names and event.button == 1:
                idx = self.findClosestToMouseClick(event)
                # separate plot
                fig, ax = subplots()
                fig.canvas.mpl_connect('close_event', self.figClose)
                ax.plot(self.cdata[idx,:,-1])
                ax.legend(['Cell {0}'.format(idx)])
                ax.set_xlabel('Time [step]')
    
                self.plotActive   = len(self.cellPlots) # set active to the latest window
                color             = cm.prism(self.plotActive)
                fig.suptitle('Selection {0}'.format(self.plotActive), color = color)
                cfg = {'pos': self.locs[[idx], :],
                       'face_color': color,
                       'size':       self.cfg['size'] * 3,
                       'edge_color': None}
                markers = scene.visuals.Markers(**cfg, parent = self.view.scene)
                
                self.cellPlots[self.plotActive] = {'fig': fig, 'markers' : markers, 'cfg' : cfg, 'indices': array(idx)} 
                print('new selection added')
                
            # plot in active window
            elif 'Control' in names and 'Alt' in names and event.button == 2:
                # find point closest to mouse
                idx      = self.findClosestToMouseClick(event)
                # check if cell is already in the subset
                if idx in self.cellPlots[self.plotActive]['indices']:
                    print('Index already in selection')
                # add it 
                else:
                    ax       = self.cellPlots[self.plotActive]['fig'].axes[0] #change this to current active
                    ax.plot(self.cdata[idx,:,-1])
                    leg      = ax.get_legend()
                    legTexts = leg.get_texts()                
                    labels   = [i.get_text() for i in legTexts]
                    labels.append('Cell {0}'.format(idx))
                    ax.legend(labels)
                    # update indices 
                    self.cellPlots[self.plotActive]['cfg']['pos'] =  vstack((self.cellPlots[self.plotActive]['cfg']['pos'], self.locs[[idx],:]))
                    self.cellPlots[self.plotActive]['markers'].set_data(**self.cellPlots[self.plotActive]['cfg'])
                    self.cellPlots[self.plotActive]['indices'] = hstack((self.cellPlots[self.plotActive]['indices'], idx))
                    # update the plot
                    draw() 
                
    def figClose(self, event):
        num = int(event.canvas.parent().windowTitle().split(' ')[-1])
        for plot, properties in self.cellPlots.items():
            if properties['fig'].number == num:
                properties['markers'].visible = False # maybe this needs a destroy
        self.view.scene.update()
            
#                ax.legend(idx.tolist())
    def findClosestToMouseClick(self, event):
        viewBoxPos = self.view.scene.transform.imap(event.pos)[:3]
        distance   = sqrt(sum((self.locs - viewBoxPos)**2,1))
        idx        = argmin(distance)
        return idx
    
    def on_key_press(self, event):
        '''
        Ctr  + s will autocomplete to first point
        Alt  + s will confirm and store selection
        '''
        # to do : key press to confirm selection
        # note : here no check for length modifier is made
        if 'Control' in event.modifiers:
            if event.key.name == 'S':
                print('Saving')
                self.roiLine.getPolygon()
        elif 'Alt' in event.modifiers:
            if event.key.name == 'S':
                self.plotRoi(self.roiLine.getData())
    
    def setPlotIdx(self, threshold):
        '''
        adjusts the plotIdx to specific alpha value threshold
        '''
        print(threshold)
        # for -1 or None plot entire range else plot value
        if  0<=threshold <= 1:
            print('wrong')
            alpha = self.cdata[:, self.getT(), -1] #current alpha values in plot
            idx   = where(alpha >= threshold)[0]          # disregard the rest
            
            # if none are bigger than threshold, don't show markers
            if len(idx) == 0:
                self.mark.visible = False
            else:
                print(self.plotIdx)
                self.plotIdx = arange(0, self.nC)[idx]              # update the plotIdx
                print(self.plotIdx)
                self.mark.visible = True
        elif threshold == None or threshold < 0:
            self.plotIdx = arange(self.nC)
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
            self.plotIdx = Delaunay(roi[:,:2]).find_simplex(self.locs[:, :2]) >= 0 #figure this out
        except scipy.spatial.qhull.QhullError: #tmp hack to assume you want slice in y direction
            print('here')
            tmp = []
            self.roi = roi
            for c in roi:
                tmp.append(logical_and(self.locs[:,0] < c[0], self.locs[:,1] < c[1]))

            idx = where(all(tmp))
            print(idx, tmp)
            self.plotIdx = idx
        #update config file
        updateCfg = {'face_color': self.cdata[self.plotIdx, self.getT(), :4],
                     'pos':  self.locs[self.plotIdx,:]}
        self.updateCfg(self.cfg, updateCfg)
        # show markers
        self.mark.set_data(**self.cfg)
#        self.mark.set_data(pos = self.locs[self.plotIdx , :],
#                                           face_color = self.cdata[self.plotIdx , self.getT(), :4])


class ZebraFishViewer (QDialog):
    def __init__(self, coordinates, data, motorAndBehavior = None):
        # init window
        QDialog.__init__(self)

#        app = QApplication([])
        screenResolution = QDesktopWidget().screenGeometry(-1)
        
        self.resize(800,screenResolution.height()- 100)
        self.setWindowTitle('zebra fish viewer')
        # attributes
        self.coordinates = coordinates
        self.dims        = coordinates.max(0)
        self.plotIdx     = range(coordinates.shape[0]) # used for plotting Z-slicing etc
        # add slider
        horSlider = QSlider(Qt.Horizontal)
        horSlider.setMinimum(0)
        horSlider.setMaximum(data.shape[-1] - 1)

        horSlider.setTickPosition(QSlider.TicksBothSides)
#        horSlider.setSingleStep(1)
        horSlider.setTickInterval(data.shape[-1] //10)
        horSlider.setGeometry(10,10, 0,0)
        self.lab = QLabel(parent = horSlider)
        self.lab.setText(str(horSlider.minimum()) * len(str(horSlider.maximum()))) #tmp hacky way; if you dont set the size the larger number will be clipped
        self.lab.setAlignment(Qt.AlignHCenter)
        self.lab.setStyleSheet("QLabel {background-color: red;}")
#        self.horSlider = horSlider
        self.scatter   = Scatter(coordinates, data)
        widgets        = [self.scatter.native]
        # add subplot lines
        if motorAndBehavior != None:
            lineCanvs  = []
            self.lines = []
            
            for label, data in motorAndBehavior.items():                
                line = Line(data, label = label)
                lineCanvs.append(line.native)
                self.lines.append(line)
            [widgets.append(l) for l in lineCanvs]
        else:
            self.scatter.size = array(self.scatter.size) * 2
        widgets.append(horSlider)
        # init control group
        self.createControlGroup('Controls')
        widgets.append(self.controlGroup)
        horSlider.valueChanged.connect(self.updateAll)
            # gridlayout (sim. matplotlib)
        layout = QGridLayout()
        # add widgets
        for row, widget in enumerate(widgets):
            if type(widget) != QGroupBox:
                coordinate = (row, 0, 1, 3)
            else:
                coordinate = (row, 0, 1, 1)
            
            layout.addWidget(widget, *coordinate)

        self.widgets    = widgets
        # initialize layout
        self.layout     = layout
        self.setLayout(layout) # needed
        self.on_flatten()
        self.show()

    def createControlGroup(self, title):
        # create group add title
        self.controlGroup = QGroupBox(title)
#        label = QLabel('Start-stop slice:')

        ZSliceControl = QLineEdit()
        ZSliceControl.objectNameChanged.connect(self.selectSlice)
        sliceButton        = QPushButton('Slice\n min {0} max{1}'.format(0, self.dims[-1]))
        sliceButton.clicked.connect(self.selectSlice)


        # flatten option
        flatButton = QPushButton('Flatten')
        flatButton.clicked.connect(self.on_flatten)
        
        # reset all cells
        resetButton = QPushButton('Reset')
        resetButton.clicked.connect(self.on_reset)
        
        alphaButton  = QPushButton('Alpha')
        alphaControl = QLineEdit()
        alphaControl.geometry().setWidth(10)
        alphaControl.objectNameChanged.connect(self.setAlpha)
        alphaButton.clicked.connect(self.setAlpha)
        
        coordinates = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,2)]
        widgets     = [ZSliceControl, sliceButton, flatButton, alphaControl, alphaButton, resetButton]
        controlLayout = QGridLayout()
        controlLayout.geometry().setWidth(10)
        controlLayout.geometry().setHeight(10)
        [controlLayout.addWidget(widget, *coordinate) for widget, coordinate in zip(widgets, coordinates)]
        self.controlGroup.setLayout(controlLayout)
        self.ZSliceControl = ZSliceControl
        self.alphaControl  = alphaControl
        
        
    def setAlpha(self):
        text = float(self.alphaControl.text())
        print(text)
        self.scatter.setPlotIdx(text)
        self.updateAll(self.scatter.getT())
    
    
    def getDigits(self, control):
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
        # get the correct dimensions
        zdims = [0, self.dims[-1]]
        # test boundary conditions
        try:
            # note this has the added benefit that if len(text) === 1, single slice is returned
            minZ, maxZ  = int(text[0]), int(text[-1])
            if minZ > maxZ:
                print('Error, minimal slice should be smaller than max')
            #update markers
            elif zdims[0] <= minZ <= zdims[-1] and zdims[0] <= maxZ <= zdims[-1] :
                self.scatter.plotIdx = where(logical_and(self.scatter.locs[...,-1] >= minZ, self.scatter.locs[...,-1] <= maxZ))[0]
#                self.scatter.mark.set_data(pos = self.scatter.locs[self.scatter.plotIdx , :],
#                                           face_color = self.scatter.cdata[self.scatter.plotIdx , self.scatter.getT(), :4])
         
                coordinates = self.scatter.locs[self.scatter.plotIdx, :]
                minmax      = array([coordinates.min(0), coordinates.max(0)]).T
                self.scatter.view.camera.set_range(x = minmax[0,:], y = minmax[1,:], z = minmax[2,:])
                self.updateAll(self.scatter.getT())

            else:
                print('Error, z-range is {0}-{1}'.format(zdims[0], zdims[1]))

        except ValueError:
            print('Input not understood\nFormat should be min Z - max Z')
    def on_reset(self):
        self.scatter.plotIdx = arange(0, self.scatter.nC)
        self.updateAll(self.scatter.getT())

    def on_flatten(self):
        self.scatter.view.camera.elevation = 180
        self.scatter.view.camera.azimuth   = 180

  

    def updateAll(self, value):
        '''
        Binds together the slider value with the corresponding graphs
        '''
        # update markers
        self.scatter.setT(value)
        updateCfg = {'pos': self.scatter.locs[self.scatter.plotIdx,:], 'face_color' :  self.scatter.cdata[self.scatter.plotIdx, self.scatter.getT(), :4]}
        
        self.scatter.updateCfg(self.scatter.cfg, updateCfg)
        self.scatter.mark.set_data(**self.scatter.cfg)

#        self.scatter.cfg = self.scatter.updateCfg(self.scatter.cfg, ['face_color', self.scatter.cdata[self.scatter.plotIdx, self.scatter.getT(), :4]])
#        self.scatter.cfg = self.scatter.updateCfg(self.scatter.cfg, ['pos', self.scatter,locs[self.scatter.plotIdx, :]])
#        self.scatter.mark.set_data(**self.scatter.cfg)
        # update slider text
        self.lab.setText(str(value))
        # update scrub lines
        [line.scrub.set_data(data = vstack(([value, value], line.ylim)).T,
                             face_color = (0,0,0,0), width = 3) for line in self.lines]
        
    def closeEvent(self, event):
        self.deleteLater()
        self.close()
        QCoreApplication().quit()

if __name__ =='__main__':

    nC = 100
    nT = 1000
    # load locations

    c = random.randint(0, 100, size = (nC,3))

#    data = sin(array([linspace(0, nT, nT) for i in range(nC)]) * 2 * pi) + random.randn(nC, nT) * 5 - 2
    d = random.randn(nC, nT) * 10 - 1
    behaviorData = random.randn(nT,1)
    motorData    = random.rand(*behaviorData.shape) * 15
    print(motorData.shape)
#    from multiprocessing import Pool, cpu_count
#    data = array(Pool(cpu_count()).map(calcCdata, data))

    app = QCoreApplication.instance()
    
   
#    app = QCoreApplication(sys.argv)
    if app is None:
        app = QtGui.QApplication(sys.argv)
        w = ZebraFishViewer(c,data, [motorData, behaviorData])
        w.show()
    else:
        print('QApplication instance already exists: %s' % str(app))
#        w = ZebraFishViewer(c,d, [motorData, behaviorData])
        w = ZebraFishViewer(c,d, None)
        w.show()
