# -*- coding: utf-8 -*-
"""
Created on Sat May 13 04:07:54 2017

@author: caspe
"""
from pylab import *


def timeAnim(data, behavior, coordinates, rows = 2):
    print(behavior.shape, data.shape)
    dimMax = tuple(coordinates.max(axis = 0)+1)
    dimMin = coordinates.min(axis = 0)


    nT      = data.shape[-1]
    shape   = dimMax[:2]
    print(shape)
    storage = zeros(shape)
#    nLayers = dimMax[-1] - dimMin[-1]
    nLayers = 1
    cols    = int(ceil(nLayers / rows))

    print(cols * rows, nLayers)

    fig = figure()
    ax = []
    for i in range(rows):
        for j in range(cols):
            ax.append(subplot2grid((rows + 1,cols), (i,j)))
            ax[-1].axes.tick_params(labelbottom = 'off', labelleft = 'off')
            tight_layout(h_pad = 0, w_pad = 0)

    # set up the graphs
    h   = [ax[i].imshow(storage, vmin = -1, vmax = 1) for i in range(nLayers)]

    bAx = subplot2grid((rows + 1,cols),(rows, 0))
#    ax[-1].plot(behavior, 'b')
#    hb  = ax[-1].plot([0,0], [behavior.min(), behavior.max()],'--r)
    bAx.plot(behavior, 'b')
    [hb] = bAx.plot([0,0], [behavior.min(), behavior.max()],'--r')

    for i in range(nT):
        for layer in range(nLayers):
#            zidx = where(coordinates[:, -1] == layer)
#            print(zidx)
            storage.fill(0)
            for ci, di in zip(coordinates, data):
                print(ci, di[i])
                storage[ci[0], ci[1]] += di[i]
#            storage[coordinates[:,:2]] = di[i]
            h[layer].set_data(storage)
        hb.set_xdata([nT,nT])
        suptitle('t = {0}'.format(i))
        draw()
        pause(1e-7)

def plotLocations(locations, idx = None, j = None,  zfactor = 4):
    # TO DO: put this in visualize
    from vispy import scene, visuals
    from vispy.color import Color
    from numpy import median, array
    from PyQt5 import QtWidgets
#    print(locations.shape); assert 0
    window = QtWidgets.QMainWindow()
    j = 'Component {0}'.format(j)
    canvas = scene.SceneCanvas(keys='interactive', title = j, size=(800, 600), show=True)

    # Set up a viewbox to display the cube with interactive arcball
    view          = canvas.central_widget.add_view()
#    view.bgcolor  = '#efefef'
    view.bgcolor  = 'black'
    view.camera   = 'turntable'
    view.camera.elevation = 180
    view.camera.azimuth   = 90


#    # center the viewbox
    locations = locations.copy()                # deep copy: prevent aliasing
    locations[:,2] =  locations[:,2] * zfactor  # prevent pancake plot
    minmax = array([(i.min(), i.max()) for i in locations.T])
    view.camera.viewbox.camera.set_range(x = minmax[0,:], y = minmax[1,:], z = [0,0])
#        view.camera.viewbox.camera.elevation = 100000

#        view.camera.azimuth   = 320
#        view.scene.opacity = 0
    ax = scene.visuals.XYZAxis(parent = view.scene)
    tr = visuals.transforms.MatrixTransform()
    ax.transform = tr
    ax.transform.scale(locations.max(0))
    if idx != None:
        idx, r = idx
        r = (r - r.min()) / (r.max() - r.min())
        print(r.max(), r.min())
        scene.visuals.Markers(pos = locations[idx,:], face_color = cm.jet(r),
                              symbol = 'disc', parent = view.scene)


    shape = array(locations.shape)

    shape[-1] += 1
    c = ones(shape)
    c[...,-1] = .02 #alpha
#        print(c.shape)
    markers = scene.visuals.Markers(pos=locations, face_color =c,
                               symbol = 'disc', parent = view.scene, edge_color = None)
    window.centralWidget = canvas.native

    return view

from vispy import scene, color
from PyQt5.QtWidgets import QDialog, QFormLayout
class fastImshow(QDialog):
    '''
    Input data :
        (nfeature x nfeature x ntime (opt))
    '''
    def __init__(self, data, cmap = 'RdBu', use_slider = False):
        QDialog.__init__(self)
        self.setWindowTitle('fastImshow')
        canvas          = scene.SceneCanvas(keys = 'interactive', size = (800,800))

        # embed the imshow
        grid            = canvas.central_widget.add_grid()
        colormap        = color.get_colormap(cmap)

        # layout
        layout          = QFormLayout()
        # add slider if there is time data present
        if use_slider or len(data.shape) > 2:
            from PyQt5.QtWidgets import QSlider
            from PyQt5.QtCore import Qt
            slider = QSlider(Qt.Horizontal)
            slider.valueChanged.connect(self.updateAll)
            slider.setMinimum(0)
            slider.setMaximum(data.shape[-1] - 1)
            data = data.T # to always loop through time
            shape = data.shape[0:]
            use_slider = True

        cdata           = zeros((*data.shape, 4))
        print('Generating colors')
        from tqdm import tqdm
        for idx, column in enumerate(tqdm(data)):
            tmp             = 1 / (1 + exp(-column)) #sigmoid to map 0 1

            # 3d matrices need other reshape
            if use_slider:
                colors          = colormap.map(tmp).reshape((*column.shape,-1))
            else:
                colors          = colormap.map(tmp)
            cdata[idx, ...] = colors
        self.cdata = cdata
        if use_slider:
            tmp = cdata[0,...]
        else:
            tmp = cdata
        self.cfg = {'data': tmp}
        view            = grid.add_view(0,1, col_span = 4, row_span = 1)
        view.camera     = 'panzoom'
        print(tmp.shape)
        self.image = scene.visuals.Image(data = tmp, parent = view.scene)
        print('still')

        
#        assert 0


        xaxis = scene.AxisWidget(orientation='bottom',
                                axis_label='Component',
                                axis_font_size=40,
                                axis_label_margin=45,
                                tick_label_margin=5)


        yaxis = scene.AxisWidget(orientation='left',
                                 axis_label='Component',
                                 axis_font_size=40,
                                 axis_label_margin=60,
                                 tick_label_margin=5)

        grid.add_widget(yaxis, 0, 0, row_span = 1)
        grid.add_widget(xaxis, 1, 1, col_span = 4)
        right_padding = grid.add_widget(row=1, col=0, col_span = 4)
        right_padding.width_max = 120
        right_padding.height_max = 120
        print('here')
        # center the camera
        xr =  array([0, data.T.shape[0]]); yr = array([0, data.T.shape[1]])

        view.camera.viewbox.camera.set_range(x = xr, y = yr)


#        assert 0
        layout.addRow(canvas.native)
        xaxis.link_view(view)
        yaxis.link_view(view)
        if use_slider:
            layout.addRow(slider)
        print('here2')

        self.setLayout(layout)
        self.view = view
        self.show()

    def updateAll(self, value):
        '''
        if scroller present, this will change the values on plot
        '''
        image = self.cdata[value,...]
        self.image.set_data(image)
        self.view.update()


#%%




if __name__ == '__main__':

    import sys
    from PyQt5.QtCore import QCoreApplication
    app = QCoreApplication.instance()

    from numpy import *
    tmp = random.rand(100000,200)
#    app = QCoreApplication(sys.argv)
    if app is None:
        app = QtGui.QApplication(sys.argv)
        w = fastImshow(tmp)
    else:
        print('QApplication instance already exists: %s' % str(app))
#        w = ZebraFishViewer(c,d, [motorData, behaviorData])
        w = fastImshow(tmp)
#        w.show()
