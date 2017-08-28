# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:37:13 2017

@author: caspe
"""

from vispy import scene, app, geometry
from numpy import * 


class EditLine(scene.visuals.Line):
    def __init__(self, pos, *args, **kwargs):
        scene.visuals.Line.__init__(self, *args, **kwargs)
        
        self.unfreeze()
        self.markers = scene.visuals.Markers(pos = pos)
        self.pos = pos
        

class Canvas(scene.SceneCanvas):
    """ A simple test canvas for testing the EditLineVisual """

    def __init__(self):
        scene.SceneCanvas.__init__(self, keys='interactive',
                                   size=(800, 800), show = True)
        self.unfreeze()
        self.view = self.central_widget.add_view()
#        self.view.camera = scene.PanZoomCamera(rect=(-100, -100, 200, 200),
#                                               aspect=1.0)
        self.view.camera = 'arcball'
        # the left mouse button pan has to be disabled in the camera, as it
        # interferes with dragging line points
        # Proposed change in camera: make mouse buttons configurable
#        self.view.camera._viewbox.events.mouse_move.disconnect(
#            self.view.camera.viewbox_mouse_event)
        
#        self.events.mouse_press.connect(self.on_mouse_press)
#        self.events.mouse_release.connect(self.on_mouse_release,'mouse_release' )
        self.events.mouse_move.connect(self.on_mouse_press, 'mouse_move')
        self.events.mouse_wheel.connect(self.on_mouse_wheel, 'mouse_wheel')
#        self.view.camera._viewbox.events.mouse_move.disconnect(self.view.camera.viewbox_mouse_event)
#        self.view.add(self.line)
        vertices, faces, outline = geometry.create_plane(width=2, height=4,
                                                 width_segments=4,
                                                 height_segments=8,
                                                 direction='+y')

        plane = scene.visuals.Plane(width=2, height=4, width_segments=4,
                            height_segments=8, direction='+y',
                            vertex_colors=vertices['color'],
                            edge_color='k',
                            parent = self.view.scene)

        self.data = random.rand(100,3)
        
        self.line = scene.visuals.Line(pos = self.data, color = 'white')
#        t = scene.visuals.Markers(pos = self.data, face_color = 'white', symbol = 'disc', size = 10, parent = self.view.scene)
        geometry.MeshData(vertices = self.data)
        self.view.add(self.line)
 
#        self.show()
        
#        scene.visuals.GridLines(parent = self.view.scene)

        self.LEFT = 1
        self.RIGHT= 2
        self.MID  = 3

    def on_mouse_press(self,event):
        modifier = event.modifiers
#        print(modifier)
#        print(event.last_event.is_dragging)
        if event.button == self.RIGHT:
            
            
            print(self.view.scene.transform.imap(event.pos), self.view.scene.transform.map(event.pos))
            pos = self.view.scene.transform.imap(event.pos)[:3]
#            print(pos.shape, self.data.shape)
            self.data = vstack((self.data, pos[None,:]))
#            print(pos)
#            self.line.set_data(pos = self.data.T)
            self.line.set_data(pos = self.data, color = 'white')
            scene.visuals.Markers(pos = pos[None,:], face_color = 'white', symbol = 'disc', size = 10, parent = self.view.scene)
#            self.view.add(scene.visuals.Markers(pos = self.data, size = 10, face_color = 'white',symbol = 'disc'))
#            print(self.data.shape)
#            plane = scene.visuals.Plane(width=random.randint(1,100), height=4, width_segments=4,
#                            height_segments=8, direction='+y',
#                            vertex_colors= self.vertices['color'],
#                            edge_color='k',
#                            parent = self.view.scene)
#            self.line.draw()
            self.view.scene.update()
#            self.show()

        
    def on_mouse_wheel(self, event):
        '''
        -1 is zoom out; 1 is zoom in
        list of 2 ; vert; hor
        '''
        print(self.view.camera.get_state())
        state = self.view.camera.get_state() # returns dict with center fov and scale factor
        self.view.camera.set_state({'scale_factor': state['scale_factor'] + .1* event.delta[1]})
        print(self.view.camera.get_state())
        self.view.camera.update()
#        print(self.view.camera.zoom)
        
    
if __name__ == '__main__':
    win = Canvas()
    app.run()
