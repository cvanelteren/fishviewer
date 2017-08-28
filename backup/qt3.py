# -*- coding: utf-8 -*-
"""
Created on Fri May 26 08:14:00 2017

@author: caspe
"""

"""
This example shows how to display 3D objects.
You should see a colored outlined spinning cube.
"""

import numpy as np
from vispy import app, gloo
from vispy.util.transforms import perspective, translate, rotate

vert = """
// Uniforms
// ------------------------------------
uniform   mat4 u_model;
uniform   mat4 u_view;
uniform   mat4 u_projection;
uniform   vec4 u_color;

// Attributes
// ------------------------------------
attribute vec3 a_position;
attribute vec4 a_color;
attribute vec3 a_normal;

// Varying
// ------------------------------------
varying vec4 v_color;

void main()
{
    v_color = a_color * u_color;
    gl_Position = u_projection * u_view * u_model * vec4(a_position,1.0);
}
"""


frag = """
uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_normal;

uniform vec3 u_light_intensity;
uniform vec3 u_light_position;

varying vec3 v_position;
varying vec3 v_normal;
varying vec4 v_color;

void main()
{
    gl_FragColor = v_color;
}
"""


# -----------------------------------------------------------------------------
def cube(num_of_cubes):
    """
    Build vertices for a colored cube.

    V  is the vertices
    I1 is the indices for a filled cube (use with GL_TRIANGLES)
    I2 is the indices for an outline cube (use with GL_LINES)
    """

    for i in range(0,num_of_cubes):
        # Vertices positions
        v = np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1],
             [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, -1]],dtype=np.float32)

        v[:,0]=v[:,0]+2.
        v[:,1]=v[:,1]+2.
        v[:,2]=v[:,2]+2.

        # Face Normals
        n =np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0],
             [-1, 0, 1], [0, -1, 0], [0, 0, -1]],dtype=np.float32)
        # Vertice colors
        c = np.array([[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1],
             [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]],dtype=np.float32)

        V_aux = np.array([(v[0], n[0], c[0]), (v[1], n[0], c[1]),
                      (v[2], n[0], c[2]), (v[3], n[0], c[3]),
                      (v[0], n[1], c[0]), (v[3], n[1], c[3]),
                      (v[4], n[1], c[4]), (v[5], n[1], c[5]),
                      (v[0], n[2], c[0]), (v[5], n[2], c[5]),
                      (v[6], n[2], c[6]), (v[1], n[2], c[1]),
                      (v[1], n[3], c[1]), (v[6], n[3], c[6]),
                      (v[7], n[3], c[7]), (v[2], n[3], c[2]),
                      (v[7], n[4], c[7]), (v[4], n[4], c[4]),
                      (v[3], n[4], c[3]), (v[2], n[4], c[2]),
                      (v[4], n[5], c[4]), (v[7], n[5], c[7]),
                      (v[6], n[5], c[6]), (v[5], n[5], c[5])]
        )
        I1_aux = np.resize(np.array([0, 1, 2, 0, 2, 3], dtype=np.uint32), 6 * (2 * 3))
        I1_aux += np.repeat(4 * np.arange(2 * 3, dtype=np.uint32), 6)

        I2_aux = np.resize(
            np.array([0, 1, 1, 2, 2, 3, 3, 0], dtype=np.uint32), 6 * (2 * 4))
        I2_aux += np.repeat(4 * np.arange(6, dtype=np.uint32), 8)


        if i==0:
            V=V_aux
            I1=I1_aux
            I2=I2_aux
        else:
            V=np.vstack((V,V_aux))
            I1=np.vstack((I1,I1_aux+i*24))
            I2=np.vstack((I2,I2_aux+i*24))




    return V, I1, I2


# -----------------------------------------------------------------------------
class Canvas(app.Canvas):

    def __init__(self):
        app.Canvas.__init__(self, keys='interactive', size=(800, 600))

        num_of_cubes=1 #number of cubes to draw
        self.V, self.filled, self.outline = cube(num_of_cubes)


        self.store_pos=np.array((0,0)) #for mouse interaction

        self.vert_data=np.vstack(self.V[:,0])
        self.V_buf=np.vstack(self.V[:,0])
        self.V_buf.dtype=[('a_position',np.float32,3)]
        self.vert_buf=gloo.VertexBuffer(self.V_buf)

        self.N_buf=np.vstack(self.V[:,1])
        self.N_buf.dtype=[('a_normal',np.float32,3)]
        self.norm_buf=gloo.VertexBuffer(self.N_buf)

        self.C_buf=np.vstack(self.V[:,2])
        self.C_buf.dtype=[('a_color',np.float32,4)]
        self.colo_buf=gloo.VertexBuffer(self.C_buf)

        self.filled_buf=gloo.IndexBuffer(self.filled.flatten())
        self.outline_buf=gloo.IndexBuffer(self.outline.flatten())

        self.program = gloo.Program(vert, frag)
        self.translate = 1

        #self.vert_buf=gloo.VertexBuffer(self.vertices.flatten())
        self.program.bind(self.vert_buf)
        self.program.bind(self.norm_buf)
        self.program.bind(self.colo_buf)


        self.view = translate((0, 0, -10))
        self.model = np.eye(4, dtype=np.float32)

        gloo.set_viewport(0, 0, self.physical_size[0], self.physical_size[1])
        self.projection = perspective(45.0, self.size[0] /
                                      float(self.size[1]), 2.0, 10.0)

        self.program['u_projection'] = self.projection

        self.program['u_model'] = self.model
        self.program['u_view'] = self.view

        self.theta = 0
        self.phi = 0

        gloo.set_clear_color('white')
        gloo.set_state('opaque')
        gloo.set_polygon_offset(1, 1)

        self._timer = app.Timer('auto', connect=self.on_timer, start=True)

        self.show()
        self.t=0
        
    def get_view_axis_in_scene_coordinates(view):
        import numpy
        tform=view.scene.transform
        w,h = view.canvas.size
        screen_center = numpy.array([w/2,h/2,0,1]) # in homogeneous screen coordinates
        d1 = numpy.array([0,0,1,0]) # in homogeneous screen coordinates
        point_in_front_of_screen_center = screen_center + d1 # in homogeneous screen coordinates
        p1 = tform.imap(point_in_front_of_screen_center) # in homogeneous scene coordinates
        p0 = tform.imap(screen_center) # in homogeneous screen coordinates
        assert(abs(p1[3]-1.0) < 1e-5) # normalization necessary before subtraction
        assert(abs(p0[3]-1.0) < 1e-5)
        return p0[0:3],p1[0:3] # 2 point representation of view axis in 3d scene coordinates

    # ---------------------------------
    def on_timer(self, event):
        self.update()

    # ---------------------------------
    def print_mouse_event(self, event, what):
        modifiers = ', '.join([key.name for key in event.modifiers])
        print('%s - pos: %r, button: %s, modifiers: %s, delta: %r' %
              (what, event.pos, event.button, modifiers, event.delta))

    def on_mouse_press(self, event):
        self.print_mouse_event(event, 'Mouse press')

        #convert to NDC
        left=event.pos[0]*2/self.size[0]-1
        bottom=(self.size[1]-event.pos[1])*2/self.size[1]-1
        
        tform= self.view.scene.transform
        w,h = self.view.canvas.size
        screen_center = event.pos # in homogeneous screen coordinates
        d1 = numpy.array([0,0,1,0]) # in homogeneous screen coordinates
        point_in_front_of_screen_center = screen_center + d1 # in homogeneous screen coordinates
        p1 = tform.imap(point_in_front_of_screen_center) # in homogeneous scene coordinates
        p0 = tform.imap(screen_center) # in homogeneous screen coordinates
        assert(abs(p1[3]-1.0) < 1e-5) # normalization necessary before subtraction
        assert(abs(p0[3]-1.0) < 1e-5)
        print(p0[0:3],p1[0:3]) # 2

        z_clip=np.linspace(-1.,1.,100)
        for val in z_clip:
            aux=np.dot(np.dot(np.linalg.inv(self.view),np.linalg.inv(self.projection)),np.array((left,bottom,val,1.)))
            pos3d=aux/aux[3]
            print(pos3d)

    def on_mouse_wheel(self, event):

        self.translate -= event.delta[1]
        self.translate = max(-1, self.translate)
        self.view[3,2]=-self.translate

        self.program['u_view'] = self.view
        self.update()



    def on_draw(self, event):
        gloo.clear()

        # Filled cube

        gloo.set_state(blend=False, depth_test=True, polygon_offset_fill=True)
        self.program['u_color'] = 1, 0, 1, 1


        self.program.draw('triangles', self.filled_buf)

        # Outline
        gloo.set_state(polygon_offset_fill=False, blend=True, depth_mask=False)
        gloo.set_depth_mask(False)
        self.program['u_color'] = 0, 0, 0, 1

        self.program.draw('lines', self.outline_buf)
        gloo.set_depth_mask(True)


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    c = Canvas()
    app.run()