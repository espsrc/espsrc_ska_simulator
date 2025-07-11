import sys
from PIL import Image, ImageChops
import numpy as np
def show_exc(e):
    exc_type, exc_obj, exc_tb = sys.exc_info()
    return ("ERROR ===:> [%s in %s:%d]: %s" % (exc_type, exc_tb.tb_frame.f_code.co_filename, exc_tb.tb_lineno, str(e)))

def show_warn(e):
    exc_type, exc_obj, exc_tb = sys.exc_info()
    return ("WARNING =:> [%s in %s:%d]: %s" % (exc_type, exc_tb.tb_frame.f_code.co_filename, exc_tb.tb_lineno, str(e)))

sign = lambda a:int(a>0)-int(a<0)

def pngHtml(path=None):
    if not path:
        import os
        path = os.path.dirname(os.path.realpath(__file__))
    print (path)

def autocrop(path, newpath=None, fmt="PNG"):
    image = Image.open(path)
    (x,y) = (image.width, image.height)
    matrix = np.array(image.getdata())
    matrix = matrix / 255
    matrix = np.array([ item[0] * item[1] * item[2] * item[3] for item in matrix[:] ])
    matrix = matrix.reshape((y,x))
    rows = []
    cols = []

    for i in range(y):
        if (not (np.all(matrix[i,:] == 1.))):
            rows.append(i)
    for j in range(x):
        if (not (np.all(matrix[:,j] == 1.))):
            cols.append(j)
    aux = image.crop((min(cols), min(rows), max(cols), max(rows)))
    if newpath is None:
        aux.save(path, fmt)
    else:
        aux.save(newpath, fmt)
    return (aux)

