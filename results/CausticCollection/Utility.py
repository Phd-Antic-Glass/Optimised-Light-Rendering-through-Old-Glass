
# FROM https://gist.github.com/arseniy-panfilov/4dc8fc5131277affe64619b1a9d00da0
# 
from PIL import Image
import OpenEXR
import Imath
import numpy
#  import numexpr as ne


FLOAT = Imath.PixelType(Imath.PixelType.FLOAT)

def exr_to_array(exrfile):
    file = OpenEXR.InputFile(exrfile)
    dw = file.header()['dataWindow']

    channels = list(file.header()['channels'].keys())
    channels_list = [c for c in ('R', 'G', 'B', 'A') if c in channels]
    size = (dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1)

    color_channels = file.channels(channels_list, FLOAT)
    channels_tuple = [numpy.frombuffer(channel, dtype='f') for channel in color_channels]
    
    return numpy.dstack(channels_tuple).reshape(size + (len(channels_tuple),))

  
def render(scene, name, jpg=False):
    # Get the scene's sensor (if many, can pick one by specifying the index)
    sensor = scene.sensors()[0]

    # Call the scene's integrator to render the loaded scene with the desired sensor
    scene.integrator().render(scene, sensor)

    # The rendered data is stored in the film
    film = sensor.film()

    # Write out data as high dynamic range OpenEXR file
    film.set_destination_file(name)
    film.develop()

    if(jpg):
        # Write out a tone-mapped JPG of the same rendering
        from mitsuba.core import Bitmap, Struct
        img = film.bitmap(raw=True).convert(
            Bitmap.PixelFormat.RGB, Struct.Type.UInt8, srgb_gamma=True)
        img.write(name + ".jpg")
    
    
def exr_to_srgb(exrfile):
    array = exr_to_array(exrfile)
    
    return array
