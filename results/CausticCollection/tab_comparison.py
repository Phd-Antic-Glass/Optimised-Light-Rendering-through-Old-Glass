##################
# MITSUBA IMPORT #
##################
import mitsuba

import os
import sys
import Caustic_collection_scene as s



from mitsuba.core import Thread
from mitsuba.core.xml import load_file
from mitsuba.core.xml import load_dict

from mitsuba.render import FermatNEE
from mitsuba.render import SMSConfig, SpecularManifold, SurfaceInteraction3f
from mitsuba.core import Vector3f, Ray3f

import numpy as np
from matplotlib import pyplot as plt
from Utility import *
import math


##########
# RENDER #
##########
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


#################
# SCENE LOADING #
#################
resolution = 1024
spp = 10000

# Absolute or relative path to the scene XML file
filename = './Caustic_collection.xml'

# Add the scene directory to the FileResolver's search path
Thread.thread().file_resolver().append(os.path.dirname(filename))
############
# # RENDER #
############
scene_dict = s.build_heightfield_scene()
scene = load_dict(scene_dict)
# output_folder = "CausticCollection/"


def gather(mf):
    res=256
    samples=100
    cond = []
    failures = 0
    success_c = 0
    epsilon=1e-6
    solution_number=0
    total_launch=0
    for i in range(0,res):
        print("gather: ",i,"/",res,end="\r")
        for j in range(0,res):
            solution_list=[]
            for r in range(samples):
                total_launch = total_launch+1
                ray, _ = camera.sample_ray_differential(0.0, sampler.next_1d(), [float(i)/float(res), float(j)/float(res)], [0.5, 0.5])
                si = scene.ray_intersect(ray)
                if(si.is_valid() and si.shape.is_caustic_receiver()):
                    ei = SpecularManifold.sample_emitter_interaction(si, scene.emitters(), sampler)
                    output=0.0
                    success, direction = mf.specular_connection_debug(si, ei, [-1,-1,-1,-1])
                    if(success==True):
                        success_c=success_c+1
                        duplicate = False
                        for s in solution_list:
                            if np.abs((direction @ s) - 1.0) < epsilon:
                                duplicate = True
                                break
                        if not duplicate:
                            solution_list.append(direction)

                        cond.append(output)
                    else:
                        failures=failures+1

            solution_number = solution_number + len(solution_list)
    print("\nDONE")
    total_iterations = res*res*samples
    ratio_success = (float(success_c)/float(total_launch))

    return np.asarray(cond), ratio_success, solution_number


scene_dict["camera_1"]["sampler"]["seed"] = 2
sampler = load_dict(scene_dict["camera_1"]["sampler"])
camera = load_dict(scene_dict["camera_1"])

sms_config = SMSConfig()
sms_config.max_iterations = 20
sms_config.solver_threshold = 1e-4
sms_config.uniqueness_threshold = 1e-5
sms_config.bounces = 2;
sms_config.max_trials = 1000

# h_i evaluations, in millions
def nb_eval(Heightfields):
    c = 0.0
    for h in Heightfields:
        c = c + h.eval_attribute_1("nbEval", SurfaceInteraction3f(), True)*1e-6
    return c

## fermat solver
mf = FermatNEE()
mf.init(scene, sampler, 1e-5, 1000, 0.5, 0.5, True, sms_config, False)
cond_fnee, success_FNEE, unique_solution_number_FNEE = gather(mf)
FNEE_eval = nb_eval(scene.caustic_casters_double_refraction())

# sms solver
mf = FermatNEE()
mf.init(scene, sampler, 1e-5, 1000, 0.5, 0.5, True, sms_config, True)
cond_sms, success_SMS, unique_solution_number_SMS = gather(mf)
SMS_eval = nb_eval(scene.caustic_casters_double_refraction())

print("\nFNEE:")
print("success rate ",success_FNEE*100, "%")
print("unique solutions discovered ",unique_solution_number_FNEE)
print("nb eval: ", FNEE_eval)
print("\nSMS:")
print("success rate ",success_SMS*100, "%")
print("unique solutions discovered ",unique_solution_number_SMS)
print("nb eval: ", SMS_eval)



# Evaluation number $h_{eval}$
# # scene_dict = s.build_mesh_scene()
# scene_dict["camera_1"]["sampler"]["sample_count"] = 100
# scene = load_dict(scene_dict)
# s.FNEE_integrator["timeout"] = float(-1)
# scene_dict["integrator"] = load_dict(s.FNEE_integrator)
# scene = load_dict(scene_dict)
# render(scene, "test.exr")

# scene_dict["camera_1"]["sampler"]["sample_count"] = 100
# scene = load_dict(scene_dict)
# s.FNEE_integrator_SMS["timeout"] = float(-1)
# scene_dict["integrator"] = load_dict(s.FNEE_integrator_SMS)
# scene = load_dict(scene_dict)
# render(scene, "test.exr")

# # ------------------------------------------------------------------- #
# s.FNEE_integrator["timeout"] = float(timeout[i])
# scene_dict["integrator"] = load_dict(s.FNEE_integrator)
# scene = load_dict(scene_dict)

# s.FNEE_integrator["timeout"] = float(timeout[i])
# scene_dict["integrator"] = load_dict(s.FNEE_integrator)
# scene = load_dict(scene_dict)



