import mitsuba
mitsuba.set_variant("scalar_rgb")
import materials
from mitsuba.core import ScalarTransform4f, Vector3f
import numpy as np

from mitsuba.core.xml import load_dict, load_string
import os

from mitsuba.core import Thread


##############
# INTEGRATOR #
##############
FNEE_integrator = {
    "type": "path_fnee",
    "max_depth": 10,
    "samples_per_pass": 1,
    "timeout": 2,
    "biased": False,
    "twostage": False,
    "halfvector_constraints": True,
    "mnee_init": False,
    "max_iterations": 20,
    "solver_threshold": 0.00001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
    "use_SMS": True,
    "crop_caustic": True,
}

SMS_integrator = {
    "type": "path_sms_ms",
    "max_depth": 8,
    "samples_per_pass": 1,
    "timeout": 2,
    "biased": False,
    "twostage": False,
    "halfvector_constraints": True,
    "mnee_init": False,
    "max_iterations": 20,
    "solver_threshold": 0.00001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
}

AOV_FNEE = {
    "type" : "aov",
    "aovs" : "N:sh_normal,Ng:geo_normal, P:position,Albedo:Albedo",
    #  "aovs" : "N:geo_normal",
    "timeout": 2,
    "samples_per_pass": 1,

    "integrator" : {
    "type": "path_fnee",
    "max_depth": 8,
    "biased": False,
    "twostage": False,
    "halfvector_constraints": True,
    "mnee_init": False,
    "max_iterations": 20,
    "solver_threshold": 0.00001,
    "uniqueness_threshold": 0.000001,
    "max_trials": 1000,
    "bounces": 2,
    "use_SMS": False
},
}

path_integrator = {
    "type": "path",
    "max_depth": 8,
    "samples_per_pass": 1,
    "timeout": 2,
}


################
# BSDF_DIFFUSE #
################
bsdf_diffuse = {
    "type": "diffuse",
    "reflectance": {
        "type": "rgb",
        "value": [0.4, 0.4, 0.4],
    }
}


Texture_manchon_1 = {
    "type": "bitmap",
    "id": "id_texture_manchon_1",
    "filename": "textures/Heightmaps/manchon_1_1024.exr",
    "raw": True
}

Texture_flat = {
    "type": "bitmap",
    "id": "id_texture_flat",
    "filename": "textures/Heightmaps/Plat.exr",
    "raw": True
}

Textures_heightfield = []
for i in range(3):
    t = {
        "type": "bitmap",
        "id": "id_texture_manchon_"+str(i),
        "filename": "textures/Heightmaps/manchon_"+str(i)+"_256.exr",
        "raw": True
    }
    Textures_heightfield.append(t)

panel_mesh_out = {"type": "obj",
                  "filename": "meshes/Vitres.obj",
                  "bsdf": materials.bsdf_window,
                  }

panel_mesh_in = {"type": "obj",
                 "filename": "meshes/Vitres_in.obj",
                 "bsdf": materials.bsdf_window,
                 }

##############
# SCENE_DICT #
##############
scene_dict = {
    "type": "scene",

   "Barreaux": {
       "type": "obj",
       "filename": "meshes/Barreaux.obj",
       "bsdf": materials.bsdf_barreaux,
       "caustic_receiver": True,
   },

    "Briques_interrieur": {
        "type": "obj",
        "filename": "meshes/Briques_interrieur.obj",
        "bsdf": materials.bsdf_briques_interrieur,
        "caustic_receiver": True,
    },
    
   "Moullures": {
       "type": "obj",
       "filename": "meshes/Moullures.obj",
       "bsdf": materials.bsdf_moulures,
       "caustic_receiver": True,
   },
    "Parquet": {
        "type": "obj",
        "filename": "meshes/Parquet.obj",
        "bsdf": materials.bsdf_parquet,
        "caustic_receiver": True,
    },
    "Parquet_haut": {
        "type": "obj",
        "filename": "meshes/Parquet_haut.obj",
        "bsdf": materials.bsdf_parquet_diffuse,
        "caustic_receiver": True,
    },
    "Poutres": {
        "type": "obj",
        "filename": "meshes/Poutres.obj",
        "bsdf": materials.bsdf_poutres,
        "caustic_receiver": True,
    },
    "Sous_basements": {
        "type": "obj",
        "filename": "meshes/Sous_bassements.obj",
        "bsdf": materials.bsdf_briques_interrieur,
        "caustic_receiver": True,
    },
     
     "light_1": {
         "type": "obj",
         "filename": "meshes/light.obj",
         "emitter": {
             "type": "area",
             "radiance": 800000.0,
             "caustic_emitter_multi": True,
             "caustic_emitter_single": True
         },        
     },
    
    "light_env": {
        "type": "envmap",
        "filename": "textures/env/sky.exr",
        "scale": 1.5,
        "caustic_emitter_multi": True,
        "caustic_emitter_single": True
    },

    "camera_1": {
        "type": "perspective",
        "fov": 40,
        "fov_axis": "smaller",
        "to_world": ScalarTransform4f.look_at(origin=[-20.206117630004883, 4.09607458114624, 0.9126691818237305],
                                              target=[-19.820100784301758, 4.041653156280518, 0.9240520000457764],
                                              up=[-0.0112719452008605, 0.0015847840113565326, 0.389833927154541]),
        "film": {
            "type": "hdrfilm",
            "rfilter": {"type": "box"},
            "width": 1024,
            "height": 1024,
        },
        "sampler": {
            "type": "independent",
            "sample_count": 100000,
        }
    },

    "integrator": path_integrator
}


def add_window_to_scene(scene_dict, window_name, id, to_world_1, to_world_2, H1, H2, texture_id_1, texture_id_2, L=2, W=2):
    Heighfield_dict = {
        "type": "heightfield",
        "to_world": ScalarTransform4f.translate([2.9, -3.26, 1.0]) * ScalarTransform4f.rotate([0, 0, 1], -90.0),
        "L": 1,
        "W": 1,
        "H": 0.00,
        "flip_normals": False,
        "caustic_bouncer": True,
        "caustic_caster_multi": False,
        "double_refraction_caster": True,
        "self_id": 1,
        "pair_id": 2,
        "partition_number": 1.0,
        "bsdf": materials.bsdf_window,
        "heightmap": Texture_manchon_1
    }

    Heighfield_dict["H"] = H1
    Heighfield_dict["to_world"] = to_world_1
    Heighfield_dict["flip_normals"] = False
    Heighfield_dict["self_id"] = id * 2
    Heighfield_dict["pair_id"] = id * 2 + 1
    Heighfield_dict["heightmap"] = texture_id_1
    Heighfield_dict["caustic_bouncer"] = True
    Heighfield_dict["caustic_caster_multi"] = False
    Heighfield_dict["double_refraction_caster"] = True
    Heighfield_dict["L"] = L
    Heighfield_dict["W"] = W
    scene_dict[window_name + "_1"] = Heighfield_dict.copy()

    Heighfield_dict["H"] = H2
    Heighfield_dict["self_id"] = id * 2 + 1
    Heighfield_dict["pair_id"] = id * 2
    Heighfield_dict["to_world"] = to_world_2
    Heighfield_dict["flip_normals"] = True
    Heighfield_dict["heightmap"] = texture_id_2
    Heighfield_dict["caustic_bouncer"] = False
    Heighfield_dict["caustic_caster_multi"] = True
    Heighfield_dict["double_refraction_caster"] = True
    Heighfield_dict["L"] = L
    Heighfield_dict["W"] = W

    scene_dict[window_name + "_2"] = Heighfield_dict.copy()

def build_heightfield_scene():
    scene_heightfield_dict = scene_dict.copy()
    Textures = [[Texture_manchon_1, Texture_flat]]
    e = 0.05
    H = [[-0.001, 0.00]]
    # H = [[0, 0.1], [0, 0.1], [0, 0.1], [0, 0.1]]
    index = 0

    origin = np.asarray([-18.4409 , 0.54889+e, 0.081797])
    backward = np.asarray([0,1,0])
    up = np.asarray([0,0,1])
    left = np.asarray([1,0,0])

    # 4x4 panels
    for i in range(0, 11):
        for j in range(3):
            index = (i+0)*3 + j
            point = origin + j*1.0*up + i*3.454*left 
            T1 = ScalarTransform4f.translate(
                point - e*backward) * ScalarTransform4f.rotate([0, 0, 1], -90.0)
            T2 = ScalarTransform4f.translate(
                point) * ScalarTransform4f.rotate([0, 0, 1], -90.0)
            w=1.0
            l = 1.28
            #if(j==2):
            #    w = 1.3
            add_window_to_scene(scene_heightfield_dict, "H"+str(index), index, T1, T2, H[0][0], H[0][1], Textures[0][0], Texture_manchon_1, L=l, W=w)
    return scene_heightfield_dict

