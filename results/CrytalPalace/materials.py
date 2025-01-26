
#---------------------------------
#           textures
#---------------------------------

tex_Chanvre_Briques_Inst_D= {
        "type": "bitmap",
        "id": "Chanvre_Briques_Inst_D",
        "filename": "textures/Briques_interrieur/Briques_Palace-SousBassement_D.bmp"
        }

tex_Chanvre_Briques_Inst_N= {
        "type": "bitmap",
        "id": "Chanvre_Briques_Inst_N",
        "filename": "textures/Briques_interrieur/Briques_Palace-SousBassement_N.bmp"
        }

tex_Mi_Parquet_D= {
        "type": "bitmap",
        "id": "Mi_Parquet_D",
        "filename": "textures/Parquet/Mi_Parquet_D.bmp"
        }

tex_Mi_Poutre_D= {
        "type": "bitmap",
        "id": "Mi_Poutre_D",
        "filename": "textures/Poutres/Mi_PoutresBois_D.bmp"
        }

#---------------------------------
#             BSDF
#---------------------------------
bsdf_window = {
    "type": "dielectric",
    "int_ior": 1.54,
    "id": "Mi_bsdf_window",
    }


#bsdf_barreaux = {
#    "type": "roughconductor",
#    "material": "Au",
#    "alpha": 0.45,
#    }

bsdf_barreaux = {
    "type": "diffuse",
    }

bsdf_briques_interrieur = {
    "type": "diffuse",
    "reflectance": tex_Chanvre_Briques_Inst_D,
    "id": "Mi_bsdf_window",
    }

# bsdf_briques_interrieur = {
#         "type" : "normalmap",
#         "filename" : "textures/Briques_interrieur/Briques_Palace-SousBassement_N.exr",
#         "bsdf" : bsdf_briques_interrieur_D,
# }

# <bsdf type="normalmap">
#     <texture name="normalmap" type="bitmap">
#         <boolean name="raw" value="true"/>
#         <string name="filename" value="textures/normalmap.jpg"/>
#     </texture>
#     <bsdf type="roughplastic"/>
# </bsdf>



bsdf_moulures = {
        "type": "diffuse",
        }


bsdf_parquet = {
        "type": "diffuse",
        "reflectance": tex_Mi_Parquet_D
        }

bsdf_parquet_diffuse = {
        "type": "diffuse",
        "reflectance": tex_Mi_Parquet_D
        }

bsdf_poutres= {
        "type": "diffuse",
        "reflectance": tex_Mi_Poutre_D,
        }


