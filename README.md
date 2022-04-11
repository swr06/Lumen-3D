# Lumen Engine

(For the record, I came up with this name much much before unreal did.) 


A **TOY ENGINE**. This project is mostly used as a testing ground to experiment with various graphics techniques. 

## Current Feature List 

- Voxelization
- Indirect Diffuse Lighting (Voxel Raytracing) 
- Indirect Specular/Reflections (Voxel Raytracing + (parallax correct) Cubemap/Screenspace reflections)
- Cook torrance BRDF for direct lighting
- Shadowmapping (with PCF filtering)
- Spatio-Temporal Variance Guided Filtering for indirect lighting
- TAA, FXAA 
- Volumetric Lighting
- Bloom

## Planned 

### Lighting 

- PCSS, CSM, Omnidirectional contact-hardening shadows
- Mesh tracing using SAH BVH (with surfel GI)

### Rendering
- Transparency/Translucency OIT

### Post process
- DOF, Auto exposure

### Volumetrics 
- Clouds

### Other 
- Physically based atmosphere rendering
