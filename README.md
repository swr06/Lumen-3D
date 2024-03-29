# Lumen Engine

(For the record, I came up with this name much much before Unreal did :P) 

An ***Unfinished*** **TOY ENGINE**. This project is mostly used as a testing ground to experiment with and learn various graphics techniques. 

## Current Feature List 

- Model Loading 
- Voxelization (with support for cascades)
- Indirect Multibounce Diffuse Lighting (Voxel Raytracing + Screenspace RTGI/SSGI) 
- Indirect Specular/Reflections (Voxel Raytracing + (parallax correct) Cubemap/Screenspace reflections)
- Support for checkerboard rendering (Heavily used for indirect lighting and other raytraced effects)
- Temporal filtering/supersampling 
- SVGF For Indirect Diffuse and specialized denoising filter for indirect specular lighting
- Cook torrance BRDF for direct lighting
- Shadowmapping (with PCF filtering)
- TAA + TAA-Upscaling
- FXAA 3.11 
- Volumetric Fog (Crepuscular rays + exponential fog)
- Post Processing pipeline (Bokeh DOF, Bloom, Chromatic Aberration, Contrast Adaptive Sharpening, Film Grain)
- Skybox/Environment Map Support

## Known issues
- Shadow mapping iffy at steep-ish sun angles 
- Noise on aggressive movement or disocclusions 

## Planned 
- Physically based sky rendering (with volumetric clouds)
- Cascaded shadow maps with PCSS filtering

## Special Thanks
- [UglySwedishFish](https://github.com/UglySwedishFish)
- [Lars](https://github.com/Ciwiel3/)

## License
- MIT (See LICENSE)

## Notice
This project is purely educational. I own none of the assets. All the rights go to their respective owners.

## Screenshots 

</br>

![s1](https://github.com/swr06/Lumen-Engine/blob/main/Screenshots/1.png)

</br>

</br>

![s2](https://github.com/swr06/Lumen-Engine/blob/main/Screenshots/2.jpg)

</br>

</br>

![s3](https://github.com/swr06/Lumen-Engine/blob/main/Screenshots/3.jpg)

</br>

</br>

![s4](https://github.com/swr06/Lumen-Engine/blob/main/Screenshots/4.png)

</br>

</br>

![s5](https://github.com/swr06/Lumen-Engine/blob/main/Screenshots/5.jpg)

</br>

</br>

![s6](https://github.com/swr06/Lumen-Engine/blob/main/Screenshots/6.jpg)

</br>

# Supporting

If you'd like to support this project, consider starring it. Thanks!
