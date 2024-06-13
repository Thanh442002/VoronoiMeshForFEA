# Polygonal Mesh for FEM

## Voronoi diagram

Illustrate polygonal elements commonly found in nature

![Natural Voronoi](fig/hexagonal.jpg)

The creation of polygonal elements or Voronoi cells is based on the library [Hull-Delaunay-Voronoi](https://github.com/Scrawk/Hull-Delaunay-Voronoi).

![voro](fig/voronoi.png) 

In which the main implementation function to filter out polygons in the problem domain from randomly planted points is [GenerateVoronoi](https://github.com/Thanh442002/VoronoiMeshForFEA/blob/29f5e425a10c2549837c27a40b05c3c5baef7018/PolygonalMesher.cs#L1103). 

![noLloyd](fig/voronoi_noloop.png)

Lloyd's algorithm is applied to improve the mesh :thumbsup:

![lloyd](fig/voronoi_loop50.png)

A few mesh images

<img src="fig/mesh1700.png" width="425"/> <img src="fig/800mesh.png" width="250"/> 

<img src="fig/1.png" width="425"/> <img src="fig/2.png" width="250"/> 
