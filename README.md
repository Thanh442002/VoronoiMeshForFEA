# Polygonal Mesh for FEM

## Voronoi diagram

Minh cho cho phần tử đa giác thường gặp trong tự nhiên 

![Natural Voronoi](fig/hexagonal.jpg)

Việc tạo các phần tử đa giác hay còn gọi là ô Voronoi dựa trên thư viện [Hull-Delaunay-Voronoi](https://github.com/Scrawk/Hull-Delaunay-Voronoi).

![voro](fig/voronoi.png) 

Trong đó hàm triển khai chính để lọc ra các đa giác trong miền bài toán từ các điểm gieo ngẫu nhiên là [GenerateVoronoi](https://github.com/Thanh442002/VoronoiMeshForFEA/blob/29f5e425a10c2549837c27a40b05c3c5baef7018/PolygonalMesher.cs#L1103). 

![noLloyd](fig/voronoi_noloop.png)

Thuật toán Lloyd's được áp dụng cải thiện lưới :thumbsup:

![lloyd](fig/voronoi_loop50.png)

Một vài hình ảnh lưới thực hiện được 

<img src="fig/mesh1700.png" width="425"/> <img src="fig/800mesh.png" width="250"/> 
