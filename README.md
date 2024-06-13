# Polygonal Mesh for FEM

## Voronoi diagram

Minh cho cho phần tử đa giác thường gặp trong tự nhiên 

![Natural Voronoi](fig/hexagonal.jpg)

Việc tạo các phần tử đa giác hay còn gọi là ô Voronoi dựa trên thư viện [Hull-Delaunay-Voronoi](https://github.com/Scrawk/Hull-Delaunay-Voronoi).

![voro](fig/voronoi.png)

Trong đó hàm triển khai chính để lọc ra các đa giác trong miền bài toán từ các điểm gieo ngẫu nhiên là [GenerateVoronoi](https://github.com/Thanh442002/VoronoiMeshForFEA/blob/29f5e425a10c2549837c27a40b05c3c5baef7018/PolygonalMesher.cs#L1103).

