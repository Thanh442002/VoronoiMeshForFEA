using DEMSoft.Drawing.Geometry;
using DEMSoft.NURBS;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace DEMSoft.PolygonalMesher
{
    public class SubDomain
    {
        private double[] BdBox;
        private List<Abstract1DParametricGeometry> bound;
        private List<Abstract1DParametricGeometry> boundReflex;
        private bool isAdd;

        public SubDomain(List<Abstract1DParametricGeometry> curves, List<Abstract1DParametricGeometry> curvesReflex, bool isAdd = true)
        {
            this.bound = curves;
            this.boundReflex = curvesReflex;
            BdBox = ComputeBoundBox();
            this.isAdd = isAdd;
        }
        public bool GetIsAdd()
        { return isAdd; }

        public double[] GetBoundBox()
        {
            return BdBox;
        }
        private double[] ComputeBoundBox()
        {
            double[] tempBox = new double[] { 999999, -999999, 999999, -999999 };//2D
            foreach (var curve in bound)
            {
                double[] bdBoxSubDomain = curve.GetBoundBox();
                tempBox[0] = Math.Min(tempBox[0], bdBoxSubDomain[0]);
                tempBox[1] = Math.Max(tempBox[1], bdBoxSubDomain[1]);
                tempBox[2] = Math.Min(tempBox[2], bdBoxSubDomain[2]);
                tempBox[3] = Math.Max(tempBox[3], bdBoxSubDomain[3]);
            }
            return tempBox;
        }
        public double[] ComputeDistance(List<double[]> P)
        {
            double[] d = new double[P.Count];
            for (int i = 0; i < P.Count; i++)
            {
                d[i] = ComputeDistance(P[i]);
            }
            return d;
        }


        public double ComputeDistance(double[] p)
        {
            double minDistance = -999999999;
            for (int i = 0; i < bound.Count; i++)
            {
                var control = bound[i].SelectAllControlPoints();
                double xiProjection = bound[i].Projection(p[0], p[1], 0)[0];
                double[] aProjection = bound[i].TangentUnitVectorCurve(xiProjection); //Vector tiep tuyen don vi cua duong cong
                double[] pointProjection = bound[i].PointAt(xiProjection);
                //double d = ((pointProjection[0] - p[0]) * aProjection[1] - (pointProjection[1] - p[1]) * aProjection[0]) *
                //    Math.Sqrt(Math.Pow(pointProjection[0] - p[0], 2) + Math.Pow(pointProjection[1] - p[1], 2));
                double d = ((pointProjection[0] - p[0]) * aProjection[1] - (pointProjection[1] - p[1]) * aProjection[0]) / Math.Sqrt(aProjection[0] * aProjection[0] + aProjection[1] * aProjection[1]);

                minDistance = Math.Max(minDistance, d);
            }
            return minDistance;
        }
        public double[][] ComputeDistanceToReflex(List<double[]> P)
        {
            double[][] d = new double[P.Count][];
            for (int i = 0; i < P.Count; i++)
            {
                d[i] = ComputeDistanceToReflex(P[i]);
            }
            return d;
        }
        public double[] ComputeDistanceToReflex(double[] p)
        {
            double minDistance = -999999999;
            double[] dis = new double[boundReflex.Count + 1];
            for (int i = 0; i < boundReflex.Count; i++)
            {
                double xiProjection = boundReflex[i].Projection(p[0], p[1], 0)[0];
                double[] aProjection = boundReflex[i].TangentUnitVectorCurve(xiProjection); //Vector tiep tuyen don vi cua duong cong
                double[] pointProjection = boundReflex[i].PointAt(xiProjection);
                double d = ((pointProjection[0] - p[0]) * aProjection[1] - (pointProjection[1] - p[1]) * aProjection[0]) / Math.Sqrt(aProjection[0] * aProjection[0] + aProjection[1] * aProjection[1]);
                minDistance = Math.Max(minDistance, d);

                dis[i] = d;

            }
            dis[boundReflex.Count] = minDistance;
            return dis;
        }

        #region multiple_edge
        public List<double[]> GenerateSingleSeedReflex(double[] p, double alpha, out double minDistance)//28_9_2023, 
        {
            minDistance = -999999999;
            List<double[]> reflexP = new List<double[]>();

            for (int i = 0; i < boundReflex.Count; i++)
            {
                double xiProjection = boundReflex[i].Projection(p[0], p[1], 0)[0];
                double[] aProjection = boundReflex[i].TangentUnitVectorCurve(xiProjection);
                double[] pointProjection = boundReflex[i].PointAt(xiProjection);
                //if (xiProjection == 0 || xiProjection == 1)
                //{
                double d = ((pointProjection[0] - p[0]) * aProjection[1] - (pointProjection[1] - p[1]) * aProjection[0]) / Math.Sqrt(aProjection[0] * aProjection[0] + aProjection[1] * aProjection[1]);

                double x = 1.0 / (aProjection[0] * aProjection[0] + aProjection[1] * aProjection[1]) * ((aProjection[0] * aProjection[0] - aProjection[1] * aProjection[1]) * p[0] + 2 * aProjection[1] * aProjection[1] * pointProjection[0] + 2 * aProjection[0] * aProjection[1] * (p[1] - pointProjection[1]));
                double y = 0;
                if (aProjection[0] != 0)
                {
                    y = aProjection[1] / aProjection[0] * (x + p[0] - 2 * pointProjection[0]) - p[1] + 2 * pointProjection[1];
                }
                else
                {
                    y = p[1];
                }

                if (d <= 0)
                {
                    if (minDistance < d)
                    {
                        minDistance = d;
                    }
                    reflexP.Add(new double[] { x, y });
                }
            }
            return reflexP;
        }

        public List<double[]> SelectSeedsNearBoundary(List<double[]> P, double alpha)
        {
            List<double[]> nearBoundP = new List<double[]>();
            double[] d = ComputeDistance(P);
            List<int> indexPointNearBoundary = FindIndexNearBoundary(d, alpha);
            for (int i = 0; i < indexPointNearBoundary.Count; i++)
            {
                nearBoundP.Add(P[indexPointNearBoundary[i]]);
            }
            return nearBoundP;
        }

        private List<int> FindIndexNearBoundary(double[] d, double alpha)
        {
            List<int> index = new List<int>();
            for (int i = 0; i < d.Length; i++)
            {
                if (Math.Abs(d[i]) < alpha)
                    index.Add(i);
            }
            return index;
        }
        /////////////Gốc////////
        public List<double[]> GenerateSeedsReflex(List<double[]> P, double alpha, out double[] d, out List<double[]> nearBoundSubDomain)
        {
            double eta = 0.9;
            List<double[]> reflexP = new List<double[]>();
            nearBoundSubDomain = SelectSeedsNearBoundary(P, alpha);///////////
            d = new double[nearBoundSubDomain.Count];
            double[] distanceToReflex = ComputeDistance(nearBoundSubDomain);
            for (int i = 0; i < nearBoundSubDomain.Count; i++)
            {
                //if (Math.Abs(distanceToReflex[i]) < alpha)
                //{
                List<double[]> singleSeedReflex = GenerateSingleSeedReflex(nearBoundSubDomain[i], alpha, out double dSingle);///15/1/2024
                double[] distance = ComputeDistance(singleSeedReflex);
                for (int j = 0; j < singleSeedReflex.Count(); j++)
                {
                    if (Math.Abs(distance[j]) <= alpha && distance[j] > 0)
                    {
                        reflexP.Add(singleSeedReflex[j]);
                    }
                }
                //reflexP.AddRange(singleSeedReflex);
                d[i] = dSingle;
                //}
            }
            return reflexP;
        }
        //////////////////////////////
        public List<double[]> GenerateSeedsReflex(List<double[]> P, double alpha)
        {
            double eta = 0.9;
            double eps = 1e-8;
            var dis = ComputeDistanceToReflex(P);
            
            var NBdrySegs = dis[1].Length - 1;
            List<double[]> reflexP = new List<double[]>();

            var sub_dis_n1 = ComputeDistanceToReflex(P.Select(s => new double[] { s[0] + eps, s[1] }).ToList());
            var sub_dis_n2 = ComputeDistanceToReflex(P.Select(s => new double[] { s[0], s[1] + eps }).ToList());

            var n1 = SubtractArrays(sub_dis_n1, dis).Select(row => row.Select(elem => elem / eps).ToArray()).ToArray();
            var n2 = SubtractArrays(sub_dis_n2, dis).Select(row => row.Select(elem => elem / eps).ToArray()).ToArray();

            var nearBound = SelectSeedsNearBoundary(P, alpha);
            double[,] I = new double[P.Count, NBdrySegs];
            List<double> disI = new List<double>();
            for (int i = 0; i < P.Count; i++)
            {
                for (int j = 0; j < NBdrySegs; j++)
                {
                    if (Math.Abs(dis[i][j]) < alpha)
                    {
                        double[] rp = new double[2];
                        I[i, j] = 1;
                        disI.Add(dis[i][j]);
                        rp[0] = P[i][0] - 2 * n1[i][j] * sub_dis_n1[i][j];
                        rp[1] = P[i][1] - 2 * n2[i][j] * sub_dis_n2[i][j];
                        reflexP.Add(rp);
                    }
                }
            }
            var d_R_P = ComputeDistanceToReflex(reflexP);

            List<double[]> R_P = new List<double[]>();

            for (int i = 0; i < reflexP.Count; i++)
            {
                if (Math.Abs(d_R_P[i][NBdrySegs]) >= eta * Math.Abs(disI[i]) && d_R_P[i][NBdrySegs] > 0)
                {
                    R_P.Add(reflexP[i]);
                }
            }

            ////////////////
            R_P = GetUniqueRows(R_P).ToList();
            return R_P;
        }
        #endregion

        #region single_edge
        //Lấy các điểm gần biên, tạo các hạt đối xứng tại 1 biên. Fail tại góc
        //public double[] GenerateSingleSeedReflex(double[] p, out double minDistance)
        //{
        //  minDistance = -999999999;
        //  double[] reflexP = new double[2];
        //  for (int i = 0; i < bound.Count; i++)
        //  {
        //    double xiProjection = bound[i].Projection(p[0], p[1], 0)[0];
        //    double[] pointProjection = bound[i].PointAt(xiProjection);
        //    double[] aProjection = bound[i].TangentUnitVectorCurve(xiProjection);
        //    double d = ((pointProjection[0] - p[0]) * aProjection[1] - (pointProjection[1] - p[1]) * aProjection[0]) * Math.Sqrt(Math.Pow(pointProjection[0] - p[0], 2) + Math.Pow(pointProjection[1] - p[1], 2));
        //    if (minDistance < d)
        //    {
        //      minDistance = d;//tính có dấu nên lấy số dương lớn nhất
        //      reflexP[0] = 2.0 * pointProjection[0] - p[0];
        //      reflexP[1] = 2.0 * pointProjection[1] - p[1];
        //    }
        //  }
        //  return reflexP;
        //}

        //public List<double[]> GenerateSeedsReflex(List<double[]> P, out double[] d)
        //{
        //  List<double[]> reflexP = new List<double[]>();
        //  d = new double[P.Count];
        //  for (int i = 0; i < P.Count; i++)
        //  {
        //    reflexP.Add(GenerateSingleSeedReflex(P[i], out double dSingle));
        //    d[i] = dSingle;
        //  }
        //  return reflexP;
        //}
        #endregion
        public static double[][] SubtractArrays(double[][] array1, double[][] array2)
        {
            if (array1.Length != array2.Length || array1[0].Length != array2[0].Length)
            {
                throw new ArgumentException("Arrays must have the same dimensions.");
            }

            return array1.Zip(array2, (row1, row2) => row1.Zip(row2, (elem1, elem2) => elem1 - elem2).ToArray()).ToArray();
        }
         public static IEnumerable<double[]> GetUniqueRows(List<double[]> inputList)
        {
            HashSet<double> seenFirstElements = new HashSet<double>();

            foreach (var row in inputList)
            {
                if (seenFirstElements.Add(row[0]))
                {
                    yield return row;
                }
            }
        }
    }
}
