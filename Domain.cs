using DEMSoft.Drawing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DEMSoft.PolygonalMesher
{
    public class Domain
    {
        private double[] BdBox;
        private List<SubDomain> listSubDomain;

        public Domain(List<SubDomain> listSubDomain)
        {
            this.listSubDomain = listSubDomain;
            this.BdBox = ComputeBoundBox();
        }
        public double[] GetBoundBox()
        {
            return BdBox;
        }
        private double[] ComputeBoundBox()
        {
            double[] tempBox = new double[] { 999999, -999999, 999999, -999999 };//2D
            foreach (var subDomain in listSubDomain)
            {
                double[] bdBoxSubDomain = subDomain.GetBoundBox();
                tempBox[0] = Math.Min(tempBox[0], bdBoxSubDomain[0]);
                tempBox[1] = Math.Max(tempBox[1], bdBoxSubDomain[1]);
                tempBox[2] = Math.Min(tempBox[2], bdBoxSubDomain[2]);
                tempBox[3] = Math.Max(tempBox[3], bdBoxSubDomain[3]);
            }
            return tempBox;
        }

        public List<double[]> GenerateSeedRandom(int nElem, double dOffset = 0)
        {
            List<double[]> P = new List<double[]>();
            double[] Px = new double[nElem];
            double[] Py = new double[nElem];
            int Ctr = 0;
            Random rand = new Random();
            double[] Yx = new double[nElem];
            double[] Yy = new double[nElem];
            while (Ctr < nElem)
            {
                double[] bdBox = GetBoundBox();

                List<double[]> Y = new List<double[]>();
                for (int i = 0; i < nElem; i++)
                {
                    double[] xx = new double[2];

                    xx[0] = (bdBox[1] - bdBox[0]) * rand.NextDouble() + bdBox[0];
                    xx[1] = (bdBox[3] - bdBox[2]) * rand.NextDouble() + bdBox[2];
                    Y.Add(xx);
                }

                //double[] d = ComputeDistance(Y); //// normal problem
                //List<int> indexPointInside = FindIndexInside(d, dOffset);

                //////////////////specific problem
                List<double[]> dist = new List<double[]>();
                Dist d = new Dist(bdBox);
                dist = d.CantileverBeam(Y);
                //dist = d.Cilyndrical(Y);
                //dist = d.LShape(Y);
                //dist = d.Photoelastic_Samples_1(Y);
                //dist = d.Photoelastic_Samples_2(Y);
                List<int> indexPointInside = FindIndexInside(dist, dOffset);///Use for class Dist

                //////////////////

                int numAdded = Math.Min(nElem - Ctr, indexPointInside.Count);
                for (int i = 0; i < numAdded; i++)
                {
                    P.Add(Y[indexPointInside[i]]);
                }
                Ctr += numAdded;
            }
            return P;
        }

        public List<int> FindIndexInside(List<double[]> dist, double dOffset)
        {
            List<int> index = new List<int>();
            var end = dist[0].Length;
            for (int i = 0; i < dist.Count; i++)
            {
                if (dist[i][end - 1] + dOffset < 0)
                    index.Add(i);
            }

            return index;
        }

        public List<double[]> GenerateSeedHoneycomb(double size, double[] firstPoint, double dOffset = 0)
        {
            List<double[]> P = new List<double[]>();
            double[] bdBox = GetBoundBox();
            double xCurrent0 = firstPoint[0];
            double yCurrent0 = firstPoint[1];
            double xCurrent = firstPoint[0];
            double yCurrent = firstPoint[1];
            bool isPositive = true;
            while (yCurrent < bdBox[3])
            {
                while (xCurrent < bdBox[1])
                {
                    double[] currentPoint = new double[] { xCurrent, yCurrent };
                    if (isInside(currentPoint, dOffset))
                    {
                        P.Add(currentPoint);
                    }
                    xCurrent += size;
                }
                if (isPositive)
                {
                    xCurrent = xCurrent0 + size / 2.0;
                    isPositive = false;
                }
                else
                {
                    xCurrent = xCurrent0;
                    isPositive = true;
                }
                yCurrent += size * Math.Cos(Math.PI / 6.0);

            }
            return P;
        }
        #region single_edge
        //Lấy các điểm gần biên, tạo các hạt đối xứng tại 1 biên. Fail tại góc
        //public List<double[]> GenerateSeedsReflex(List<double[]> P, double alpha)
        //{
        //  List<double[]> nearBoundP = SelectSeedsNearBoundary(P, alpha);
        //  List<double[]> reflexP = new List<double[]>();
        //  double[] d = null;
        //  foreach (var sub in listSubDomain)
        //  {
        //    List<double[]> reflexPSub = sub.GenerateSeedsReflex(nearBoundP, out double[] dd);


        //    if (d == null)
        //    {
        //      d = dd;
        //      for (int i = 0; i < d.Length; i++)
        //      {
        //        reflexP.Add(reflexPSub[i]);
        //      }
        //    }
        //    else
        //    {
        //      for (int i = 0; i < d.Length; i++)
        //      {
        //        if (d[i] < dd[i])
        //        {
        //          d[i] = dd[i];
        //          reflexP.RemoveAt(i);
        //          reflexP.Insert(i, reflexPSub[i]);
        //        }
        //      }
        //    }
        //  }
        //  return reflexP;
        //}
        #endregion

        #region multiple_edge
        //public List<double[]> GenerateSeedsReflex(List<double[]> P, double alpha, int nElem, out List<double[]> nearBoundP)
        //{
        //    double eps = 1e-8;
        //    double eta = 0.9;
        //    List<double[]> reflexP = new List<double[]>();

        //    nearBoundP = new List<double[]>();

        //    foreach (var sub in listSubDomain)
        //    {
        //        List<double[]> subReflexP = sub.GenerateSeedsReflex(P, alpha);
        //        //nearBoundP.AddRange(nearBoundSubDomain);
        //        reflexP.AddRange(subReflexP);
        //    }

        //    return reflexP;
        //}
        ///////// Gốc///////
        public List<double[]> GenerateSeedsReflex(List<double[]> P, double alpha, int nElem, out List<double[]> nearBoundP)
        {
            //double eps = 1e-8;
            //double eta = 0.9;
            List<double[]> reflexP = new List<double[]>();
            nearBoundP = new List<double[]>();
            foreach (var sub in listSubDomain)
            {
                List<double[]> subReflexP = sub.GenerateSeedsReflex(P, alpha, out double[] d, out List<double[]> nearBoundSubDomain);
                nearBoundP.AddRange(nearBoundSubDomain);
                reflexP.AddRange(subReflexP);
            }

            double[] dd = ComputeDistance(reflexP);

            //List<int> indexInside = FindIndexInside(dd, 0);
            //for (int i = reflexP.Count - 1; i >= 0; i--)
            //{
            //    if (indexInside.Contains(i))
            //        reflexP.RemoveAt(i);
            //}
            return reflexP;
        }
        ////////////////////////////////////////////////////////////////////////////////////////
        //public List<double[]> GenerateSeedsReflex(List<double[]> P, double alpha)
        //{
        //  List<double[]> nearBoundP = SelectSeedsNearBoundary(P, alpha);
        //  List<double[]> reflexP = new List<double[]>();
        //  double[] d = null;
        //  foreach (var sub in listSubDomain)
        //  {
        //    List<double[]> reflexPSub = sub.GenerateSeedsReflex(nearBoundP, out double[] dd);

        //    if (d == null)
        //    {
        //      d = dd;
        //      for (int i = 0; i < d.Length; i++)
        //      {
        //        reflexP.Add(reflexPSub[i]);
        //      }
        //    }
        //    else
        //    {
        //      for (int i = 0; i < d.Length; i++)
        //      {
        //        if (d[i] < dd[i])
        //        {
        //          d[i] = dd[i];
        //          reflexP.RemoveAt(i);
        //          reflexP.Insert(i, reflexPSub[i]);
        //        }
        //      }
        //    }
        //  }
        //  return reflexP;
        //}
        #endregion


        public double[] ComputeDistance(List<double[]> P)
        {
            double[] d = null;
            foreach (var sub in listSubDomain)
            {
                double[] dSub = sub.ComputeDistance(P);
                if (d == null)
                    d = dSub;
                else
                {
                    for (int i = 0; i < d.Length; i++)
                    {
                        d[i] = Math.Max(d[i], dSub[i]);
                    }
                }
            }
            return d;
        }
        /////////////////////////////////

        public List<double[]> SelectPointsInside(List<double[]> P, double offset)
        {
            double[] d = ComputeDistance(P);
            
            List<double[]> result = new List<double[]>();
            for (int i = 0; i < P.Count; i++)
            {
                if (d[i] + offset < 0)
                {
                    result.Add(P[i]);
                }
            }
            return result;
        }
        public bool isInside(double[] p, double offset)
        {
            List<double[]> PP = new List<double[]>();
            PP.Add(p);
            if (!((p[0] < BdBox[0]) || (p[0] > BdBox[1]) || (p[1] < BdBox[2]) || (p[1] > BdBox[3])))
            {
                double[] d = ComputeDistance(PP);////

                return (d[0] + offset) <= 0 ? true : false;
            }
            return false;
        }
        public bool isNearBound(double[] p, double offset)//11/9
        {
            List<double[]> PP = new List<double[]>();
            PP.Add(p);
            if (!((p[0] < BdBox[0]) || (p[0] > BdBox[1]) || (p[1] < BdBox[2]) || (p[1] > BdBox[3])))
            {
                double[] d = ComputeDistance(PP);

                return (Math.Round(d[0] + offset, 3) <= 0.0) ? true : false;
                //return (d[0] + offset) <= 0.0 ? true : false;
            }

            return false;

        }
        public bool isInsideBoundBox(double[] p, double offset)
        {
            double dx = BdBox[1] - BdBox[0];
            double dy = BdBox[3] - BdBox[2];

            return !((p[0] < BdBox[0] - offset * dx) || (p[0] > BdBox[1] + offset * dx) || (p[1] < BdBox[2] - offset * dy) || (p[1] > BdBox[3] + offset * dy));
        }
        private List<int> FindIndexInside(double[] d, double dOffset)
        {
            List<int> index = new List<int>();
            for (int i = 0; i < d.Length; i++)
            {
                if (d[i] + dOffset < 0)
                    index.Add(i);
            }
            return index;
        }
        private List<int> FindIndexOutside(double[] d, double dOffset)
        {
            List<int> index = new List<int>();
            for (int i = 0; i < d.Length; i++)
            {
                if (d[i] + dOffset > 0)
                    index.Add(i);
            }
            return index;
        }

    }
}
