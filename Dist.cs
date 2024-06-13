using DEMSoft.Drawing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DEMSoft.PolygonalMesher
{
    public class Dist
    {
        private List<double[]> P;
        private Domain domain;
        private double[] bdBox;
        public Dist(double[] Bdbox)
        {
            this.bdBox = Bdbox;
        }
        ////////Photoelastic_Samples_1/////////

        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////

        public List<double[]> Photoelastic_Samples_1(List<double[]> P)
        {

            List<double[]> dist = new List<double[]>();

            var d1 = dRectangle(P, bdBox[0], bdBox[1], bdBox[2], bdBox[3]);
            var d2 = dRectangle(P, 15, 45, 0, 35);
            dist = dDiff(d1, d2);

            return dist;
        }
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////////////////////////////////////
        ///////// /arc           arc/ /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////
        /////////                     /////////

        public List<double[]> Photoelastic_Samples_2(List<double[]> P)
        {
            List<double[]> dist = new List<double[]>();
            var d1 = dRectangle(P, bdBox[0], bdBox[1], bdBox[2], bdBox[3]);

            var d2 = dCircle(P, 22.5, 27.5, 7.5);
            var d3 = dCircle(P, 37.5, 27.5, 7.5);
            var d4 = dRectangle(P, 15, 45, 0, 27.5);
            var d5 = dRectangle(P, 22.5, 37.5, 27.5, 35);

            dist = dDiff(d1, dUnion(d4, (dUnion(d3, dUnion(d2, d5)))));
            return dist;

        }
        public List<double[]> Cilyndrical(List<double[]> P)
        {
            List<double[]> dist = new List<double[]>();
            var d1 = dCircle(P, 0, 0, 2);
            var d2 = dCircle(P, 0, 0, 1);
            var d3 = dRectangle(P, 0, 2, 0, 2);
            dist = dDiff(d3, dUnion (d1, d2));
            return dist;
        }

        public List<double[]> CantileverBeam(List<double[]> P)
        {
            List<double[]> dist = new List<double[]>();
            //dist = dRectangle(P, 0, 3, 0, 1);
            dist = dRectangle(P, 0, 3, -0.5, 0.5);
            return dist;
        }

        public List<double[]> LShape(List<double[]> P)
        {
            List<double[]> dist = new List<double[]>();
            var d1 = dRectangle(P, bdBox[0], bdBox[1], bdBox[2], bdBox[3]);
            var d2 = dRectangle(P, 12, 31, 12, 31);
            dist = dDiff(d1, d2);
            return dist;

        }
        public List<double[]> dRectangle(List<double[]> P, double x1, double x2, double y1, double y2)
        {
            List<double[]> d = new List<double[]>();

            foreach (var point in P)
            {
                double[] row = new double[5];
                row[0] = x1 - point[0];
                row[1] = point[0] - x2;
                row[2] = y1 - point[1];
                row[3] = point[1] - y2;

                // Calculate the maximum value in the row
                row[4] = row.Take(4).Max();

                d.Add(row);
            }

            return d;
        }
        public List<double[]> dCircle(List<double[]> P, double xc, double yc, double r)
        {
            List<double[]> d = new List<double[]>();

            foreach (var point in P)
            {
                double[] row = new double[2];
                double distance = Math.Sqrt(Math.Pow(point[0] - xc, 2) + Math.Pow(point[1] - yc, 2)) - r;
                row[0] = distance;
                row[1] = distance;
                d.Add(row);
            }
            return d;
        }
        public List<double[]> dDiff(List<double[]> d1, List<double[]> d2)
        {
            List<double[]> d = new List<double[]>();

            for (int i = 0; i < d1.Count; i++)
            {
                double[] row = d1[i].Take(d1[i].Length - 1)
                    .Concat(d2[i].Take(d2[i].Length - 1))
                    .ToArray();

                double maxVal = Math.Max(d1[i].Last(), -d2[i].Last());
                row = row.Concat(new double[] { maxVal }).ToArray();

                d.Add(row);
            }

            return d;
        }
        public List<double[]> dUnion(List<double[]> d1, List<double[]> d2)
        {
            List<double[]> d = new List<double[]>();

            for (int i = 0; i < d1.Count; i++)
            {
                double[] row = d1[i].Take(d1[i].Length - 1)
                    .Concat(d2[i].Take(d2[i].Length - 1))
                    .ToArray();

                double minVal = Math.Min(d1[i].Last(), d2[i].Last());
                row = row.Concat(new double[] { minVal }).ToArray();

                d.Add(row);
            }

            return d;
        }
        public List<double[]> dLine(List<double[]> P, double x1, double y1, double x2, double y2)
        {
            double[] a = { x2 - x1, y2 - y1 };
            double normA = Math.Sqrt(a[0] * a[0] + a[1] * a[1]);
            a[0] /= normA;
            a[1] /= normA;

            List<double[]> d = new List<double[]>();

            foreach (var point in P)
            {
                double[] b = { point[0] - x1, point[1] - y1 };
                double distance = b[0] * a[1] - b[1] * a[0];
                d.Add(new double[] { distance, distance });
            }

            return d;
        }
        public List<double[]> dIntersect(List<double[]> d1, List<double[]> d2)
        {
            List<double[]> d = new List<double[]>();

            for (int i = 0; i < d1.Count; i++)
            {
                double[] row = d1[i].Take(d1[i].Length - 1)
                    .Concat(d2[i].Take(d2[i].Length - 1))
                    .ToArray();

                double maxVal = Math.Max(d1[i].Last(), d2[i].Last());
                row = row.Concat(new double[] { maxVal }).ToArray();

                d.Add(row);
            }

            return d;
        }
        public List<double[]> dEllipse(List<double[]> P, double xc, double yc, double a, double b)
        {
            List<double[]> d = new List<double[]>();

            foreach (var point in P)
            {
                double distance = (Math.Sqrt(Math.Pow((point[0] - xc) / a, 2) + Math.Pow((point[1] - yc) / b, 2)) - 1) / 2;
                d.Add(new double[] { distance, distance });
            }

            return d;
        }

    }
}
