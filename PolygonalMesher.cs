using DEMSoft.Drawing;
using DEMSoft.Drawing.Geometry;
using HullDelaunayVoronoi.Primitives;
using HullDelaunayVoronoi.Voronoi;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Reflection.Emit;
using System.Windows.Forms.VisualStyles;
using System.Xml.Linq;
using Color = System.Drawing.Color;
using System.Text.Json;
using System.IO;

namespace DEMSoft.PolygonalMesher
{
    public class PolygonalMesher
    {
        private int nElem;//number of seeds or elements
        private int maxIter;//max iterater to get a good mesh
        private Domain domain;
        private Domain domainReflex;
        private List<double[]> P;
        private List<double[]> PFixxed;
        private List<double[]> Pc;
        private List<double[]> nearBoundaryP;
        private List<double[]> reflexP;
        private double dOffset;
        private VoronoiMesh2 voronoi;
        private List<double[]> Node { get; set; }
        private List<int[]> Element { get; set; }
        private List<int[]> ElementCollapse;
        private List<double[]> NodeCollapse;
        private double c;
        private double Tol;
        private List<double[]> pointFix;

        public PolygonalMesher(Domain domain, int nElem, int MaxIter, List<double[]> PFixxed, double Tol, double c, List<double[]> pointFix = null, double dOffset = 0)
        {
            this.domain = domain;
            this.nElem = nElem;
            this.maxIter = MaxIter;
            this.PFixxed = PFixxed;
            this.dOffset = dOffset;
            this.c = c;
            this.Tol = Tol;
            this.pointFix = pointFix;

        }
        public List<double[]> DomainElement(string demand, List<double[]> Arg, List<double[]> Co)
        {
            double[] BdBox = { Co.Min(p => p[0]), Co.Max(p => p[0]), Co.Min(p => p[1]), Co.Max(p => p[1]) };

            switch (demand)
            {
                case "Dist":
                    return DistFnc(Arg, BdBox, Co);
                case "BdBox":
                    return new List<double[]>() { domain.GetBoundBox() };

                default:
                    throw new ArgumentException("Invalid demand.");
            }
        }
        public List<double[]> DistFnc(List<double[]> P, double[] BdBox, List<double[]> Co)
        {
            int nnel = Co.Count;
            int[,] TRI = new int[nnel, 2];
            for (int i = 0; i < nnel; i++)
            {
                TRI[i, 0] = i;
                TRI[i, 1] = (i + 1) % nnel;
            }
            Dist d = new Dist(BdBox);
            List<List<double[]>> Distance = new List<List<double[]>>();
            for (int i = 0; i < nnel; i++)
            {
                var ddd = d.dLine(P, Co[TRI[i, 0]][0], Co[TRI[i, 0]][1], Co[TRI[i, 1]][0], Co[TRI[i, 1]][1]);
                Distance.Add(ddd);
            }
            var dis = d.dIntersect(Distance[1], Distance[0]);
            for (int i = 2; i < nnel; i++)
            {
                dis = d.dIntersect(Distance[i], dis);
            }
            return dis;
        }

        public void GenerateMeshForRandomDomain(DenseMatrix Coordinate, out List<double[]> Node, out List<int[]> Element)
        {
            int numPointFixxed;
            if (this.PFixxed == null)
            {
                this.P = domain.GenerateSeedRandom(nElem, dOffset);
            }
            else
            {
                numPointFixxed = PFixxed.Count;
                int numPointNeedBeInserted = nElem - numPointFixxed;
                if (numPointNeedBeInserted >= 0)
                {
                    P = new List<double[]>();
                    P.AddRange(PFixxed);
                    P.AddRange(domain.GenerateSeedRandom(numPointNeedBeInserted, dOffset));
                }
            }

            this.nElem = P.Count;
            int It = 0;
            double Err = 1;
            double[] bdBox = domain.GetBoundBox();
            double area = (bdBox[1] - bdBox[0]) * (bdBox[3] - bdBox[2]);
            Pc = P.Select(s => new double[] { s[0], s[1] }).ToList();

            bool flag = true;
            Node = null;
            Element = null;

            while ((It <= maxIter) && (Err > 5e-4))
            {
                double alpha = c * Math.Sqrt(area / nElem);
                //if (PFixxed != null)
                //{
                //    var di = DomainElement("Dist", Pc, Coordinate.ToRowArrays().ToList());
                //    List<int> indInside = domain.FindIndexInside(di, alpha);///Use for class Dist
                //    List<double[]> PInside = new List<double[]>();//except fixxed points/// 
                //    for (int i = 0; i < indInside.Count; i++)
                //    {
                //        PInside.Add(Pc[indInside[i]]);
                //    }
                //    P = new List<double[]>();
                //    P.AddRange(PFixxed);
                //    P.AddRange(PInside);
                //}
                //else
                //{
                //    P = Pc.Select(s => new double[] { s[0], s[1] }).ToList();// Lloyd's update
                //}

                P = Pc.Select(s => new double[] { s[0], s[1] }).ToList();// Lloyd's update

                reflexP = GenerateReflexPoints(P, alpha, Coordinate);
                //reflexP = domain.GenerateSeedsReflex(P, alpha, nElem, out nearBoundaryP);   // Generate the reflections

                PolyMshr_FixedPoints(P, reflexP, pointFix);

                List<Vertex2> vertices = new List<Vertex2>();
                for (int i = 0; i < P.Count; i++)
                {
                    vertices.Add(new Vertex2((float)P[i][0], (float)P[i][1]));
                }
                for (int i = 0; i < reflexP.Count; i++)
                {
                    vertices.Add(new Vertex2((float)reflexP[i][0], (float)reflexP[i][1]));
                }

                ///////Generate Voronoi From Matlab function///////
                GenerateVoronoiFromMatlabFunction(vertices, out Node, out Element);
                //GenerateVoronoi(vertices, out Node, out Element);

                ComputeCndPoly(Element, Node, nElem, out Pc, out double[] A);

                area = A.Select(x => Math.Abs(x)).Sum();
                Err = CalculateErr(A, Pc, P, nElem, area);
                ///////////////////////
                //area = 0;
                //ComputeCentriodPolygon(Node, Element, out Pc, out double[] A);
                //List<int> Pin = new List<int>();
                //int countElement = Element.Count;
                //for (int iP = 0; iP < P.Count; iP++)
                //{

                //    for (int i = 0; i < countElement; i++)
                //    {
                //        int nv = Element[i].Length;
                //        double[] vx = new double[nv];
                //        double[] vy = new double[nv];

                //        for (int j = 0; j < nv; j++)
                //        {
                //            double[] vx_vy = Node[Element[i][j]];
                //            vx[j] = vx_vy[0];
                //            vy[j] = vx_vy[1];
                //        }
                //        if (IsPointInPolygon(P[iP][0], P[iP][1], vx, vy))
                //        {
                //            Pin.Add(i);

                //        }
                //    }
                //}

                //List<double[]> PcInside = new List<double[]>();
                //double[] AA = new double[Pin.Count];
                //for (int i = 0; i < Pin.Count; i++)
                //{
                //    AA[i] = A[Pin[i]];
                //    area += Math.Abs(A[Pin[i]]);
                //    PcInside.Add(Pc[Pin[i]]);
                //}
                //Err = CalculateErr(AA, PcInside, P, nElem, area);
                ////////////

                //for (int i = 0; i < Pc.Count; i++)
                //{
                //    if (domain.isInside(Pc[i], 0))
                //    {
                //        PcInside.Add(Pc[i]);
                //    }
                //}
                List<double[]> dist = new List<double[]>();
                dist = DomainElement("Dist", Pc, Coordinate.ToRowArrays().ToList());

                List<double[]> PcInside = new List<double[]>();

                List<int> indexPointInside = domain.FindIndexInside(dist, 0);///Use for class Dist
                for (int i = 0; i < indexPointInside.Count; i++)
                {
                    PcInside.Add(Pc[indexPointInside[i]]);
                }

                this.Pc = PcInside;
                this.Element = Element;
                this.Node = Node;

                It++;

            }
            List<double[]> pc = new List<double[]>();

            Node = this.Node;
            Element = this.Element;
            Element = SelectElementInside(Node, Element, out pc, out Node);

            this.Node = Node;
            this.Element = Element;

            NodeCollapse = this.Node;
            ElementCollapse = this.Element;

            PolyMshr_CllpsEdgs(Node, Element, Tol);

            Element = this.Element;
            Node = this.Node;

            ReOrderNodeAndElement();
            Node = this.Node;
            Element = this.Element;

            this.Node = Node;
            this.Element = Element;
        }

        private void ComputeCndPoly(List<int[]> element, List<double[]> node, int nElem, out List<double[]> Pc, out double[] A)
        {
            Pc = new List<double[]>();
            A = new double[nElem];
            for (int i = 0; i < nElem; i++)
            {
                int nv = element[i].Length;
                double[] vx = new double[nv];
                double[] vy = new double[nv];

                for (int j = 0; j < nv; j++)
                {
                    double[] vx_vy = node[element[i][j]];
                    vx[j] = vx_vy[0];
                    vy[j] = vx_vy[1];
                }
                double[] vxS = new double[nv];
                double[] vyS = new double[nv];

                for (int j = 0; j < nv; j++)
                {
                    vxS[j] = vx[(j + 1) % nv];
                    vyS[j] = vy[(j + 1) % nv];
                }

                double[] temp = Enumerable.Range(0, nv).Select(j => vx[j] * vyS[j] - vy[j] * vxS[j]).ToArray();

                A[i] = 0.5 * temp.Sum();

                double sum_vx_vxS_temp = 0.0;
                double sum_vy_vyS_temp = 0.0;

                for (int j = 0; j < nv; j++)
                {
                    sum_vx_vxS_temp += (vx[j] + vxS[j]) * temp[j];
                    sum_vy_vyS_temp += (vy[j] + vyS[j]) * temp[j];
                }
                double[] PcElement = new double[] { (1.0 / (6 * A[i])) * sum_vx_vxS_temp, (1.0 / (6 * A[i])) * sum_vy_vyS_temp };
                Pc.Add(PcElement);
            }

        }

        public void RunGenerateMesh(out List<double[]> Node, out List<int[]> Element, out List<int[]> ElementCollapse, out List<double[]> NodeCollapse)
        {
            int numPointFixxed;
            if (this.PFixxed == null)
            {
                this.P = domain.GenerateSeedRandom(nElem, dOffset);
            }
            else
            {
                numPointFixxed = PFixxed.Count;
                int numPointNeedBeInserted = nElem - numPointFixxed;
                if (numPointNeedBeInserted >= 0)
                {
                    P = new List<double[]>();
                    P.AddRange(PFixxed);
                    P.AddRange(domain.GenerateSeedRandom(numPointNeedBeInserted, dOffset));
                }
            }
            this.nElem = P.Count;
            int It = 0;
            double Err = 1;
            double[] bdBox = domain.GetBoundBox();
            double area = (bdBox[1] - bdBox[0]) * (bdBox[3] - bdBox[2]);

            Pc = P.Select(s => new double[] { s[0], s[1] }).ToList();
            bool flag = true;
            Node = null;
            Element = null;

            while ((It <= maxIter) && (Err > 5e-4))
            {
                double alpha = c * Math.Sqrt(area / nElem);
                if (PFixxed != null)
                {
                    List<double[]> PInside = domain.SelectPointsInside(Pc, alpha);//except fixxed points/// offset = alpha if use Pfix is contour
                    P = new List<double[]>();
                    P.AddRange(PFixxed);
                    P.AddRange(PInside);
                }
                else
                {
                    P = Pc.Select(s => new double[] { s[0], s[1] }).ToList();// Lloyd's update
                }

                reflexP = GenerateReflexPoints(P, alpha);
                //reflexP = domain.GenerateSeedsReflex(P, alpha, nElem, out nearBoundaryP);   // Generate the reflections

                PolyMshr_FixedPoints(P, reflexP, pointFix);

                List<Vertex2> vertices = new List<Vertex2>();
                for (int i = 0; i < P.Count; i++)
                {
                    vertices.Add(new Vertex2((float)P[i][0], (float)P[i][1]));
                }
                for (int i = 0; i < reflexP.Count; i++)
                {
                    vertices.Add(new Vertex2((float)reflexP[i][0], (float)reflexP[i][1]));
                }

                ///////Generate Voronoi From Matlab function///////
                //GenerateVoronoiFromMatlabFunction(vertices, out Node, out Element);

                GenerateVoronoi(vertices, out Node, out Element);

                ComputeCentriodPolygon(Node, Element, out Pc, out double[] A);
                List<double[]> PcInside = new List<double[]>();
                /////////////
                //for (int i = 0; i < Pc.Count; i++)
                //{
                //    if (domain.isInside(Pc[i], 0))
                //    {
                //        PcInside.Add(Pc[i]);
                //    }
                //}
                //////////////Dist///////////////
                List<double[]> dist = new List<double[]>();
                Dist d = new Dist(bdBox);
                dist = d.CantileverBeam(Pc);
                List<int> indexPointInside = domain.FindIndexInside(dist, dOffset);///Use for class Dist
                for (int i = 0; i < indexPointInside.Count; i++)
                {
                    PcInside.Add(Pc[indexPointInside[i]]);
                }
                ////////
                this.Pc = PcInside;
                this.Element = Element;
                this.Node = Node;
                ///////
                It++;

            }
            List<double[]> pc = new List<double[]>();


            Node = this.Node;
            Element = this.Element;
            Element = SelectElementInside(Node, Element, out pc, out Node);

            this.Node = Node;
            this.Element = Element;

            NodeCollapse = this.Node;
            ElementCollapse = this.Element;

            PolyMshr_CllpsEdgs(Node, Element, Tol);

            Element = this.Element;
            Node = this.Node;

            ReOrderNodeAndElement();
            Node = this.Node;
            Element = this.Element;

            this.Node = Node;
            this.Element = Element;
        }

        public List<double[]> GenerateReflexPoints(List<double[]> p, double alpha, DenseMatrix Coordinate = null)
        {
            if (Coordinate != null)
            {
                var dis = DomainElement("Dist", p, Coordinate.ToRowArrays().ToList());

                double eta = 0.9;
                double eps = 1e-8;

                double[][] diss = dis.Select(arr => arr.ToArray()).ToArray();


                var NBdrySegs = dis[0].Length - 1;
                List<double[]> reflexP = new List<double[]>();

                var sub_dis_n1 = DomainElement("Dist", p.Select(s => new double[] { s[0] + eps, s[1] }).ToList(), Coordinate.ToRowArrays().ToList());
                var sub_dis_n2 = DomainElement("Dist", p.Select(s => new double[] { s[0], s[1] + eps }).ToList(), Coordinate.ToRowArrays().ToList());
                double[][] result1 = sub_dis_n1.Select(arr => arr.ToArray()).ToArray();
                double[][] result2 = sub_dis_n2.Select(arr => arr.ToArray()).ToArray();
                var n1 = SubDomain.SubtractArrays(result1, diss).Select(row => row.Select(elem => elem / eps).ToArray()).ToArray();
                var n2 = SubDomain.SubtractArrays(result2, diss).Select(row => row.Select(elem => elem / eps).ToArray()).ToArray();


                double[,] I = new double[p.Count, NBdrySegs];
                List<double> disI = new List<double>();
                for (int i = 0; i < p.Count; i++)
                {
                    for (int j = 0; j < NBdrySegs; j++)
                    {
                        if (Math.Abs(dis[i][j]) < alpha)
                        {
                            double[] rp = new double[2];
                            I[i, j] = 1;
                            disI.Add(dis[i][j]);
                            rp[0] = p[i][0] - 2 * n1[i][j] * dis[i][j];
                            rp[1] = p[i][1] - 2 * n2[i][j] * dis[i][j];
                            reflexP.Add(rp);
                        }
                    }
                }
                var d_R_P = DomainElement("Dist", reflexP, Coordinate.ToRowArrays().ToList());


                List<double[]> R_P = new List<double[]>();

                for (int i = 0; i < reflexP.Count; i++)
                {
                    if (Math.Abs(d_R_P[i][NBdrySegs]) >= eta * Math.Abs(disI[i]) && d_R_P[i][NBdrySegs] > 0)
                    {
                        R_P.Add(reflexP[i]);
                    }
                }

                ////////////////
                R_P = SubDomain.GetUniqueRows(R_P).ToList();
                R_P.Sort((x, y) => x[0].CompareTo(y[0]));
                return R_P;

            }
            else
            {
                Dist dist = new Dist(domain.GetBoundBox());

                double eta = 0.9;
                double eps = 1e-8;
                var dis = dist.CantileverBeam(p);
                //var dis = dist.Cilyndrical(p);
                //var dis = dist.LShape(p);
                //var dis = dist.Photoelastic_Samples_1(p);
                //var dis = dist.Photoelastic_Samples_2(p);
                double[][] diss = dis.Select(arr => arr.ToArray()).ToArray();


                var NBdrySegs = dis[0].Length - 1;
                List<double[]> reflexP = new List<double[]>();

                var sub_dis_n1 = dist.CantileverBeam(p.Select(s => new double[] { s[0] + eps, s[1] }).ToList());
                var sub_dis_n2 = dist.CantileverBeam(p.Select(s => new double[] { s[0], s[1] + eps }).ToList());
                double[][] result1 = sub_dis_n1.Select(arr => arr.ToArray()).ToArray();
                double[][] result2 = sub_dis_n2.Select(arr => arr.ToArray()).ToArray();
                var n1 = SubDomain.SubtractArrays(result1, diss).Select(row => row.Select(elem => elem / eps).ToArray()).ToArray();
                var n2 = SubDomain.SubtractArrays(result2, diss).Select(row => row.Select(elem => elem / eps).ToArray()).ToArray();


                double[,] I = new double[p.Count, NBdrySegs];
                List<double> disI = new List<double>();
                for (int i = 0; i < p.Count; i++)
                {
                    for (int j = 0; j < NBdrySegs; j++)
                    {
                        if (Math.Abs(dis[i][j]) < alpha)
                        {
                            double[] rp = new double[2];
                            I[i, j] = 1;
                            disI.Add(dis[i][j]);
                            rp[0] = p[i][0] - 2 * n1[i][j] * dis[i][j];
                            rp[1] = p[i][1] - 2 * n2[i][j] * dis[i][j];
                            reflexP.Add(rp);
                        }
                    }
                }
                var d_R_P = dist.CantileverBeam(reflexP);

                List<double[]> R_P = new List<double[]>();

                for (int i = 0; i < reflexP.Count; i++)
                {
                    if (Math.Abs(d_R_P[i][NBdrySegs]) >= eta * Math.Abs(disI[i]) && d_R_P[i][NBdrySegs] > 0)
                    {
                        R_P.Add(reflexP[i]);
                    }
                }

                ////////////////
                R_P = SubDomain.GetUniqueRows(R_P).ToList();
                R_P.Sort((x, y) => x[0].CompareTo(y[0]));
                return R_P;
            }
        }

        private void ReOrderNodeAndElement()
        {
            List<int[]> ints = new List<int[]>();
            for (int i = 0; i < Node.Count; i++)
            {
                List<int> list = new List<int>();
                for (int j = 0; j < Node.Count && j != i; j++)
                {
                    var d = Math.Sqrt(Math.Pow(Node[i][0] - Node[j][0], 2) + Math.Pow(Node[i][1] - Node[j][1], 2));
                    if (d < 0.001)
                    {
                        list.Add(i);
                        list.Add(j);
                        list = list.Select(n => n).OrderBy(n => n).ToList();
                        ints.Add(list.ToArray());
                    }
                }
            }

            int[] cNode = Enumerable.Range(0, Node.Count).ToArray();
            foreach (var edge in ints)
            {
                cNode[edge[1]] = cNode[edge[0]];
            }
            RebuildList(Node, Element, cNode);
        }
        public static void SaveListToJson<T>(List<T> list, string filePath)
        {
            string json = JsonSerializer.Serialize(list);
            File.WriteAllText(filePath, json);
        }
        private void PolyMshr_FixedPoints(List<double[]> P, List<double[]> reflexP, List<double[]> pointFix)
        {
            if (pointFix != null)
            {
                List<double[]> PP = P.Concat(reflexP).ToList();
                for (int i = 0; i < pointFix.Count; i++)
                {
                    double[] B = new double[PP.Count];
                    for (int j = 0; j < PP.Count; j++)
                    {
                        B[j] = Math.Sqrt(Math.Pow(PP[j][0] - pointFix[i][0], 2) + Math.Pow(PP[j][1] - pointFix[i][1], 2));
                    }
                    int[] indices = Enumerable.Range(0, B.Length).ToArray();
                    Array.Sort(indices, (k, j) => B[k].CompareTo(B[j]));

                    for (int j = 1; j <= 3; j++)
                    {
                        int index = indices[j];
                        double[] n = new double[] { PP[index][0] - pointFix[i][0], PP[index][1] - pointFix[i][1] };
                        double normN = Math.Sqrt(n[0] * n[0] + n[1] * n[1]); // Norm
                        n[0] /= normN; // Normalize n
                        n[1] /= normN; // Normalize n

                        // Update PP
                        PP[index][0] -= n[0] * (B[indices[j]] - B[indices[0]]);
                        PP[index][1] -= n[1] * (B[indices[j]] - B[indices[0]]);
                    }
                }
                this.P = PP.GetRange(0, P.Count);
                this.reflexP = PP.GetRange(P.Count, PP.Count - P.Count);
            }
        }

        private void PolyMshr_ExtrNds(int nElem, List<double[]> node0, List<int[]> element0)
        {
            List<int> mapE = new List<int>();
            for (int i = 0; i < nElem; i++)
            {
                var ele = element0[i];
                mapE.AddRange(ele);
            }
            var map = mapE.Distinct().OrderBy(num => num).ToArray();
            int[] cNode = Enumerable.Range(0, node0.Count).Select(i => (int)i).ToArray();
            var diffElements = cNode.Except(map);

            // Check duplicate
            if (diffElements.Any())
            {              
                int maxMapValue = map.Max();
                foreach (var index in diffElements)
                {
                    cNode[index] = maxMapValue;
                }
            }
            RebuildList(node0, element0.GetRange(0, nElem), cNode);
        }
        private List<double[]> SelectInside(List<double[]> node1, List<int[]> element1, out List<double[]> finalNode, out List<int[]> element2)
        {
            List<double[]> node2 = new List<double[]>();
            element2 = new List<int[]>();
            List<double[]> distinctNodes = new List<double[]>();
            List<double[]> roundedNode = new List<double[]>();
            finalNode = new List<double[]>();
            List<int> index = new List<int>();
            foreach (var item in node1)
            {
                double roundedX = Math.Round(item[0], 2);
                double roundedY = Math.Round(item[1], 2);
                roundedNode.Add(new double[] { roundedX, roundedY });
            }
            for (int i = 0; i < roundedNode.Count; i++)
            {
                distinctNodes.Add(roundedNode[i]);
                index.Add(i);
            }
            node2 = distinctNodes.GroupBy(node => (node[0], node[1])).Select(group => group.First()).ToList();

            List<int> originalIndices = new List<int>();
            for (int j = 0; j < distinctNodes.Count; j++)
            {
                if (node2.Contains(distinctNodes[j]))
                {
                    originalIndices.Add(j);
                }
            }
            for (int i = 0; i < originalIndices.Count; i++)
            {
                finalNode.Add(node1[index[originalIndices[i]]]);
            }
            foreach (var ele in element1)
            {
                List<int> subElement = new List<int>();
                foreach (var edge in ele)
                {
                    for (int i = 0; i < node2.Count; i++)
                    {
                        if (roundedNode[edge][0] == node2[i][0] && roundedNode[edge][1] == node2[i][1])
                        {
                            if (!subElement.Contains(i))
                                subElement.Add(i);
                        }
                    }
                }
                element2.Add(subElement.ToArray());
            }
            return finalNode;
        }
        private void PolyMshr_CllpsEdgs(List<double[]> Node0, List<int[]> Element0, double tol)
        {
            while (true)
            {
                List<List<int>> cEdge = new List<List<int>>();

                for (int i = 0; i < Element0.Count; i++)
                {
                    int countNodeInElement = Element0[i].Length;
                    List<double> vx = new List<double>();
                    List<double> vy = new List<double>();
                    if (countNodeInElement < 4)
                        continue;

                    for (int j = 0; j < countNodeInElement; j++)
                    {
                        double[] vx_vy = Node0[Element0[i][j]];// Duyệt qua từng node của phần tử
                        vx.Add(vx_vy[0]);
                        vy.Add(vx_vy[1]);
                    }
                    int nv = vx.Count;

                    List<double> beta = new List<double>();
                    for (int k = 0; k < nv; k++)
                    {
                        var vySum = vy.Sum() / nv;
                        var vxSum = vx.Sum() / nv;
                        double betaValue = Math.Atan2(vy[k] - vySum, vx[k] - vxSum);
                        beta.Add(betaValue);
                    }

                    List<double> betaMod = new List<double>();
                    for (int ii = 0; ii < beta.Count; ii++)
                    {
                        int nextIndex = (ii + 1) % beta.Count; // Circle index
                        double diff = beta[nextIndex] - beta[ii];
                        double diffMod = (diff % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI);
                        betaMod.Add(diffMod);
                    }

                    double betaIdeal = 2 * Math.PI / countNodeInElement;

                    List<List<int>> edge = Element0[i].Zip(Element0[i].Skip(1).Concat(new[] { Element0[i][0] }),
                        (node1, node2) => new List<int> { node1, node2 }).ToList();

                    //Edge with small angle
                    for (int h = 0; h < nv; h++)
                    {
                        if (betaMod[h] < tol * betaIdeal)
                        {
                            cEdge.Add(edge[h]);
                        }
                    }
                }
                if (cEdge.Count == 0) { break; }

                Node = new List<double[]>();
                Element = new List<int[]>();
              
                cEdge = cEdge
                    .Select(edge => edge.OrderBy(nodeIndex => nodeIndex).ToList()) 
                    .OrderBy(edge => edge[0]) 
                    .Distinct()
                    .ToList();

                int[] cNode = Enumerable.Range(0, Node0.Count).ToArray();
                foreach (var edge in cEdge)
                {
                    cNode[edge[1]] = cNode[edge[0]];
                }
                RebuildList(Node0, Element0, cNode);

                break;
            }
        }
        private void RebuildList(List<double[]> Node0, List<int[]> Element0, int[] cNode)
        {
            Node = new List<double[]>();
            Element = new List<int[]>();

            ///////////////Build node//////////////////
            var foo = cNode.Distinct().ToArray(); // Unique value
            var ix = foo;
            var jx = cNode.Select(value => Array.IndexOf(foo, value)).ToArray();

            if (Node0.Count > ix.Length) { ix[ix.Length - 1] = cNode.Max(); }
            for (int i = 0; i < ix.Length; i++)
            {
                Node.Add(Node0[ix[i]]);
            }

            ///////////////Build element////////////////
            for (int j = 0; j < Element0.Count; j++)
            {
                var element = Element0[j];
                int[] newElement = new int[element.Length];

                foreach (var item2 in element.Select((ele, index) => new { ele, index }))
                {
                    var ele = item2.ele;
                    var index = item2.index;
                    for (int i = 0; i < jx.Length; i++)
                    {
                        if (ele == i)
                        {
                            newElement[index] = jx[i];
                            continue;
                        }
                    }
                }
                var uniqueElement = newElement.Distinct().OrderBy(x => x).ToArray();

                List<double> vx = new List<double>();
                List<double> vy = new List<double>();
                for (int i = 0; i < uniqueElement.Length; i++)
                {
                    double[] vx_vy = Node[uniqueElement[i]];
                    vx.Add(vx_vy[0]);
                    vy.Add(vx_vy[1]);
                }
                int nv = vx.Count;

                var iix = Enumerable.Range(0, vy.Count) 
                .OrderBy(index => Math.Atan2(vy[index] - vy.Average(), vx[index] - vx.Average()))
                .ToList();

                uniqueElement = iix.Select(index => uniqueElement[index]).ToArray(); 
                Element.Add(uniqueElement);
            }
        }
        private void GenerateVoronoiFromMatlabFunction(List<Vertex2> vertices, out List<double[]> nodes, out List<int[]> element)
        {
            nodes = new List<double[]>();
            element = new List<int[]>();
            double[,] vertice = new double[vertices.Count, 2];

            for (int i = 0; i < vertices.Count; i++)
            {
                vertice[i, 0] = vertices[i].X;
                vertice[i, 1] = vertices[i].Y;
            }
            string[,] options = { { "Qbb", "Qz" } };
            MLApp.MLApp matlab = new MLApp.MLApp();
            matlab.Execute(@"cd C:\Users\dinhd\Desktop\PolygonalMesher\voronoi");

            object result = null;
            matlab.Feval("voronoin", 2, out result, vertice, options);

            object[] res = result as object[];

            var a = res[0];
            var b = res[1];
            double[,] GetNode = a as double[,];
            for (int i = 0; i < GetNode.GetLength(0); i++)
            {
                nodes.Add(new double[] { GetNode[i, 0], GetNode[i, 1] });
            }
            object[,] ele = b as object[,];
            for (int i = 0; i < ele.Length; i++)
            {
                List<int> subE = new List<int>();
                var e = ele[i, 0];
                double[,] getE = e as double[,];
                List<int> resultList = getE.Cast<double>().Select(el => (int)(el - 1)).ToList();
                element.Add(resultList.ToArray());
            }
        }
        static double CalculateErr(double[] A, List<double[]> Pc, List<double[]> P, int NElem, double Area)
        {

            double[,] res = new double[NElem, 2];
            for (int i = 0; i < Pc.Count; i++)
            {
                for (int j = 0; j < Pc[i].Length; j++)
                {
                    res[i, j] = Math.Pow(Pc[i][j] - P[i][j], 2);
                }

            }
            double[] rowSums = Enumerable.Range(0, res.GetLength(0))
                                      .Select(i => Enumerable.Range(0, res.GetLength(1))
                                                             .Sum(j => res[i, j]))
                                      .ToArray();

            var squaredSum = A.Select(x => x * x).ToArray();
            double result = 0;
            for (int i = 0; i < A.Length; i++)
            {
                result += squaredSum[i] * rowSums[i];
            }
            result = Math.Sqrt(result) * NElem / Math.Pow(Area, 1.5);
            return result;
        }
        static bool IsPointInPolygon(double x, double y, double[] polyX, double[] polyY)
        {
            int n = polyX.Length;
            int crossings = 0;

            for (int i = 0, j = n - 1; i < n; j = i++)
            {
                if (((polyY[i] > y) != (polyY[j] > y)) &&
                    (x < (polyX[j] - polyX[i]) * (y - polyY[i]) / (polyY[j] - polyY[i]) + polyX[i]))
                {
                    crossings++;
                }
            }
            return crossings % 2 != 0; 
        }
        private List<double[]> FixListNodeAndElement()
        {
            List<double[]> newNode = new List<double[]>();
            List<int> indexInsertedNode = new List<int>();
            for (int i = 0; i < Element.Count; i++)
            {
                int[] arrayIndexNode = GetElement(i);
                int[] newArrayIndexNode = new int[arrayIndexNode.Length];
                for (int j = 0; j < arrayIndexNode.Length; j++)
                {
                    int indexOf = indexInsertedNode.IndexOf(arrayIndexNode[j]);
                    if (indexOf == -1)
                    {
                        indexInsertedNode.Add(arrayIndexNode[j]);
                        newNode.Add(Node[arrayIndexNode[j]]);
                        newArrayIndexNode[j] = newNode.Count - 1;
                    }
                    else
                    {
                        for (int k = 0; k < newNode.Count; k++)
                        {
                            double[] currentNode = Node[arrayIndexNode[j]];
                            double[] positionNewNode = newNode[k];
                            if ((Math.Abs(currentNode[0] - positionNewNode[0]) < 1e-6)
                              && (Math.Abs(currentNode[1] - positionNewNode[1]) < 1e-6))
                            {
                                newArrayIndexNode[j] = k;
                                break;
                            }
                        }
                    }
                }
                Element[i] = newArrayIndexNode;
            }
            return newNode;
        }

        public void RunGenerateHoneycombMesh(double size, double[] firstPoint, out List<double[]> Node, out List<int[]> Element)
        {
            this.P = domain.GenerateSeedHoneycomb(size, firstPoint, dOffset);//PolyMshr_RndPtSet();
            if (this.PFixxed != null)
            {
                P.AddRange(PFixxed);
            }
            //int numPointFixxed; 
            //if (this.PFixxed == null)
            //{
            //  this.P = domain.GenerateSeedHoneycomb(size, firstPoint);//PolyMshr_RndPtSet();
            //}
            //else
            //{
            //  numPointFixxed = PFixxed.Count;
            //  int numPointNeedBeInserted = nElem - numPointFixxed;
            //  if (numPointNeedBeInserted > 0)
            //  {
            //    P = new List<double[]>();
            //    P.AddRange(PFixxed);
            //    P.AddRange(domain.GenerateSeedRandom(numPointNeedBeInserted, dOffset));
            //  }
            //}
            this.nElem = P.Count;

            double Tol = 1e-15;
            int It = 0;
            double Err = 1;
            double c = 5.5;
            double[] bdBox = domain.GetBoundBox();
            double area = (bdBox[1] - bdBox[0]) * (bdBox[3] - bdBox[2]);
            double[] pointFix = new double[2] { 0.1, 0.1 };

            Pc = P.Select(s => new double[] { s[0], s[1] }).ToList();
            bool flag = true;
            Node = null;
            Element = null;

            while ((It <= maxIter) && (Err > Tol))
            {
                double alpha = c * Math.Sqrt(area / nElem);
                if (PFixxed != null)
                {
                    List<double[]> PInside = domain.SelectPointsInside(Pc, alpha);//except fixxed points
                    P = new List<double[]>();
                    P.AddRange(PFixxed);
                    P.AddRange(PInside);
                }
                else
                {
                    P = Pc.Select(s => new double[] { s[0], s[1] }).ToList();// Lloyd's update
                }

                reflexP = domain.GenerateSeedsReflex(P, alpha, nElem, out nearBoundaryP);   // Generate the reflections
                //PolyMshr_FixedPoints(P, reflexP, pointFix);

                List<Vertex2> vertices = new List<Vertex2>();
                for (int i = 0; i < P.Count; i++)
                {
                    vertices.Add(new Vertex2((float)P[i][0], (float)P[i][1]));
                }
                for (int i = 0; i < reflexP.Count; i++)
                {
                    vertices.Add(new Vertex2((float)reflexP[i][0], (float)reflexP[i][1]));
                }

                GenerateVoronoi(vertices, out Node, out Element);
                ComputeCentriodPolygon(Node, Element, out Pc, out double[] A);

                Element = SelectElementInside(Node, Element, out Pc, out Node);
                //RebuildList(Node, Element, out Node, out Element);
                this.Node = Node;
                this.Element = Element;

                //double dPPc = 0;
                //for (int i = 0; i < P.Count; i++)
                //{
                //  dPPc += Math.Pow(Pc[i][0] - P[i][0], 2) + Math.Pow(Pc[i][1] - P[i][1], 2);
                //}
                //for (int i = 0; i < reflexP.Count; i++)
                //{
                //  dPPc += Math.Pow(Pc[i + P.Count - 1][0] - reflexP[i][0], 2) + Math.Pow(Pc[i + P.Count - 1][1] - reflexP[i][1], 2);
                //}
                //Err = Math.Sqrt(dPPc);

                //double[] sumPc = new double[A.Length];
                //for (int i = 0; i < A.Length; i++)
                //{
                //  sumPc[i] += Math.Pow(Pc[i][0] - P[i][0], 2) + Math.Pow(Pc[i][1] - P[i][1], 2);
                //}
                //area = 0;
                //double sumErr = 0;
                //for (int i = 0; i < A.Length; i++)
                //{
                //  area += A[i];
                //  sumErr += (A[i] * A[i]) * sumPc[i];
                //}
                //Err = Math.Sqrt(sumErr) * nElem / Math.Pow(area, 1.5);
                It++;
                //break;
            }
            Node = this.Node;
            Element = this.Element;
            PolyMshr_CllpsEdgs(Node, Element, Tol);
            Node = this.Node;
            Element = this.Element;


            SelectInside(Node, Element, out Node, out Element);
            this.Node = Node;
            this.Element = Element;


        }
        private List<int[]> SelectElementInside(List<double[]> node0, List<int[]> Element0, out List<double[]> PcInside, out List<double[]> nodeInside)//2_11_2023
        {
            List<int[]> selectedElementInsideDomain = new List<int[]>();
            List<int[]> elementInside = new List<int[]>();
            List<int[]> selectedElementOutsideDomain = new List<int[]>();
            PcInside = new List<double[]>();
            List<double[]> PcOutside = new List<double[]>();
            nodeInside = new List<double[]>();
            List<double[]> node2 = new List<double[]>();
            List<int[]> element1 = new List<int[]>();
            List<int> indexInside = new List<int>();
            List<int> indexOutside = new List<int>();
            ComputeCentriodPolygon(node0, Element0, out List<double[]> Pc, out double[] A);
            for (int i = 0; i < Pc.Count; i++)
            {
                if (domain.isInside(Pc[i], 0))
                {
                    selectedElementInsideDomain.Add(Element0[i]);
                    PcInside.Add(Pc[i]);
                    indexInside.Add(i);
                }
                else
                {
                    PcOutside.Add(Pc[i]);
                    selectedElementOutsideDomain.Add(Element0[i]);
                    indexOutside.Add(i);
                }
            }
            node2.AddRange(node0);
            List<int> nodeOut = new List<int>();
            foreach (var ele in selectedElementOutsideDomain)
            {
                var numEdge = ele.Length;
                for (int i = 0; i < numEdge; i++)
                {
                    nodeOut.Add(ele[i]);
                }
            }
            nodeOut.Distinct();
            List<int> nodeIn = new List<int>();
            foreach (var ele in selectedElementInsideDomain)
            {
                var numEdge = ele.Length;
                for (int i = 0; i < numEdge; i++)
                {
                    nodeIn.Add(ele[i]);
                }
            }
            nodeIn.Distinct();
            List<int> uniqueIndicesInNodeOut = nodeOut.Except(nodeIn).ToList();
            nodeInside = node2.Where((item, index) => !uniqueIndicesInNodeOut.Contains(index)).ToList();
            foreach (var ele in selectedElementInsideDomain)
            {
                List<int> subElement = new List<int>();
                foreach (var edge in ele)
                {
                    for (int i = 0; i < nodeInside.Count; i++)
                    {
                        if (node2[edge][0] == nodeInside[i][0] && node2[edge][1] == nodeInside[i][1])
                        {
                            subElement.Add(i);
                        }
                    }
                }
                elementInside.Add(subElement.ToArray());
            }
            return elementInside;
        }

        private void GenerateVoronoi(List<Vertex2> vertices, out List<double[]> node, out List<int[]> Element)
        {
            voronoi = new VoronoiMesh2();
            voronoi.Generate(vertices);



            node = new List<double[]>();
            Element = new List<int[]>();
            for (int i = 0; i < voronoi.Cells.Count; i++)
            {
                node.Add(new double[] { voronoi.Cells[i].CircumCenter.X, voronoi.Cells[i].CircumCenter.Y });
            }
            for (int i = 0; i < voronoi.Regions.Count; i++)
            {
                if (voronoi.Regions[i].Edges.Count < 3)
                    continue;

                List<double[]> vc = GetNodeOfElement(i);
                int nv = vc.Count;
                double vxc = 0;
                double vyc = 0;
                for (int j = 0; j < nv; j++)
                {
                    vxc += vc[j][0];
                    vyc += vc[j][1];
                }
                double[] PcElement = new double[] { vxc / nv, vyc / nv };

                double[] anglePoint = new double[nv];
                for (int j = 0; j < nv; j++)
                {
                    anglePoint[j] = Math.Atan2(vc[j][1] - PcElement[1], vc[j][0] - PcElement[0]);
                }
                double[] sortedAnglePoint = (double[])anglePoint.Clone();
                Array.Sort(sortedAnglePoint);
                int[] indexSort = new int[nv];
                List<int> arrayNode = new List<int>();
                List<double> listAngleAdded = new List<double>();
                for (int j = 0; j < nv; j++)
                {
                    indexSort[j] = Array.IndexOf(anglePoint, sortedAnglePoint[j]);
                    int indexSortj = Array.IndexOf(anglePoint, sortedAnglePoint[j]);
                    bool isFlag = true;
                    for (int k = 0; k < listAngleAdded.Count; k++)
                    {
                        if (Math.Abs((listAngleAdded[k] - sortedAnglePoint[j]) / listAngleAdded[k]) < 1e-4)
                        {
                            isFlag = false;
                            break;
                        }
                    }
                    if (isFlag)
                    {
                        arrayNode.Add(voronoi.Regions[i].Cells[indexSortj].CircumCenter.Id);
                        listAngleAdded.Add(sortedAnglePoint[j]);
                    }
                }
                if (arrayNode.Count < 3)
                    continue;
                //arrayNode.Add(voronoi.Regions[i].Cells[indexSort[0]].CircumCenter.Id);
                //for (int j = 0; j < nv; j++)
                //{
                //  arrayNode.Add(voronoi.Regions[i].Cells[indexSort[j]].CircumCenter.Id);
                //}

                //List<int[]> allPairEdge = new List<int[]>();
                //for (int j = 0; j < voronoi.Regions[i].Edges.Count; j++)
                //{
                //  allPairEdge.Add(new int[] { voronoi.Regions[i].Edges[j].From.CircumCenter.Id, voronoi.Regions[i].Edges[j].To.CircumCenter.Id });
                //}
                //arrayNode.Add(allPairEdge[0][0]);
                //arrayNode.Add(allPairEdge[0][1]);
                //allPairEdge.RemoveAt(0);
                //while (allPairEdge.Count > 0)
                //{
                //  bool flag = false;
                //  bool isLast = true;
                //  int j = 0;
                //  for (j = 0; j < allPairEdge.Count; j++)
                //  {
                //    if (arrayNode.Contains(allPairEdge[j][0]) && arrayNode.Contains(allPairEdge[j][1]))
                //    {
                //      allPairEdge.RemoveAt(j);
                //      break;
                //    }
                //    if (allPairEdge[j].Contains(arrayNode[arrayNode.Count - 1]))
                //    {
                //      isLast = true;
                //      flag = true;
                //      break;
                //    }
                //    if (allPairEdge[j].Contains(arrayNode[0]))
                //    {
                //      isLast = false;
                //      flag = true;
                //      break;
                //    }
                //  }
                //  if (flag)
                //  {
                //    if (isLast)
                //    {
                //      if (allPairEdge[j][0] == arrayNode[arrayNode.Count - 1])
                //        arrayNode.Add(allPairEdge[j][1]);
                //      else
                //        arrayNode.Add(allPairEdge[j][0]);
                //    }
                //    else
                //    {
                //      if (allPairEdge[j][0] == arrayNode[0])
                //        arrayNode.Insert(0, allPairEdge[j][1]);
                //      else
                //        arrayNode.Insert(0, allPairEdge[j][0]);
                //    }
                //    allPairEdge.RemoveAt(j);
                //  }
                //}
                Element.Add(arrayNode.ToArray());
            }
        }
        public List<double[]> GetSeeds()
        { return P; }
        public List<double[]> GetCentriodPolygon()
        { return Pc; }
        public List<double[]> GetSeedsNearBoundary()
        { return nearBoundaryP; }
        public List<double[]> GetSeedsReflex()
        { return reflexP; }
        public VoronoiMesh2 GetVoronoi()
        { return voronoi; }

        #region DrawTestPoly
        //public void Draw(ViewerForm viewer)
        //{
        //    double offset = 0.1;
        //    //foreach (VoronoiRegion<Vertex2> region in voronoi.Regions)
        //    //{
        //    //  bool draw = false;
        //    //  foreach (DelaunayCell<Vertex2> cell in region.Cells)
        //    //  {
        //    //    //if (!domain.isInside(new double[] { cell.CircumCenter.X, cell.CircumCenter.Y }, -0.001))
        //    //    if (domain.isInsideBoundBox(new double[] { cell.CircumCenter.X, cell.CircumCenter.Y }, offset))
        //    //    {
        //    //      draw = true;
        //    //      break;
        //    //    }
        //    //  }
        //    //  if (!draw) continue;

        //    //  //for (int i = 0; i < region.Edges.Count; i++)
        //    //  //{
        //    //  //  double x1 = region.Edges[i].From.CircumCenter.X;
        //    //  //  double y1 = region.Edges[i].From.CircumCenter.Y;
        //    //  //  double x2 = region.Edges[i].To.CircumCenter.X;
        //    //  //  double y2 = region.Edges[i].To.CircumCenter.Y;

        //    //  //  Line l = new Line(x1, y1, 0, x2, y2, 0);
        //    //  //  viewer.AddObject3D(l);
        //    //  //}
        //    //  for (int i = 0; i < region.Cells.Count; i++)
        //    //  {
        //    //    double x1 = 0;
        //    //    double y1 = 0;
        //    //    double x2 = 0;
        //    //    double y2 = 0;
        //    //    if (i == region.Cells.Count - 1)
        //    //    {
        //    //      x1 = region.Cells[i].CircumCenter.X;
        //    //      y1 = region.Cells[i].CircumCenter.Y;
        //    //      x2 = region.Cells[0].CircumCenter.X;
        //    //      y2 = region.Cells[0].CircumCenter.Y;
        //    //    }
        //    //    else
        //    //    {
        //    //      x1 = region.Cells[i].CircumCenter.X;
        //    //      y1 = region.Cells[i].CircumCenter.Y;
        //    //      x2 = region.Cells[i + 1].CircumCenter.X;
        //    //      y2 = region.Cells[i + 1].CircumCenter.Y;
        //    //    }
        //    //    Line l = new Line(x1, y1, 0, x2, y2, 0);
        //    //    viewer.AddObject3D(l);
        //    //  }

        //    //  List<double> xx = new List<double>();
        //    //  List<double> yy = new List<double>();
        //    //  for (int i = 0; i < voronoi.Cells.Count; i++)
        //    //  {
        //    //    double x = voronoi.Cells[i].CircumCenter.X;
        //    //    double y = voronoi.Cells[i].CircumCenter.Y;
        //    //    xx.Add(x);
        //    //    yy.Add(y);
        //    //    //if (domain.isInside(new double[] { x, y }, -0.001))
        //    //    //if (domain.isInsideBoundBox(new double[] { x, y }, offset))
        //    //    //{
        //    //    //  xx.Add(x);
        //    //    //  yy.Add(y);
        //    //    //}
        //    //  }
        //    //  PointSet ps = new PointSet(xx.ToArray(), yy.ToArray(), new double[xx.Count]);
        //    //  ps.SetPointSize(5);
        //    //  viewer.AddObject3D(ps);
        //    //}

        //    foreach (int[] region in Element)
        //    {
        //        bool draw = false;
        //        foreach (int nodeIndex in region)
        //        {
        //            //if (!domain.isInside(new double[] { cell.CircumCenter.X, cell.CircumCenter.Y }, -0.001))
        //            if (domain.isInsideBoundBox(new double[] { Node[nodeIndex][0], Node[nodeIndex][1] }, offset))
        //            {
        //                draw = true;
        //                break;
        //            }
        //        }
        //        if (!draw) continue;

        //        for (int i = 0; i < region.Length; i++)
        //        {
        //            double x1 = Node[region[i]][0];
        //            double y1 = Node[region[i]][1];
        //            double x2, y2;
        //            if (i == region.Length - 1)
        //            {
        //                x2 = Node[region[0]][0];
        //                y2 = Node[region[0]][1];
        //            }
        //            else
        //            {
        //                x2 = Node[region[i + 1]][0];
        //                y2 = Node[region[i + 1]][1];
        //            }
        //            Line l = new Line(x1, y1, 0, x2, y2, 0);
        //            viewer.AddObject3D(l);
        //        }
        //    }
        //    List<double> xx = new List<double>();
        //    List<double> yy = new List<double>();
        //    for (int i = 0; i < Node.Count; i++)
        //    {
        //        double x = voronoi.Cells[i].CircumCenter.X;
        //        double y = voronoi.Cells[i].CircumCenter.Y;
        //        xx.Add(x);
        //        yy.Add(y);
        //    }
        //    PointSet ps = new PointSet(xx.ToArray(), yy.ToArray(), new double[xx.Count]);
        //    ps.SetPointSize(5);
        //    viewer.AddObject3D(ps);
        //}
        #endregion

        public void Draw(ViewerForm viewer, List<double[]> Node, List<int[]> Element) //Check collapse
        {
            int c = 0;
            foreach (int[] region in Element)
            {
                for (int i = 0; i < region.Length; i++)
                {
                    double x1 = Node[region[i]][0];
                    double y1 = Node[region[i]][1];
                    double x2, y2;
                    if (i == region.Length - 1)
                    {
                        x2 = Node[region[0]][0];
                        y2 = Node[region[0]][1];
                    }
                    else
                    {
                        x2 = Node[region[i + 1]][0];
                        y2 = Node[region[i + 1]][1];
                    }
                    Line l = new Line(x1, y1, 0, x2, y2, 0);
                    l.ColorObject = Color.Blue;
                    l.SetWidth(2);
                    viewer.AddObject3D(l);
                }
                c++;
                //if (c == 2)
                //  break;
            }
            List<double> xx = new List<double>();
            List<double> yy = new List<double>();
            for (int i = 0; i < Node.Count; i++)
            {
                xx.Add(Node[i][0]);
                yy.Add(Node[i][1]);
            }
            PointSet ps = new PointSet(xx.ToArray(), yy.ToArray(), new double[xx.Count]);
            ps.SetPointSize(5);
            ps.ColorObject = Color.Black;
            viewer.AddObject3D(ps);
        }
        public void Draw(ViewerForm viewer)
        {
            int c = 0;
            foreach (int[] region in Element)
            {
                for (int i = 0; i < region.Length; i++)
                {
                    double x1 = Node[region[i]][0];
                    double y1 = Node[region[i]][1];
                    double x2, y2;
                    if (i == region.Length - 1)
                    {
                        x2 = Node[region[0]][0];
                        y2 = Node[region[0]][1];
                    }
                    else
                    {
                        x2 = Node[region[i + 1]][0];
                        y2 = Node[region[i + 1]][1];
                    }
                    Line l = new Line(x1, y1, 0, x2, y2, 0);
                    l.ColorObject = Color.Black;
                    l.SetWidth(2);
                    viewer.AddObject3D(l);
                }
                c++;
                //if (c == 2)
                //  break;
            }
            List<double> xx = new List<double>();
            List<double> yy = new List<double>();
            for (int i = 0; i < Node.Count; i++)
            {
                xx.Add(Node[i][0]);
                yy.Add(Node[i][1]);
            }
            PointSet ps = new PointSet(xx.ToArray(), yy.ToArray(), new double[xx.Count]);
            ps.SetPointSize(1);
            ps.ColorObject = Color.Black;

            viewer.AddObject3D(ps);
        }
        private void ComputeCentriodPolygon(List<double[]> Node, List<int[]> Element, out List<double[]> Pc, out double[] A)
        {

            Pc = new List<double[]>();
            int countElement = Element.Count;
            A = new double[countElement];
            for (int i = 0; i < countElement; i++)
            {
                int nv = Element[i].Length;
                double[] vx = new double[nv];
                double[] vy = new double[nv];

                for (int j = 0; j < nv; j++)
                {
                    double[] vx_vy = Node[Element[i][j]];
                    vx[j] = vx_vy[0];
                    vy[j] = vx_vy[1];
                }
                double[] vxS = new double[nv];
                double[] vyS = new double[nv];

                for (int j = 0; j < nv; j++)
                {
                    vxS[j] = vx[(j + 1) % nv];
                    vyS[j] = vy[(j + 1) % nv];
                }

                double[] temp = Enumerable.Range(0, nv).Select(j => vx[j] * vyS[j] - vy[j] * vxS[j]).ToArray();

                A[i] = 0.5 * temp.Sum();

                double sum_vx_vxS_temp = 0.0;
                double sum_vy_vyS_temp = 0.0;

                for (int j = 0; j < nv; j++)
                {
                    sum_vx_vxS_temp += (vx[j] + vxS[j]) * temp[j];
                    sum_vy_vyS_temp += (vy[j] + vyS[j]) * temp[j];
                }
                double[] PcElement = new double[] { (1.0 / (6 * A[i])) * sum_vx_vxS_temp, (1.0 / (6 * A[i])) * sum_vy_vyS_temp };
                Pc.Add(PcElement);
            }
        }
        private void ComputeCentriodPolygonRand(List<double[]> Node, List<int[]> Element, out List<double[]> Pc, out double[] A)
        {

            Pc = new List<double[]>();
            int countElement = Element.Count;
            A = new double[countElement];
            for (int i = 0; i < countElement; i++)
            {
                int nv = Element[i].Length;
                double[] vx = new double[nv];
                double[] vy = new double[nv];

                for (int j = 0; j < nv; j++)
                {
                    double[] vx_vy = Node[Element[i][j]];
                    vx[j] = vx_vy[0];
                    vy[j] = vx_vy[1];
                }
                double[] vxS = new double[nv];
                double[] vyS = new double[nv];

                for (int j = 0; j < nv; j++)
                {
                    vxS[j] = vx[(j + 1) % nv];
                    vyS[j] = vy[(j + 1) % nv];
                }

                double[] temp = Enumerable.Range(0, nv).Select(j => vx[j] * vyS[j] - vy[j] * vxS[j]).ToArray();

                A[i] = 0.5 * temp.Sum();

                double sum_vx_vxS_temp = 0.0;
                double sum_vy_vyS_temp = 0.0;

                for (int j = 0; j < nv; j++)
                {
                    sum_vx_vxS_temp += (vx[j] + vxS[j]) * temp[j];
                    sum_vy_vyS_temp += (vy[j] + vyS[j]) * temp[j];
                }
                double[] PcElement = new double[] { (1.0 / (6 * A[i])) * sum_vx_vxS_temp, (1.0 / (6 * A[i])) * sum_vy_vyS_temp };
                Pc.Add(PcElement);
            }
        }
        public double[] ComputeCentriodElement(int idElement)
        {
            double vxc = 0;
            double vyc = 0;
            int nv = Element[idElement].Length;
            for (int j = 0; j < nv; j++)
            {
                int indexNode = Element[idElement][j];
                vxc += Node[indexNode][0];
                vyc += Node[indexNode][1];
            }
            return new double[] { vxc / nv, vyc / nv };
        }
        private List<double[]> GetNodeOfElement(int idElement, List<double[]> Node, List<int[]> Element)
        {
            List<double[]> vc = new List<double[]>();
            int[] listNode = Element[idElement];
            for (int i = 0; i < listNode.Length; i++)
            {
                vc.Add(new double[] { Node[listNode[i]][0], Node[listNode[i]][1] });
            }
            return vc;
        }
        private List<double[]> GetNodeOfElement(int idElement)
        {
            List<double[]> vc = new List<double[]>();
            VoronoiRegion<Vertex2> region = voronoi.Regions[idElement];
            for (int i = 0; i < region.Cells.Count; i++)
            {
                double x = region.Cells[i].CircumCenter.X;
                double y = region.Cells[i].CircumCenter.Y;
                vc.Add(new double[] { x, y });
            }
            return vc;
        }
        private int CountRegionPolygon()
        { return voronoi.Regions.Count; }
        public int[] GetElement(int index)
        {
            return Element[index];
        }
        public int CountElement()
        {
            return Element.Count;
        }
        public int CountNode() { return Node.Count; }
        public double[] GetNode(int index)
        {
            return Node[index];
        }
        public void UpdateListNode(List<double[]> list)
        {
            Node = list;
        }
        public void UpdateListElement(List<int[]> list)
        {
            Element = list;
        }
        public List<int[]> GetListElement() { return Element; }
        public List<double[]> GetListNode() { return Node; }


    }
}
