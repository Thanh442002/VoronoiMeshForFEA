using DEMSoft.Drawing;
using DEMSoft.NURBS;
using DEMSoft.PolygonalMesher;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;


namespace DemoMesh
{
    internal static class Program
    {
        static void Main() 
        {
            //If use another problem, go to class dist in library PolygonalMesher to change the distance function//
            //Replace dis function before by your function to copmpute the distance //
            /////////////////Geomertry BEAM///////////////////
            List<Abstract1DParametricGeometry> curves = new List<Abstract1DParametricGeometry>();
            curves.Add(GeometryCreator.CreateStraightNURBSCurve(0, -0.5, 0, 0, 0.5, 0));
            curves.Add(GeometryCreator.CreateStraightNURBSCurve(0, 0.5, 0, 3, 0.5, 0));
            curves.Add(GeometryCreator.CreateStraightNURBSCurve(3, 0.5, 0, 3, -0.5, 0));
            curves.Add(GeometryCreator.CreateStraightNURBSCurve(3, -0.5, 0, 0, -0.5, 0));
            NURBSCurve mergedCurve = (NURBSCurve)Abstract1DParametricGeometry.MergeGeometry(curves);
            List<Abstract1DParametricGeometry> curveOutside = new List<Abstract1DParametricGeometry>();
            curveOutside.Add(mergedCurve);
            SubDomain subDomain = new SubDomain(curves, curves);


            ////////////Draw Geomertry////////////
            ViewerForm viewGeo = new ViewerForm(true);
            List<double[]> pointContuor = new List<double[]>();

            foreach (var cur in curves)
            {
                cur.resolution = 30;
                cur.isDrawControlNet = false;
                cur.isDrawControlPoint = false;
                cur.isDrawKnot = false;
                cur.colorCurve = Color.Yellow;
                cur.Draw(viewGeo);
            }

            viewGeo.UpdateCamera();
            viewGeo.Run();

            //////////////Domain and SubDomain for meshing///////////////
            List<SubDomain> listSub = new List<SubDomain>();
            listSub.Add(subDomain);

            Domain domain = new Domain(listSub);
            var bdbox = domain.GetBoundBox();                    
            //////////////Meshing////////////////
            double tol = 0.2;
            double c = 2.5;
            List<double[]> pointFix = new List<double[]>();
            pointFix.Add(new double[] { 3, 0 });////beam
            PolygonalMesher mesh = new PolygonalMesher(domain, 200, 50, null, tol, c, pointFix, 0);
            mesh.RunGenerateMesh(out List<double[]> Node, out List<int[]> Element, out List<int[]> ElementCollapse, out List<double[]> NodeCollapse);

            Console.WriteLine(Node.Count);
            Console.WriteLine(Element.Count);
            /////////////Plot Mesh///////////////
            ViewerForm viewMesh = new ViewerForm(true);
            mesh.Draw(viewMesh);
            viewMesh.UpdateCamera();
            viewMesh.Run();


        }
    }
}
