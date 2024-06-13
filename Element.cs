using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DEMSoft.PolygonalMesher
{
  internal class Element
  {
    public List<Node> Nodes { get; set; }
    public int CountNode()
    { return Nodes.Count; }

    public double[] Pc { get; }
    public void ComputeCentriod()
    {
      int nv = CountNode();
      double vxc = 0;
      double vyc = 0;
      for (int i = 0; i < nv; i++)
      {
        vxc += vc[i][0];
        vyc += vc[i][1];
      }
      double[] PcElement = new double[] { vxc / nv, vyc / nv };
    }
  }
}
