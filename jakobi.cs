using math=System.Math;
//using matrix=dmatrix;

public class jacobi{

  public static int eigen(matrix A, vector e)
  {
    return eigen(A,e,null);
  }
  public static int eigen(matrix A, vector e, matrix V)
  {
    int n=A.rows;
    bool V_is_given = (V != null);
    for(int i=0;i<n;i++)
    {
      e[i]=A[i,i];
      if(V_is_given)
      {
        V[i,i]=1;
        for(int j=i+1;j<n;j++) V[i,j]=V[j,i]=0;
      }
    }
    bool changed; int sweeps=0;
    do {
      sweeps++; changed=false;
      for(int q,p=0;p<n;p++)for(q=p+1;q<n;q++)
        {
          double app = e[p], aqq = e[q], apq = A[p,q];
          double phi=0.5*math.Atan2(2*apq,aqq-app);
          double s=math.Sin(phi), c=math.Cos(phi);
          /*					
double c,s,cot2f,tanf;
if(apq==0){s=0;c=1;}
else{
  cot2f=(aqq-app)/2/apq;
  if(math.Abs(cot2f)>1) tanf=cot2f*(math.Sqrt(1+1/cot2f/cot2f)-1);
  else tanf=math.Sign(cot2f)*(math.Sqrt(cot2f*cot2f+1)-math.Abs(cot2f));
  c=1/math.Sqrt(1+tanf*tanf); s=c*tanf;
}
*/
          double app1=c*c*app-2*s*c*apq+s*s*aqq;
          double aqq1=s*s*app+2*s*c*apq+c*c*aqq;
          if(app1 != app || aqq1 != aqq)
          {
            changed=true;
            e[p]=app1;
            e[q]=aqq1;
            A[p,q]=0.0;
            for(int i=0;i<p;i++)
            {
              double aip=A[i,p], aiq=A[i,q];
              A[i,p]=c*aip-s*aiq;
              A[i,q]=c*aiq+s*aip;
            }
            for(int i=p+1;i<q;i++)
            {
              double api=A[p,i], aiq=A[i,q];
              A[p,i]=c*api-s*aiq;
              A[i,q]=c*aiq+s*api;
            }
            for(int i=q+1;i<n;i++)
            {
              double api=A[p,i], aqi=A[q,i];
              A[p,i]=c*api-s*aqi;
              A[q,i]=c*aqi+s*api;
            }
            if(V_is_given)for(int i=0;i<n;i++)
              {
                double vip=V[i,p], viq=V[i,q];
                V[i,p]=c*vip-s*viq;
                V[i,q]=c*viq+s*vip;
              }
          }
        }
    }while(changed);
    return sweeps;
  }
}
