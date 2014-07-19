using System;

public class matrix{

  private int nrows, ncols;
  private double[] data;

  public matrix(int _nrows, int _ncols){
    nrows=_nrows; ncols=_ncols;
    data = new double[nrows*ncols];
  }

  public matrix(string s){
    string[] rows = s.Split(';');
    nrows = rows.Length;
    ncols = rows[0].Split(' ',',').Length;
    data = new double[nrows*ncols];
    for(int i=0;i<nrows;i++){
      string[] ws = rows[i].Split(' ',',');
      if(ws.Length!=ncols)
        throw new System.ArgumentException("matrix: bad init-string");
      for(int j=0; j<ncols; j++){
        this[i,j]=System.Convert.ToDouble(ws[j]);
      }
    }
  }

  public double this[int r,int c]{
    get{return data[r+c*nrows];}
    set{data[r+c*nrows]=value;}
  }

  public int rows{get{return nrows;}set{}}
  public int cols{get{return ncols;}set{}}
  public int size1{get{return nrows;}set{}}
  public int size2{get{return ncols;}set{}}

  public vector this[string c,int i]{
    get{vector v = new vector();
      v.size=nrows;
      v.first=i*nrows;
      v.data=data;
      return v;}
    set{for(int k=0;k<nrows;k++)this[k,i]=value[k];}
  }

  public vector this[int i]{
    get{
      vector v = new vector();
      v.size=nrows;
      v.first=i*nrows;
      v.data=data;
      return v;
    }
    set{for(int k=0;k<nrows;k++)this[k,i]=value[k];}
  }

  public vector col(int i){
//	get{
    vector v = new vector();
    v.size=nrows; v.first=i*nrows; v.data=data; return v;
//	}
//	set{for(int k=0;k<nrows;k++)data[k+i*nrows]=value[k];}
  }

  public static matrix operator+ (matrix a, matrix b){
    matrix c = new matrix(a.rows,a.cols);
    for(int i=0;i<a.rows;i++)
      for(int j=0;j<a.cols;j++)
        c[i,j]=a[i,j]+b[i,j];
    return c;
  }

  public static matrix operator- (matrix a, matrix b){
    System.Diagnostics.Debug.Assert(a.size1==b.size1 && a.size2==b.size2);
    matrix c = new matrix(a.size1,a.size2);
    for(int i=0;i<a.size1;i++)
      for(int j=0;j<a.size2;j++)
        c[i,j]=a[i,j]-b[i,j];
    return c;
  }

  public static matrix operator* (matrix a, matrix b){
    int n=a.rows, m=b.cols;
    if(a.cols != b.rows) System.Console.Error.Write("matrix mismatch\n");
    var c = new matrix(n,m);
    for (int k=0;k<a.cols;k++) for (int j=0;j<m;j++){
        double tmp=b.data[k+j*b.rows];
//		for (int i=0;i<n;i++){
        for (int jn=j*n, kn=k*n; jn<j*n+n;){
//       c.data[i+j*n]+=a.data[i+k*n]*tmp;
          c.data[jn++]+=a.data[kn++]*tmp;
//			c[i,j]+=a[i,k]*tmp;
        }
      }
    return c;}

  public static matrix operator^ (matrix a, matrix b){
    int ac=a.cols, ar=a.rows, br=b.rows;
    matrix c = new matrix(ac,br);
    for(int ir=0;ir<ac;ir++)
      for(int ic=0;ic<br;ic++)
      {
        double s=0;
        for(int irar=ir*ar, icbr=ic*br; irar<ir*ar+ar;)
//		for(int irar=ir*ar, icbr=ic*br, k=0;k<ar; k++)
          s+=a.data[irar++]*b.data[icbr++];
//			s+=a.data[k+ir*ar]*b.data[k+ic*br];
        c.data[ir+ic*ar]=s;
      }
    /*			
    {
  double s=0; int ar=a.rows;
  for(int k=0;k<ar;k++)
    s+=a[k,ir]*b[k,ic];
  c[ir,ic]=s;
  }
*/
    return c;
  }

  public static matrix operator* (double a, matrix b) {
    int bc = b.cols, br = b.rows;
    matrix c = b.copy();
    for(int ir = 0; ir < br; ir++){
      for(int ic = 0; ic < bc; ic++) {
        c[ir,ic] *= a;
      }
    }
    return c;
  }

  /*		
public int rows(){return data[0].size;}
public int cols(){return data.GetLength(0);}
*/

  public matrix(matrix b) {
    nrows=b.rows; ncols=b.cols;
    data = new double[nrows*ncols];
    System.Array.Copy(b.data,data,b.data.Length);
  }


  /*		
public vector this[int col]{
get{return data[col];}
set{data[col]=value;}
}
*/

  internal double ValueAt(int row, int col){
    return data[row+col*nrows];
  }
  internal void ValueAt(int row, int col, double value){
    data[row+col*nrows] = value;
  }

  public matrix copy(){
    matrix c = new matrix(this);
    return c;
  }

  public static vector operator*(matrix a, vector x){
    vector y = new vector(x.size);
    for(int i=0;i<x.size;i++){
      y[i]=0; for(int ac=a.cols,k=0;k<ac;k++) y[i]+=a[i,k]*x[k];
    }
    return y;
  }

  public static vector operator^ (matrix a, vector x){
    vector y = new vector(a.cols);
    for(int i=0;i<a.cols;i++){
      y[i]=0; for(int k=0;k<a.rows;k++) y[i]+=a[k,i]*x[k];
    }
    return y;
  }

  public matrix T{
    get{return this.transpose();}
    set{}}

  public matrix transpose(){
//	int ncols=this.cols, nrows=this.rows;
    matrix c = new matrix(ncols,nrows);
    for(int ir=0;ir<nrows;ir++)
      for(int ic=0;ic<ncols;ic++) c[ic,ir]=this[ir,ic];
    return c;
  }

  public double trace() {
    /**
     * @Precond: Matrix is square
     */
    double trace = 0;

    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        if (i == j) {
          trace += ValueAt(i,j);
        }
      }
    }

    return trace;
  }


  public void print(){print("");}
  public void print(string s){
    System.Console.WriteLine(s);
    for(int ir=0;ir<this.rows;ir++){
      for(int ic=0;ic<this.cols;ic++)
        System.Console.Write("{0:F5} ",this[ir,ic]);
      System.Console.WriteLine();
    }
  }

}
