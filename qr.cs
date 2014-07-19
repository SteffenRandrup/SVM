namespace Projekt {
	public class QRdecomposition{

		public matrix Q, R;

		// Giver problemer hvis en søjle består af faktorer af samme fx 3,6,9
		public QRdecomposition(matrix A){
			int m = A.cols;
			R = new matrix(m,m);
			Q = A.copy();
			for(int i=0;i<m;i++){
				R[i,i]=Q[i].norm();
				Q[i] /= R[i,i];
				for(int j=i+1;j<m;j++){
					R[i,j]=Q[i]^Q[j];
					Q[j]-=Q[i]*R[i,j];}}
		}

		public vector qrbak(vector b){ return runbak(Q^b); }

		public vector runbak(vector c){
			int m=R.cols;
			vector x = new vector(m);
			for(int i=m-1;i>=0;i--){
				double s=0;
				for(int k=i+1;k<m;k++) s+=R[i,k]*x[k];
				x[i]=(c[i]-s)/R[i,i]; }
			return x;
		}

		public matrix inverse(){
			return rinverse()*Q.T;
		}

		public matrix rinverse(){
			int m=R.cols;
			var Ri = new matrix(m,m);
			var b = new vector(m);
			for(int i=0;i<m;i++)
			{
				b[i]=1;
				Ri[i]=runbak(b);
				b[i]=0;
			}
			return Ri;
		}

	}
}
