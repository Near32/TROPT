#include "Mat.h"
#include "OO.h"

template<typename T>
Mat<T> computeTransformation( Mat<T> I1, Mat<T> I2, Mat<T> depthMapI1, Mat<T> K );
double static_lambda(void);

template<typename T>
class ONSC
{

	private :
	double alpha;	//learning rate controlled...
	vector<Mat<T> > x_k;					//variables
	int k;
	//vector<Mat<T> > params;				//argument to the function.	
		
	double (*lambda_k)(void);	//argmin f(x_k - lambda_k*direction_k)
	vector<Mat<T> > Delta;					//Error terms.
	vector<Mat<T> > grad;					//Jacobian to store.
	
	Mat<T> (*ptrF)(Mat<T>);		//ptr to the function. @param variable (@param params ? any ?)
	
	
	
	public:
	
	/*Constructor
	 * @param Ff ; pointer to the function to be optimized. OUGHT TO BE IMPLEMENTED
	 * @param init; pointer to the function which initializes and handle the variable x_k (give dimension and so on...) OUGHT TO BE IMPLEMENTED
	 * @param iteration ; number of iteration to do ; by default there will be a stopping criteria applied.
	 * @param lambda ; pointer to the function which handles the computation of the lambda_k at each step. OUGHT TO BE IMPLEMENTED.
	 **/
	ONSC(  Mat<T> (*Ff)(Mat<T>), /*vector<Mat<T> > params_,*/, Mat<T> (*init)(void), int iteration = 1, double (*lambda)(void) = static_lambda)
	{
		ptrF = Ff;
		x_k.insert( x_k.end(), init() );
		
		
		//params = params_;
		alpha = 10e7;
		this.lambda_k = lambda;				

		cout << "ONSC Initialization : DONE." << endl;	
		vector<Mat<double> > cost;
		cost.insert(cost.end(), (ptrF(x_k[0])) );
		
		for(k=0;k<=iteration;k++)
		{
			cout << "///////////////////////////\nGNA : Running : iteration " << k << endl;
			grad.insert( grad.end(), computeGrad() );			
			cout << "ONSC : gradient : computed : "<< endl;
#ifdef verbose			
			grad[k].afficher();
#endif			
			
			/*
			Mat<T> temp(invGJ(transpose(J)*J));
			//TODO : use SVD decomposition ; test if faster ?
			temp = temp*transpose(J);
			cout << "GNA : inverse of J : computed." << endl;
#ifdef verbose			
			temp.afficher();
#endif						
			
			Delta = temp * (realValue - ptrF(params, Beta));
			*/
			Delta.insert( Delta.end(), ((T)(alpha*lambda_k()))* ((T)1.0/norme(grad[k])) * grad[k] );
			cout << "ONSC : Delta : computed." << endl;
			
#ifdef verbose						
			Delta[k].afficher();
#endif			
						
			//Beta = Beta - ((T)(alpha*lambda_k()))*Delta;
			x_k.insert( x_k.end(), x_k - Delta[k] );
			cout <<  "ONSC : X_k+1 : updated : "<< endl;
			x_k[k+1].afficher();
			
			/*----------------------------------------------------------*/
			//Handling of the parameter alpha :
			cost.insert( cost.end(), ptrF(x_k[k+1]) - cost[k] );
	    	        cost[k+1].afficher();
		        if(!(cost[k+1] <= Mat<double>((double)0, cost[k+1].getLine(),cost[k+1].getColumn()) ) )
			    alpha = (alpha >= numeric_limits<double>::epsilon()*10e2 ? alpha/10 : alpha);
		        else
			    alpha = alpha*10;

		        cout << "alpha =" << alpha << endl;

		        cost[k+1] = ptrF(x_k[k+1]);
		        /*----------------------------------------------------------*/
		        
		}

	}
	
	
	~ONSC()	{};
	
	Mat<T> getX()	const
	{
		return x_k;
	}
	
	Mat<T> computeGrad()
	{
		int n = x_k[k].getLine();
		T h = numeric_limits<T>::epsilon();
		Mat<T> delta((T)0, n,1);
		delta.set(h, 1,1);
		Mat<T> grad(ptrF(x_k[k]+delta) - ptrF(x_k[k]-delta));
		
		for(int i=2;i<=n;i++)
		{
			delta.set((T)0,i-1,1);
			delta.set(h,i,1);
			grad = operatorL(grad, ptrF(x_k[k]+delta) - ptrF(x_k[k]-delta));		
		}
		
		return ((T)(1.0/(h*h)))*grad;
	}

};



template<typename T>
Mat<T> computeTransformation( Mat<T> I1, Mat<T> I2, Mat<T> depthMapI1, Mat<T> K )
{
	int iteration = 5;
	OO<T> instanceOO(I1,I2, K, depthMapI1);
	ONSC<T> instance(instanceOO.energy, instanceOO.initX, iteration);
	
	return instance.getBeta();
	

}


double static_lambda(void)
{
	return (double)1;
}

