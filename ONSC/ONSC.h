#ifndef ONSC_H
#define ONSC_H
#include "../MATv2/Mat.h"
//#define verbose
//#define verbose_mini
//#define verbose_optimstep
#include "OO.h"

template<typename T>
Mat<T> computeTransformation( Mat<T> I1, Mat<T> I2, Mat<T> depthMapI1, Mat<T> K );
template<typename T>
T static_lambda(void);

template<typename T>
class ONSC
{

    protected :

    int it;

    vector<Mat<T> > cost;
	T alpha;	//learning rate controlled...
	vector<Mat<T> > x_k;					//variables
	vector<T> mu_k;
	int k;	
    bool optima;
		
	//T (*lambda_k)(void);	//argmin f(x_k - lambda_k*direction_k)
	//cf optimStep.
	vector<Mat<T> > Delta;					//Error terms.
	vector<Mat<T> > grad;					//Jacobians to store.
	
	
	Mat<T> (*ptrF)(Mat<T>);		//ptr to the function being minimized
	OO<T>* instance;		//convenient object comprising of the objective function and its parameters
	
	
	T optimStep()
	{
        int itMu = 5;
        T mu = 1.0/(norme2(x_k[k])+((T)(rand()%1000)));
		while(mu == 0.0){	mu = norme2(x_k[k])/100+((T)(rand()%100)/100); }
		T varMu = 0;
		
		for(int i=0;i<itMu;i++)
		{		
			varMu = deltaMu(mu);
			if(fabs(varMu) <= sqrt(numeric_limits<T>::epsilon())*mu && varMu != 0.0)
				i = itMu;				
				
				
			if(mu + varMu <= (T)0)
			{				
				mu = norme2(x_k[k])/10000+((T)(rand()%100)/100);
			}
			else if(fabs_(varMu) > mu)
			{
				varMu = sqrt(numeric_limits<T>::epsilon())*mu*mu;
				mu = mu+varMu;
			}
			else
				mu = mu + varMu;
				
			
	
#ifdef verbose_optimstep			
			cout << "varMu = " << varMu << endl;	
			cout << " MU = " << mu << endl;			
#endif			
		}
	
		return mu;
	}
	
	T deltaMu(T mu)
	{
		T var = sqrt(numeric_limits<T>::epsilon())*mu;		
		T delta;	
		delta = (instance->energy(x_k[k]-(mu+var)*Delta[k]) - instance->energy(x_k[k]-(mu-var)*Delta[k])).get(1,1);
#ifdef verbose_optimstep		
		cout << "Delta = "<< delta << endl;
		cout << " h = " << var << endl;
#endif		
		
		return (T)delta/(2*var);
	}
	
	public:	
	
	/*Constructor
	 * @param OO object 
	 * ooInst ; object which contains the objective function ready to be used with its parameters. OUGHT TO BE IMPLEMENTED
	 * init; pointer to the function which initializes and handle the variable x_k (give dimension and so on...) OUGHT TO BE IMPLEMENTED
	 * @param iteration ; number of iteration to do ; TODO by default there will be a stopping criteria applied.
	 **/
    ONSC(  OO<T>* ooInst, int it = 1, bool optimstep = false)
	{		
        optima = false;
        this->it = it;
		this->instance = ooInst;
		bool OO = true;
		x_k.insert( x_k.end(), instance->init() );
		
		
		//alpha = 10e-3; //before the implementation of optimStep, this was the default value.
		alpha = (T)1;

#ifdef verbose
		cout << "ONSC Initialization : DONE." << endl;	
#endif
				
        cost.insert(cost.end(), (instance->energy(x_k[0])) );
		Mat<T> variation((T)0*cost[0]);
		
        for(k=0; (k<=it && !optima);k++)
		{
            		callback();
#ifdef verbose_mini		
			cout << "///////////////////////////\n ONSC : Running : iteration " << k << " COST : " << endl;
			cost[k].afficher();
#endif			
			grad.insert( grad.end(), computeGrad(OO) );			
#ifdef verbose						
			cout << "ONSC : gradient : computed : "<< endl;

			(grad[k]).afficher();
#endif			
			
			T factor = (T)(1);
			T normeGrad = norme2(grad[k]);
			
            if(normeGrad == (T)0)
			{
				/*minimum local...*/
				/*random gradient... trying something... luck...*/                
                grad[k] = numeric_limits<T>::epsilon()*Mat<T>(grad[k].getLine(), grad[k].getColumn(), (char)1);
                normeGrad = norme2(grad[k]);

                /*We stop here, it should mean that we have reached an optima*/
                //optima = true;
			}
						
	     		
            if(!optima)
            {
                factor = 1.0/normeGrad;
#ifdef verbose
                cout << "FACTOR = " << factor << endl;
#endif			
                Delta.insert( Delta.end(), factor * grad[k] );
                alpha = 10e-3;
                if(optimstep)
                {
                    alpha = optimStep();
	     		}
            
                mu_k.insert(mu_k.end(),alpha);
                Delta[k] = alpha*Delta[k];
#ifdef verbose									
                cout << "ONSC : Delta : computed : size vector : " << Delta.size()  << endl;
                (Delta[k]).afficher();
#endif			
									
                x_k.insert( x_k.end(), x_k[k] - Delta[k] );
#ifdef verbose_mini		
                cout <<  "ONSC : X_k+1 : updated : "<< endl;
                x_k[k+1].afficher();
#endif			
			
                /*----------------------------------------------------------*/
                //Handling of the parameter alpha :
                cost.insert( cost.end(), instance->energy(x_k[k+1]) );//- cost[k] );
                variation = variation + cost[k+1]-cost[k];
#ifdef verbose_var
                cout << "///////// VARIATION : " << endl;
                variation.afficher();
#endif		        
			  
#ifdef verbose
		        cout << "alpha =" << alpha << endl;
#endif		        
		        /*----------------------------------------------------------*/                

            }
            else
            {
                /*warping and getting out : */
                k--;
            }
		        
		}

	}
	
	
	~ONSC()	{};
	
	Mat<T> getX(int rank = -1)	const
	{
        return x_k[ (rank < 0 ? k : (rank <= k ? rank : k) )];
	}
	
	T getMu(int rank = -1)	const
	{
		return mu_k[ (rank < 0 ? k : (rank <= k ? rank : k) )];
	}
	
	Mat<T> computeGrad(bool OO = false)
	{
		int n = x_k[k].getLine();
        	T h = sqrt(numeric_limits<T>::epsilon())*norme2(x_k[k]);
		Mat<T> delta((T)0, n,1);
        	delta.set(h, 1,1);
        	Mat<T> grad( (!OO ? ptrF(x_k[k]+delta) - ptrF(x_k[k]-delta)  : instance->energy(x_k[k]+delta) - instance->energy(x_k[k]-delta)  )  );
		
		for(int i=2;i<=n;i++)
		{
			delta.set((T)0,i-1,1);
		    delta.set(h,i,1);
		    grad = operatorC(grad, (!OO ? ptrF(x_k[k]+delta) - ptrF(x_k[k]-delta)  : instance->energy(x_k[k]+delta) - instance->energy(x_k[k]-delta)  )  );
		}
		
		return ((T)(1.0/(2*h)))*grad;
	}


    void callback()
    {
        //DO NOTHING HERE...
        int l = 4;
        int pas = it/l;
        while(pas==0)
        {
            l--;
            pas = it/l;
        }

        if(instance->params.size() < 1)
            instance->params.insert(instance->params.begin(), Mat<T>((T)( (k/pas <= l-1 ? k/pas : l-1)),1,1) );
        else
            instance->params[0] = Mat<T>((T)( (k/pas <= l-1 ? k/pas : l-1) ),1,1);
    }

};



template<typename T>
Mat<T> computeTransformation( Mat<T> I1, Mat<T> I2, Mat<T> depthMapI1, Mat<T> K )
{
    int it = 5;
	OO<T> instanceOO(I1,I2, K, depthMapI1);
    ONSC<T> instance(instanceOO.energy, instanceOO.initX, it);
	
	return instance.getBeta();
	

}

template<typename T>
T static_lambda(void)
{
	return (T)10;
}

#endif
