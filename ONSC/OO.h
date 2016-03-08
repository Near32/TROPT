#ifndef OO_H
#define OO_H

#include "../MATv2/Mat.h"
#include "../MVG/Frame.h"
#include <vector>

//#define debugOptim
template<typename T>
Mat<T> pyramidSubsampling(const Mat<T>& idMap, const Mat<T>& varIdMap)
{
	Mat<T> ridMap((T)0,idMap.getLine()/2,idMap.getColumn()/2);
	
	for(int i=ridMap.getLine();i--;)
	{	
		for(int j=ridMap.getColumn();j--;)
		{
			//T varMean = (T)0;
			T coeffMean = (T)0;
			T id = (T)0;
			int usedValue = 0;
			
			for(int ii=2;ii--;)
			{
				for(int jj=2;jj--;)
				{
					T currentId = idMap.get(2*i+ii+1,2*j+jj+1); 
					
					if(currentId != (T)0)
					{
						//if there is an actual value, then we incorporate
						usedValue++;
						
						T currentVar = varIdMap.get(2*i+ii+1,2*j+jj+1);
						if(currentVar == (T)0)
							currentVar = numeric_limits<T>::epsilon();
						
						coeffMean += (T)(1.0f/currentVar);
						
						id += (1.0f/currentVar)*currentId;
						//varMean += currentVar;
					}
						
				}
			}
			
			if(usedValue)
			{
				ridMap.set( id/coeffMean, i+1,j+1);
				//if there was an actual pixel used here, then we incorporate the value, else we set 0 :
			}
			else
			{
				ridMap.set( (T)0, i+1,j+1);
			}
		}
	}
	
	return ridMap;
}

template<typename T>
class OO
{
	protected :	
	
	public :

    vector<Mat<T> > params;
	
    OO()
    {

    }

	OO(vector<Mat<T> > params)
	{
		 this->params = params;
	
	}
	
	~OO()
	{
	
	}
	
    virtual Mat<T> energy(Mat<T> X) = 0;
	
	virtual Mat<T> init() = 0;

};


template<typename T>
class OOpoly : public OO<T>
{
	public :
	
	bool verbosetest;
	
	OOpoly(vector<Mat<T> > params) : OO<T>(params)
	{
		this->params = params;
		verbosetest = false;
		/*@params :
		 * 0 : A
		 * 1 : B
		 * 2 : Kinitial...
		 **/
	
	}
	
	~OOpoly()
	{
	
	}
	
	/* returns the value of a polynome whose coefficient are given in params.
	 * @params : 
	 * X : variable vector
	 **/
    Mat<T> energy(Mat<T> X)
	{
	
		T value = (T)0;
        Mat<T> temp((T)1,X.getColumn(),1);
		for(int i=0;i<params[0].getLine();i++)
		{											
			value += ( Line(params[0],i+1)*temp ).get(1,1);
            temp = X%temp;
		}		
		
		return (T)value*Mat<T>((T)1,1,1);
	}
	
	/* returns the initial value of the matrix in a vector form
	 **/
	Mat<T> init()
	{		
		return Mat<T>((T)200,1,1);//(int)params[0].getColumn(),1,(char)1);		
	}
	
	private :
	
	vector<Mat<T> > params;

};



template<typename T>
class OODirectSlam : public OO<T>
{

	private :
	
    Frame<T>* I1;
    Frame<T>* I2;
	Mat<T>* K;
	Mat<T>* invK;
    Mat<T>* depthMapI1;
    Mat<T>* depthMapI2;
	
	public :
	
    OODirectSlam( Frame<T>* I1, Frame<T>* I2, Mat<T>* K, Mat<T>* invK, Mat<T>* depthMapI1) : OO<T>()
	{
	
        this->I1 = I1;
		this->I2 = I2;
		this->K = K;
		this->invK = invK;
        this->depthMapI1 = depthMapI1;
        this->depthMapI2 = depthMapI1;
		
	}
	
	~OODirectSlam()
	{
		
	}
	
    Mat<T> energy(Mat<T> X)
	{

        Mat<T> pos1((T)1,3,1);
        Mat<T> pos2((T)0,3,1);
        Mat<T> temp(I1->get(pos1));
        Mat<T> rot(expW( extract(X,1,1,3,1)) );
        Mat<T> t(extract(X,4,1,6,1) );
        //Mat<T> krot((*K)*rot);
		
        for(int i=depthMapI1->getLine();i--;)
		{
			pos1.set(i+1,1,1);
			
            for(int j=depthMapI1->getColumn();j--;)
			{                
                if( depthMapI1->get(i+1,j+1) > (T)0)
				//only compute for the pixels that are in front obviously...
				{
					pos1.set(j+1,2,1);

                    //pos1.afficher();
					pos2 =(*invK)*pos1;                    
                    //pos2.afficher();                    
                    pos2 = (1.0/depthMapI1->get(i+1,j+1))*pos2;
                    //pos2.afficher();                    
                    pos2 = (*K) *(rot*pos2+t);
                    //pos2.afficher();

                    temp = temp + ( I2->get( pos2) - I1->get(pos1));
				}
			}
		}
        T norm = norme2(temp);
        return norm*norm*Mat<T>((T)1,1,1);
	}
	
    Mat<T> init()
	{
        return (T)(1.0/100)*Mat<T>(6,1,(char)1);//(T)numeric_limits<T>::epsilon(),6,1);
	}
	
};

template<typename T>
class OODirectSlamSIM3 : public OO<T>
{

    private :

    Frame<T>* I1;
    Frame<T>* I2;
    Mat<T>* K;
    Mat<T>* invK;
    Mat<T>* depthMapI1;
    Mat<T>* depthMapI2;
    vector<Frame<T> > pyramidI1;
    vector<Frame<T> > pyramidI2;
    vector<Mat<T> > pyramidD1;
    //vector<Mat<T> > pyramidD2;



    public :

    int pyramidDepth;
    int h;
    int w;

    OODirectSlamSIM3( Frame<T>* I1, Frame<T>* I2, Mat<T>* K, Mat<T>* invK, Mat<T>* depthMapI1) : OO<T>()
    {

        this->I1 = I1;
        this->I2 = I2;
        this->K = K;
        this->invK = invK;
        this->depthMapI1 = depthMapI1;
        this->depthMapI2 = depthMapI1;

        computePyramid();
        pyramidDepth = 0;
        this->params.insert(this->params.begin(), Mat<T>((T)0,1,1));

    }

    ~OODirectSlamSIM3()
    {

    }

    Mat<T> energy(Mat<T> X)
    {
        pyramidDepth = this->params[0].get(1,1);
        Mat<T> pos1((T)1,3,1);
        Mat<T> pos2((T)0,3,1);
        Mat<T> temp(I1->get(pos1));
        T s = X.get(7,1);
        Mat<T> rot(s*expW( extract(X,1,1,3,1) ) );
        Mat<T> t(extract(X,4,1,6,1) );
        Mat<T> inc((T)1,1,1);

        //Mat<T> krot((*K)*rot);
        h = pyramidD1[pyramidDepth].getLine();
        w = pyramidD1[pyramidDepth].getColumn();

        for(int i=h;i--;)
        {
            pos1.set(i+1,1,1);

            for(int j=w;j--;)
            {
                if( pyramidD1[pyramidDepth].get(i+1,j+1) > (T)0)
                //only compute for the pixels that are in front obviously...
                {
                    pos1.set(j+1,2,1);

                    //pos1.afficher();
                    pos2 =(*invK)*pos1;
                    //Here is hC = ( hx hy hz) = RCW * hW
                    pos2 = rot*pos2+t;
                    //rot = RWC --> pos2 = hW + rWC = m(theta,phi)_unnormalized...
                    //normalization
                    T norme = norme2(pos2);
                    pos2 = (T)(1.0/norme)*pos2;
                    //make it have the correct distance :
                    pos2 = (T)(1.0/pyramidD1[pyramidDepth].get(i+1,j+1))*pos2;
                    //Here was Xfeatures in the camera frame of the second image
                    //now we just have to project it again :
                    pos2 = (*K)*pos2;
                    homogeneousNormalization(&pos2);
                    inc = ( pyramidI2[pyramidDepth].get(pos2) - pyramidI1[pyramidDepth].get(pos1));
                    temp = temp + inc%inc;
                }
            }
        }
        T norm = norme2(temp);
        return norm*Mat<T>((T)1,1,1);
    }

    Mat<T> init()
    {
        //return (T)(1.0/100)*Mat<T>(7,1,(char)1);//
        //return Mat<T>((T)numeric_limits<T>::epsilon(),7,1);
        Mat<T> r(7,1,(char)1);
        r = (T)numeric_limits<T>::epsilon()*r;
        //return Mat<T>((T)1,7,1);

        //Mat<T> r((T)0,7,1);
        r.set((T)1,7,1);
        return r;
    }

    void computePyramid(int nbrL = 4)
    {
        Mat<T> kernel((T)1.0/4,2,2);
        pyramidI1.insert(pyramidI1.begin(), *I1);
        pyramidI2.insert(pyramidI2.begin(), *I2);
        pyramidD1.insert(pyramidD1.begin(), *depthMapI1);

        for(int i=1;i<=nbrL-1;i++)
        {
            pyramidI1.insert(pyramidI1.begin(), Frame<T>( pooling(pyramidI1[0],kernel), I1->getChannel()) );
            pyramidI2.insert(pyramidI2.begin(), Frame<T>( pooling(pyramidI2[0],kernel), I2->getChannel()) );
            pyramidD1.insert(pyramidD1.begin(), pooling(pyramidD1[0], kernel) );

            //afficherMat(&pyramidI1[0],&pyramidI2[0],I1,false,5.0);
        }

    }

};

template<typename T>
class OODirectSlamSE3 : public OO<T>
{

    protected :

    Frame<T>* I1;
    Frame<T>* I2;
    Frame<T>* gradI2X;
    Frame<T>* gradI2Y;
    Mat<T>* K;
    Mat<T>* invK;
    Mat<T>* depthMapI1;
    Mat<T>* depthMapI2;
    Mat<T>* var_depthMapI1;
    Frame<T>* pyramidI1;
    Frame<T>* pyramidI2;
    Frame<T>* pyramidGradI2X;
    Frame<T>* pyramidGradI2Y;
    bool grad;
    Mat<T>* pyramidD1;
    Mat<T>* pyramidK;
    Mat<T>* pyramidKinv;
    Mat<T>* pyramidvar_depthMapI1;
    //vector<Mat<T> > pyramidD2;



    public :

    int pyramidDepth;
    int h;
    int w;

    OODirectSlamSE3( Frame<T>* I1, Frame<T>* I2, Mat<T>* K, Mat<T>* invK, Mat<T>* depthMapI1, Mat<T>* var_depthMapI1, int pyramidDepth = 4) : OO<T>()
    {
	grad = false;
        this->I1 = I1;
        this->I2 = I2;
        this->gradI2X = NULL;
        this->gradI2Y = NULL;
        this->K = K;
        this->invK = invK;
        this->depthMapI1 = depthMapI1;
        this->depthMapI2 = depthMapI1;
        this->var_depthMapI1 = var_depthMapI1;
        this->pyramidDepth = pyramidDepth;

        computePyramid(pyramidDepth);
        pyramidDepth = 0;
        this->params.insert(this->params.begin(), Mat<T>((T)0,1,1));

    }
    
    OODirectSlamSE3( Frame<T>* I1, Frame<T>* I2, Frame<T>* gradI2X, Frame<T>* gradI2Y, Mat<T>* K, Mat<T>* invK, Mat<T>* depthMapI1, Mat<T>* var_depthMapI1, int pyramidDepth = 4) : OO<T>()
    {
    	
    	

	grad = true;
        this->I1 = I1;
        this->I2 = I2;
        this->gradI2X = gradI2X;
        this->gradI2Y = gradI2Y;
        this->K = K;
        this->invK = invK;
        this->depthMapI1 = depthMapI1;
        this->depthMapI2 = depthMapI1;
        this->var_depthMapI1 = var_depthMapI1;
        this->pyramidDepth = pyramidDepth;
#ifdef debugOptim
	clock_t timepyramid = clock();
#endif
	
        computePyramid(pyramidDepth);
        pyramidDepth = 0;
        this->params.insert(this->params.begin(), Mat<T>((T)0,1,1));
        
#ifdef debugOptim	
    	cout << "LA CREATION DE LA PYRAMIDE A PRISE : " << ((float)(clock()-timepyramid)/CLOCKS_PER_SEC) << " secondes." << endl;
#endif  
  	
    }

    ~OODirectSlamSE3()
    {
	delete[] pyramidI1;
	delete[] pyramidI2;
	if(grad)
	{
		delete[] pyramidGradI2X;
		delete[] pyramidGradI2Y;
	}
	delete[] pyramidD1;
	delete[] pyramidK;
	delete[] pyramidKinv;
	delete[] pyramidvar_depthMapI1;
    }

    Mat<T> energy(Mat<T> X)
    {
        pyramidDepth = this->params[0].get(1,1);
        Mat<T> pos1((T)1,3,1);
        Mat<T> pos2((T)0,3,1);
        Mat<T> temp(I1->get(pos1));
        Mat<T> rot(expW( extract(X,1,1,3,1) ) );
        Mat<T> t(extract(X,4,1,6,1) );
        Mat<T> inc((T)1,1,1);

        //Mat<T> krot((*K)*rot);
        h = pyramidD1[pyramidDepth].getLine();
        w = pyramidD1[pyramidDepth].getColumn();

        for(int i=h;i--;)
        {
            pos1.set(i+1,1,1);

            for(int j=w;j--;)
            {
                if( pyramidD1[pyramidDepth].get(i+1,j+1) > (T)0)
                //only compute for the pixels that are in front obviously...
                {
                    pos1.set(j+1,2,1);

                    //pos1.afficher();
                    //pos2 =(*invK)*pos1;
                    pos2 = pyramidKinv[pyramidDepth]*pos1;
                    //Here is hC = ( hx hy hz) = RCW * hW
                    pos2 = rot*pos2+t;
                    //rot = RWC --> pos2 = hW + rWC = m(theta,phi)_unnormalized...
                    
                    //normalization                    
                    /*
                    T norme = norme2(pos2);
                    pos2 = (T)((T)1.0/norme)*pos2;*/
                    
                    //make it have the correct distance :
                    pos2 = (T)((T)1.0/pyramidD1[pyramidDepth].get(i+1,j+1))*pos2;
                    //Here was Xfeatures in the camera frame of the second image
                    //now we just have to project it again :
                    //pos2 = (*K)*pos2;
                    pos2 = pyramidK[pyramidDepth]*pos2;
                    
                    homogeneousNormalization(&pos2);
                    
                    inc = ( pyramidI2[pyramidDepth].get(pos2) - pyramidI1[pyramidDepth].get(pos1));
                    temp = temp + inc%inc;
                }
            }
        }
        T norm = norme2(temp);
        return norm*Mat<T>((T)1,1,1);
    }

    Mat<T> init()
    {
        //return (T)(1.0/100)*Mat<T>(7,1,(char)1);//
        //return Mat<T>((T)numeric_limits<T>::epsilon(),7,1);
        Mat<T> r(6,1,(char)1);
        r = (T)numeric_limits<T>::epsilon()*r;
        //return Mat<T>((T)1,7,1);

        //Mat<T> r((T)0,7,1);
        r.set((T)1,7,1);
        return r;
    }

    void computePyramid(int nbrL = 4)
    {
#ifdef debugOptim
	clock_t time = clock();
#endif
        Mat<T> kernel((T)1.0,2,2);
        pyramidI1 = new Frame<T>[nbrL];
        pyramidI1[nbrL-1] = *I1;
        pyramidI2 = new Frame<T>[nbrL];
        pyramidI2[nbrL-1] = *I2;
        
        
        if(grad)
        {
		pyramidGradI2X = new Frame<T>[nbrL];
		pyramidGradI2X[nbrL-1] = *gradI2X;
		pyramidGradI2Y = new Frame<T>[nbrL];
		pyramidGradI2Y[nbrL-1] = *gradI2Y;
	}
	else
	{
		pyramidGradI2X = NULL;
		pyramidGradI2Y = NULL;
	}
	
        pyramidD1 = new Mat<T>[nbrL];
        pyramidD1[nbrL-1] = *depthMapI1;
        pyramidK = new Mat<T>[nbrL];
        pyramidK[nbrL-1] = *K;
        pyramidKinv = new Mat<T>[nbrL];
        pyramidKinv[nbrL-1] = *invK;
        pyramidvar_depthMapI1 = new Mat<T>[nbrL];
        pyramidvar_depthMapI1[nbrL-1] = *var_depthMapI1;
        /*
        pyramidI1.insert(pyramidI1.begin(), *I1);
        pyramidI2.insert(pyramidI2.begin(), *I2);
        pyramidD1.insert(pyramidD1.begin(), *depthMapI1);
        pyramidK.insert(pyramidK.begin(), *K);
        pyramidKinv.insert(pyramidKinv.begin(), *invK);
        pyramidvar_depthMapI1.insert(pyramidvar_depthMapI1.begin(), *var_depthMapI1);
        */
        
#ifdef debugOptim	
    	cout << "LA CREATION DE LA PYRAMIDE : INITIALISATION MALLOC :: " << ((float)(clock()-time)/CLOCKS_PER_SEC) << " secondes." << endl;
#endif 
        Mat<T> tempK(*K);
#ifdef debugOptim        
        clock_t cumulPool = 0;
#endif        

	if(grad)
	{
		for(int i=nbrL-1;i--;)
		{
#ifdef debugOptim        
			clock_t incPool = clock();
#endif 		
		    //pyramidI1[i] = Frame<T>( pooling( &pyramidI1[i+1], &kernel), I1->getChannel() );
		    pyramidI1[i] = Frame<T>( subsampled( &pyramidI1[i+1], &kernel), I1->getChannel() );
#ifdef debugOptim        
        		cumulPool += clock()-incPool;
#endif 		    
		    //pyramidGradI2X[i] = Frame<T>( pooling( &pyramidGradI2X[i+1], &kernel), gradI2X->getChannel() );
		    pyramidGradI2X[i] = Frame<T>( subsampled( &pyramidGradI2X[i+1], &kernel), gradI2X->getChannel() );
		    //pyramidGradI2Y[i] = Frame<T>( pooling( &pyramidGradI2Y[i+1], &kernel), gradI2Y->getChannel() );
		    pyramidGradI2Y[i] = Frame<T>( subsampled( &pyramidGradI2Y[i+1], &kernel), gradI2Y->getChannel() );
		    //pyramidI2[i] = Frame<T>( pooling( &pyramidI2[i+1], &kernel), I2->getChannel() );
		    pyramidI2[i] = Frame<T>( subsampled( &pyramidI2[i+1], &kernel), I2->getChannel() );
		    
		    //pyramidvar_depthMapI1[i] = pooling( &pyramidvar_depthMapI1[i+1], &kernel,(int)2 );
		    pyramidvar_depthMapI1[i] = subsampled( &pyramidvar_depthMapI1[i+1], &kernel);
		    
		    //pyramidD1[i] = pooling( &pyramidD1[i+1], &kernel, (int)2/*type : 3 : mean 2 :maxabs*/ );
		    //pyramidD1[i] = subsampled( &pyramidD1[i+1], &kernel);
		    pyramidD1[i] = pyramidSubsampling( pyramidD1[i+1], pyramidvar_depthMapI1[i+1]);
		    
		    //TODO : test les projections et retroprojections avec les K correspondant à la pyramide...?
		    tempK.set( (T)tempK.get(1,3)/2, 1,3);
		    tempK.set( (T)tempK.get(2,3)/2, 1,3);
		    tempK.set( (T)tempK.get(1,1)/2, 1,1);
		    tempK.set( (T)tempK.get(2,2)/2, 2,2);		    
		    pyramidK[i] = tempK;
		    pyramidKinv[i] = invGJ(tempK);

		    //afficherMat(&pyramidI1[i],&pyramidI2[i], &pyramidD1[i],false,1.0);
		}
	}
	else
	{
		for(int i=nbrL-1;i--;)
		{
		    pyramidI1[i] = Frame<T>( pooling(pyramidI1[i+1],kernel), I1->getChannel() );
		  
		    pyramidI2[i] = Frame<T>( pooling(pyramidI2[i+1],kernel), I2->getChannel() );
		    pyramidD1[i] = pooling(pyramidD1[i+1], kernel, (int)2/*type : 3 : mean 2 :maxabs*/ );
		    pyramidvar_depthMapI1[i] = pooling(pyramidvar_depthMapI1[i+1], kernel, (int)2 );
		    
		    //TODO : test les projections et retroprojections avec les K correspondant à la pyramide...?
		    tempK.set( (T)tempK.get(1,3)/2, 1,3);
		    tempK.set( (T)tempK.get(2,3)/2, 2,3);
		    tempK.set( (T)tempK.get(1,1)/2, 1,1);
		    tempK.set( (T)tempK.get(2,2)/2, 2,2);		    
		    pyramidK[i] = tempK;
		    
		    pyramidKinv[i] = invGJ(tempK);

		    //afficherMat(&pyramidI1[0],&pyramidI2[0], &pyramidD1[0],false,5.0);
		}
	}
	
#ifdef debugOptim	
    	cout << "LA CREATION DE LA PYRAMIDE : CUMUL POOL :: " << ((float)(cumulPool)/CLOCKS_PER_SEC) << " secondes." << endl;
#endif 

    }

};

#endif
