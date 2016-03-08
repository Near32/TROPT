#ifndef LSOODIRECTSLAMSE3_H
#define LSOODIRECTSLAMSE3_H
#include "LSOO.h"
#include "../RAND/rand.h"
#include "../MATv2/Mat.h"
#include "../MVG/Frame.h"
#include <vector>

template<typename T>
class LSOODirectSLAMSE3 : public LSOO<T>
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
    
    Mat<T> offset;
    int level;	//level of the pyramid scheme
    int nbrL;
    
    Mat<T>* Jacobian;
    Mat<T>* cost;
    Mat<T>* y;

    LSOODirectSLAMSE3( Frame<T>* I1, Frame<T>* I2, Mat<T>* K, Mat<T>* invK, Mat<T>* depthMapI1, Mat<T>* var_depthMapI1, int pyramidDepth = 4) : LSOO<T>()
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
        nbrL = pyramidDepth;

        computePyramid(pyramidDepth);
        this->params.insert(this->params.begin(), Mat<T>((T)0,1,1));
        
        //------------------
        //OFFSET :
        T offsetx = K->get(1,3);
	T offsety = K->get(2,3);
	offset = Mat<T>(3,1);
	offset.set( (T)offsetx, 1,1);
	offset.set( (T)offsety, 2,1);
	offset.set( (T)0, 3,1);
	
	//------------------
	
        //------------------
        Jacobian = NULL;
        cost = NULL;
        y = NULL;

    }
    
    LSOODirectSLAMSE3( Frame<T>* I1, Frame<T>* I2, Frame<T>* gradI2X, Frame<T>* gradI2Y, Mat<T>* K, Mat<T>* invK, Mat<T>* depthMapI1, Mat<T>* var_depthMapI1, int pyramidDepth = 4) : LSOO<T>()
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

        computePyramid(pyramidDepth);
        pyramidDepth = 0;
        this->params.insert(this->params.begin(), Mat<T>((T)0,1,1));
        
        //------------------
        //OFFSET :
        T offsetx = K->get(1,3);
	T offsety = K->get(2,3);
	offset = Mat<T>(3,1);
	offset.set( (T)offsetx, 1,1);
	offset.set( (T)offsety, 2,1);
	offset.set( (T)0, 3,1);
	
        //------------------
        Jacobian = NULL;
        cost = NULL;
        y = NULL;

    }

    ~LSOODirectSLAMSE3()
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
	
	//------------------
        //
        //------------------
        if(Jacobian != NULL)
        	delete Jacobian;
        if(cost != NULL)
        	delete cost;
        if(y != NULL)
        	delete y;
    }

    Mat<T> energy(Mat<T> X)
    {
    	//INCREMENTAL UPDATE 
	//compute error :

	error = (T)0;
	Mat<T> xi((T)1,3,1);
	Mat<T> xioffset(3,1);
	Mat<T> tempsol(sol+x);
	Mat<T> warpedxi(3,1);

	for(int ii= this->pyramidD1[level].getLine();ii--;)
	{
		xi.set((T)ii+1, 1,1);

		for(int jj= this->pyramidD1[level].getColumn();jj--;)
		{
			xi.set((T)jj+1, 2,1);					

			if( (this->pyramidD1[level]).get( ii+1, jj+1) != (T)0)
			{										

				xioffset = xi+offset;    
				xioffset = retroproj(&xioffset);						
				//Forward Additive Algorithm :
				warpedxi = rigidMotion( &tempsol, &xioffset);

				T val = (T)1.0/this->pyramidvar_depthMapI1[level].get( ii+1, jj+1);                     
				//T val = (T)1.0;
				//val *= val;	                                 
				//----- RESIDUAL -----
				//Forward Additive Algorithm :
				T residual = -(this->pyramidI1[level]).mat[ii][jj] + (this->pyramidI2[level]).get( proj(&warpedxi) ).get(1,1);
				val *= pow(residual,2);
				//HUBER :
			
				T huberK = 1.345e-2;

				if(fabs_(val) <= (T)huberK)
				{
					val = (T)val*val/(2*huberK);//*this->pyramidvar_depthMapI1[level].get(i+1,j+1);
				}
				else
				    	val = fabs_(val)- huberK/2;	    
			

				//HUBER (kerl):
				/*
				T huberK = 1000.345;
				if(fabs_(val) < (T)huberK)
				{
					val = (T)1;
				}
				else
					val = huberK/fabs_(val);
				*/

				//TUKEY :
				/*
				T tukeyB = 4.685;
				if( fabs_(val) < tukeyB)
				{
					countVAR2SMALL++;
					val = pow( (1- pow(val/tukeyB, (float)2)), (float)2);
				}
				else
					val = 0;
				*/

				error += val;
			}
		}
	}
	error /= (this->pyramidD1[level].getColumn())*(this->pyramidD1[level].getLine());
        
    }

    Mat<T> init()
    {
        //Mat<T> r(6,1,(char)1);
        //r = (T)numeric_limits<T>::epsilon()*r;
        Mat<T> r((T)0,6,1);
        return r;
    }

    void computePyramid(int nbrL = 4)
    {
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
        Mat<T> tempK(*K);
        //clock_t timer = clock();

	if(grad)
	{
		for(int i=nbrL-1;i--;)
		{
		    pyramidI1[i] = Frame<T>( pooling( &pyramidI1[i+1], &kernel), I1->getChannel() );
		    pyramidGradI2X[i] = Frame<T>( pooling( &pyramidGradI2X[i+1], &kernel), gradI2X->getChannel() );
		    pyramidGradI2Y[i] = Frame<T>( pooling( &pyramidGradI2Y[i+1], &kernel), gradI2Y->getChannel() );
		    pyramidI2[i] = Frame<T>( pooling( &pyramidI2[i+1], &kernel), I2->getChannel() );
		    pyramidD1[i] = pooling( &pyramidD1[i+1], &kernel);// (int)2/*type : 3 : mean 2 :maxabs*/ );
		    pyramidvar_depthMapI1[i] = pooling( &pyramidvar_depthMapI1[i+1], &kernel );
		    
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
		    pyramidvar_depthMapI1[i] = pooling(pyramidvar_depthMapI1[i+1], kernel );
		    
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
	
	//cout << " La construction des pyramides a prise : " << (float)(clock()-timer)/CLOCKS_PER_SEC << " secondes." << endl;

    }    
    
    
    //--------------------------------------
    //
    //		Helper Functions
    //
    //--------------------------------------
    
    inline Mat<T> r(Mat<T>* khi,Mat<T>* x)
    {        
        return (this->pyramidI1[level]).get(*x) - (this->pyramidI2[level]).get( warp(khi,x) );
        //1x1 in grayscale.
    	//watch out the level.
    }
    
    inline Mat<T> warp(Mat<T>* khi, Mat<T>* x)
    {
    	Mat<T> rx(retroproj(x));
    	rx = rigidMotion(khi, &rx);
    	
    	return proj( &rx );
    }
    
    inline Mat<T> proj(Mat<T>* P)
    {        
        Mat<T> ret( (this->pyramidK[level])* (*P) );        
        homogeneousNormalization(&ret);        
        ret -= ((T)1)*offset;
        
        return ret;
        
    }
    
    inline Mat<T> retroproj(Mat<T>* x)
    {
        T zx = Z(x);
        Mat<T> rim( zx*((this->pyramidKinv[level])*(*x)) );
        return rim;
    }
    
    inline T Z(Mat<T>* x)
    {        
        return     (T)((float)1.0/(float)(this->pyramidD1[level].get( x->get(1,1)-offset.get(1,1),x->get(2,1)-offset.get(2,1) )) );
    	
    	//watch out the level.
    }
    
    inline Mat<T> rigidMotion(Mat<T>* khi,Mat<T>* P)
    {
        Mat<T> khi_hat(expM(khi) );
    	Mat<T> rot(extract(&khi_hat, 1,1, 3,3) );
    	Mat<T> t( extract(&khi_hat, 1,4, 3,4) );
    	
    	return rot*(*P)+t;
    }
    
    inline Mat<T> rigidMotionHAT(Mat<T>* khi_hat,Mat<T>* P)
    {
    	Mat<T> rot(extract(khi_hat, 1,1, 3,3) );
    	Mat<T> t( extract(khi_hat, 1,4, 3,4) );
    	
    	return rot*(*P)+t;
    }
    
    //----------------------------------------
    //
    //		END Helper Functions
    //
    //----------------------------------------
    
    
    //--------------------------------------
    //
    //		Accessors
    //
    //--------------------------------------
    
    
    Mat<T> getJacobian()
    {
    	return Jacobian;
    }
    
    Mat<T> getCost()
    {
    	return cost;
    }
    
    Mat<T> getY()
    {
    	return y;
    }
    
};

#endif

