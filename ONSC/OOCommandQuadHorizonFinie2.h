#ifndef OOCOMMANDQUAD_H
#define OOCOMMANDQUAD_H
#include "OO.h"
#include "../RAND/rand.h"

template<typename T>
class OOCommandQuadHorizonFinie : public OO<T>
{
    private :
    Mat<T> Q;    /*positive semi definite matrix : state difference  cost*/
    Mat<T> R;    /*positive semi definite matrix : command cost*/
    double dt;    /*period between two state regenaration */
    double Tech;    /*period between two command generation */
    double horizonT;    /*time within which the goal is to be reached */
    Mat<T> (*Ftrans)(Mat<T> x, Mat<T> u, double dt);
    Mat<T> dX;
    
    int nbrCommand;
    int nbrStatePerBOZ;
    int nbrState;
    T bound; 	/*set the boundary to the commands computed.*/
    bool BOZ;
    
    
    public :
    
    OOCommandQuadHorizonFinie()
    {
    
    }
    
    OOCommandQuadHorizonFinie(int nbrCommand, Mat<T> (*F)(Mat<T>,Mat<T>,double), Mat<T> dX, double dt, double horizonT, T bound, bool BOZ = false )
    {
        Ftrans = F;
        this->dt = dt;
        this->Tech = (double)horizonT/(double)nbrCommand;
        this->horizonT = horizonT;
        this->bound = bound;
        this->dX = dX;
        this->BOZ = BOZ;
                
        
        /*Quadratic cost Matrix*/
        Q = Mat<T>((T)0,5,5);
        for(int i=Q.getLine();i--;)   Q.set((T)10,i+1,i+1);
        Q.set((T)10000,3,3);
        Q.set((T)0,4,4);
        Q.set((T)0,5,5);
        
        R = Mat<T>((T)0,2,2);
        for(int i=R.getLine();i--;)   R.set((T)1,i+1,i+1);
            
        this->nbrCommand = nbrCommand;
        this->nbrStatePerBOZ = (int) (Tech/dt);
        //nbrState = (int) (horizonT/dt);
        this->nbrState = nbrCommand*nbrStatePerBOZ;
        
    }
    
    ~OOCommandQuadHorizonFinie()
    {
    
    }
    
    Mat<T> energy(Mat<T> X)
    {              
        Mat<T>* state = new Mat<T>[nbrCommand*nbrStatePerBOZ];
        // #state = nbrCommand * nbrStateperBOZ
        Mat<T>* u = new Mat<T>[nbrCommand];
        // #u = nbrCommand
        
        for(int i=0;i<nbrCommand;i++)
        {            
            if( fabs_(X.get(i*2,1)) > bound)
            {
            	X.set( (X.get(i*2+1,1)>=(T)0 ? (T)1 : (T)-1)*bound,i*2,1);
            }   
            
            if( fabs_(X.get(i*2+1,1)) > bound)
            {
            	X.set( (X.get(i*2+1,1)>=(T)0 ? (T)1 : (T)-1)*bound,i*2+1,1);
            }
            
            u[i] = extract(X, i*2,1, i*2+1, 1);
        }
        
        state.insert( state.begin(), Mat<T>((T)0,5,1) );
        Mat<T> stepCommand(u[0]);
        Mat<T> offset((T)0*stepCommand);
        
        for(int i=0;i<nbrCommand;i++)
        {
            if(!BOZ)
            {
            	if(i<nbrCommand-1)
	            	stepCommand = u[i+1]-u[i];
	        else
	        	stepCommand = ((T)-1)*u[i];
	        	
	        stepCommand = (T)(1.0/nbrStatePerBOZ)*stepCommand;
	        offset = (T)0*stepCommand;
            }
            
            for(int j=0;j<nbrStatePerBOZ;j++)
            {            
            	if(BOZ)
	                state.insert(state.end(), Ftrans( state[i*nbrStatePerBOZ+j], u[i], dt) );
	        else
	        {
	        	state.insert(state.end(), Ftrans( state[i*nbrStatePerBOZ+j], u[i]+offset, dt) );
	        	offset = offset + stepCommand;
	        }
               
            }
        }
        
        Mat<T> costState(transpose(state[0])*(Q*state[0]) );
        for(int i=1;i<nbrState;i++)    costState = costState + transpose( state[i]-dX)*(Q * (state[i]-dX) );
        
        Mat<T> costCommand(transpose(u[0])*(R*u[0]) );
        for(int i=1;i<nbrCommand;i++)   costCommand = costCommand + transpose( u[i])*(R*u[i]);
        
        return (T)(1.0/2)*(costCommand+costState);
    }
    
    Mat<T> init()
    {
        //return Mat<T>((T)numeric_limits<T>::epsilon(),nbrCommand*2,1);
        //return Mat<T>((T)3,nbrCommand*2,1);
        //return sqrt((T)numeric_limits<T>::epsilon())*Mat<T>(nbrCommand*2,1,(char)1);
        //return sqrt((T)bound/2)*Mat<T>(nbrCommand*2,1,(char)1);
        
        Mat<T> rim(nbrCommand*2,1);
        NormalRand gen(bound/2,bound,(long)100);
        
        for(int i=nbrCommand*2;i--;)	rim.set( gen.dev(), i+1,1);
        
        return rim;
    }
    
    
    vector<Mat<T> > getCommands(Mat<T> X)
    {
    	vector<Mat<T> > u;
    	for(int i=0;i<nbrCommand;i++)
        {            
            u.insert(u.end(), extract(X, i*2,1, i*2+1, 1) );
        }
        
    	vector<Mat<T> > commands;
    	Mat<T> stepCommand(u[0]);
        Mat<T> offset((T)0*stepCommand);
    	for(int i=0;i<nbrCommand;i++)
        {
            if(!BOZ)
            {
            	if(i<nbrCommand-1)
	            	stepCommand = u[i+1]-u[i];
	        else
	        	stepCommand = ((T)-1)*u[i];
	        	
	        stepCommand = (T)(1.0/nbrStatePerBOZ)*stepCommand;
	        offset = (T)0*stepCommand;
            }
            
            for(int j=0;j<nbrStatePerBOZ;j++)
            {
            	if(BOZ)
	                commands.insert(commands.end(), u[i] );
	        else
	        {
	        	commands.insert(commands.end(), u[i]+offset );
	        	offset = offset + stepCommand;
	        }
            }
        }
        
    	return commands;
    }
};

#endif
