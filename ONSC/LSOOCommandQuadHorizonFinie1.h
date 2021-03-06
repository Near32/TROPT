#ifndef LSOOCOMMANDQUAD_H
#define LSOOCOMMANDQUAD_H
#include "LSOO.h"
#include "../RAND/rand.h"

template<typename T>
class LSOOCommandQuadHorizonFinie : public LSOO<T>
{
    private :
    Mat<T> Q;    /*positive semi definite matrix : state difference  cost*/
    Mat<T> R;    /*positive semi definite matrix : command cost*/
    double dt;    /*period between two state regenaration */
    double Tech;    /*period between two command generation */
    double horizonT;    /*time within which the goal is to be reached */
    
    Mat<T> (*Ftrans)(Mat<T> x, Mat<T> u, double dt);	/*give (A,B) stacked in line with respect to x and u.*/
    Mat<T> dX;
    
    int nbrCommand;
    int nbrStatePerBOZ;
    int nbrState;
    //T bound; 	/*set the boundary to the commands computed.*/
    bool BOZ;
    
    Mat<T> cost;
    Mat<T> Jacobian;
    Mat<T> y;
    Mat<T>* A;
    Mat<T>* B;
    Mat<T>* state;
    Mat<T>* u;
    
    
    public :
    
    LSOOCommandQuadHorizonFinie()
    {
    
    }
    
    LSOOCommandQuadHorizonFinie(int nbrCommand, Mat<T> (*F)(Mat<T>,Mat<T>,double), Mat<T> dX, double dt, double horizonT, T bound, bool BOZ = false )
    {
        Ftrans = F;
        this->dt = dt;
        this->Tech = (double)horizonT/(double)nbrCommand;
        this->horizonT = horizonT;
        this->bound = bound;
        this->dX = dX;
        this->BOZ = BOZ;
                
        
        /*Quadratic cost Matrix*/
        Q = Mat<T>((T)0,3,3);
        for(int i=Q.getLine();i--;)   Q.set((T)1e-2,i+1,i+1);
        Q.set((T)1e-8,3,3);
        Q.set((T)0,4,4);
        Q.set((T)0,5,5);
        
        R = Mat<T>((T)0,2,2);
        for(int i=R.getLine();i--;)   R.set((T)1e-1,i+1,i+1);
            
        this->nbrCommand = nbrCommand;
        this->nbrStatePerBOZ = (int) (Tech/dt);
        //nbrState = (int) (horizonT/dt);
        this->nbrState = nbrCommand*nbrStatePerBOZ;
        
        Jacobian = Mat<T>(1, nbrCommand*R.getLine());
        y = Mat<T>((T)0,1,1);
        state = NULL;
        u = NULL;
        A = NULL;
        B = NULL;
        
    }
    
    ~LSOOCommandQuadHorizonFinie()
    {
    	if(u != NULL)
    		delete[] u;
    	if(state != NULL)
    		delete[] state;
    	if(A != NULL)
    		delete[] A;
    	if(B != NULL)
    		delete[] B;
    }
    
    Mat<T> energy(Mat<T> X)
    {   
    	if(u != NULL)
    		delete[] u;
    		           
        A = new Mat<T>[nbrState];
        B = new Mat<T>[nbrState];
        state = new Mat<T>[nbrState+1];
        u = new Mat<T>[nbrState];
        
        y = Mat<T>((T)0,1,1);
             
        //let us construct the commands :   
        for(int i=0;i<nbrCommand;i++)
        {            
            if( fabs_(X.get(i*2,1)) > this->bound)
            {
            	X.set( (X.get(i*2+1,1)>=(T)0 ? (T)1 : (T)-1)*this->bound,i*2,1);
            }   
            
            if( fabs_(X.get(i*2+1,1)) > this->bound)
            {
            	X.set( (X.get(i*2+1,1)>=(T)0 ? (T)1 : (T)-1)*this->bound,i*2+1,1);
            }
            
            u[i*nbrStatePerBOZ] = extract(X, i*2+1,1, i*2+2, 1);
            
            for(int j=1;j<nbrStatePerBOZ;j++)	u[i*nbrStatePerBOZ+j] = u[i*nbrStatePerBOZ+j-1];            
        }
        
        state[0] = Mat<T>((T)0,3,1);
        Mat<T> stepCommand(u[0]);
        Mat<T> offset((T)0, stepCommand.getLine(), 1);
        
        //let us construct the state of the computed trajectory :
        for(int i=0;i<nbrCommand;i++)
        {
            if(!BOZ)
            {
            	if(i<nbrCommand-1)
	            	stepCommand = u[(i+1)*nbrStatePerBOZ]-u[i*nbrStatePerBOZ];
	        else
	        	stepCommand = ((T)-1)*u[i*nbrStatePerBOZ];
	        	
	        stepCommand = (T)(1.0/nbrStatePerBOZ)*stepCommand;
	        offset = Mat<T>((T)0, stepCommand.getLine(), 1);
            }
            
            for(int j=0;j<nbrStatePerBOZ;j++)
            {            
            	if(BOZ)
            	{
            		Mat<T> temp( Ftrans( state[i*nbrStatePerBOZ+j], u[i], dt) );
            		A[i*nbrStatePerBOZ+j] = extract(&temp, 1,1, Q.getLine(),Q.getColumn());
            		B[i*nbrStatePerBOZ+j] = extract(&temp, 1, Q.getColumn()+1, temp.getLine(), temp.getColumn() );
	                state[i*nbrStatePerBOZ+j+1] =  A[i*nbrStatePerBOZ+j]*state[i*nbrStatePerBOZ+j]+B[i*nbrStatePerBOZ+j]*u[i*nbrStatePerBOZ+j];
	        }
	        else
	        {
	        	Mat<T> temp( Ftrans( state[i*nbrStatePerBOZ+j], u[i], dt) );
            		A[i*nbrStatePerBOZ+j] = extract(&temp, 1,1, Q.getLine(),Q.getColumn());
            		B[i*nbrStatePerBOZ+j] = extract(&temp, 1, Q.getColumn()+1, temp.getLine(), temp.getColumn() );
            		
            		u[i*nbrStatePerBOZ+j] += offset;
	                
	                state[i*nbrStatePerBOZ+j+1] =  A[i*nbrStatePerBOZ+j]*state[i*nbrStatePerBOZ+j]+B[i*nbrStatePerBOZ+j]*u[i*nbrStatePerBOZ+j];
	        	
	        	offset += stepCommand;
	        }
               
            }
        }
        
        
        Mat<T> costState(transpose(state[0])*(Q*state[0]) );
        y += costState;
        for(int i=1;i<nbrState;i++)
        {
        	/*----------------------------*/        
		/*Quadratic cost Matrix*/
		/*----------------------------*/
		/*
		for(int j=Q.getLine();j--;)   Q.set((T)10,j+1,j+1);
		Q.set((T)100*((T)1.0/((T)i/100)),1,1);
		Q.set((T)100*((T)1.0/((T)i/100)),2,2);
		T alpha = atan2( dX.get(2,1), dX.get(1,1) );
		
		
		Q.set((T)(alpha-state[i].get(3,1))*(alpha-state[i].get(3,1)) / ( (state[i].get(3,1)-dX.get(3,1))*(state[i].get(3,1)-dX.get(3,1))+sqrt(numeric_limits<T>::epsilon())), 3,3);
		Q.set((T)0,4,4);
		Q.set((T)0,5,5);
		*/
		/*----------------------------*/
	        
	        Mat<T> temp(transpose( dX-state[i])*(Q * (dX -state[i]) ));
	        costState += temp;
	        
	        if(i%nbrStatePerBOZ == 0)
	        	y += temp;
	}
        
        Mat<T> costCommand(transpose(u[0])*(R*u[0]) );
        y += costCommand;
        for(int i=1;i<nbrCommand*nbrStatePerBOZ;i++)   
        {
        	Mat<T> temp(transpose( u[i])*(R*u[i]) );
        	costCommand += temp;
        	
        	if(i%nbrStatePerBOZ == 0)
	        	y += temp;
	}
        	
        
        cost = (T)(1.0/2)*(costCommand+costState);
        
        
        //-----Computation of the Jacobian---------------
        Mat<T> I((T)0,R.getLine(),R.getLine());
        for(int i=I.getLine();i--;)	I.set((T)1,i+1,i+1);
        
        Mat<T> tempxi;
        Mat<T> tempui;
        //reinitialization :
        costState = Mat<T>((T)0, 1, R.getLine()*nbrCommand/*watchout : only command...*/);
        costCommand = Mat<T>((T)0, 1, R.getLine()*nbrCommand/*watchout : only command...*/);
        
        //computed only for the state,command that are part of the vector being optimized, not for the BOZ-or not states...
        // but still, many states are related to the optimized commands... :
        for(int i=0;i<nbrState;i++)
        {
        
        	/*----------------------------*/        
		/*Quadratic cost Matrix*/
		/*----------------------------*/
		/*
		for(int i=Q.getLine();i--;)   Q.set((T)10,i+1,i+1);
		Q.set((T)100*((T)1.0/((T)j/100)),1,1);
		Q.set((T)100*((T)1.0/((T)j/100)),2,2);
		T alpha = atan2( dX.get(2,1), dX.get(1,1) );
		
		
		Q.set((T)(alpha-state[j].get(3,1))*(alpha-state[j].get(3,1)) / ( (state[j].get(3,1)-dX.get(3,1))*(state[j].get(3,1)-dX.get(3,1))+sqrt(numeric_limits<T>::epsilon())), 3,3);
		Q.set((T)0,4,4);
		Q.set((T)0,5,5);
		*/
		/*----------------------------*/
		
        	tempxi = transpose(dX - state[i])*( ((T)-1)*(transpose(Q)+Q) );
        	Mat<T> tempA((T)0,state[0].getLine(),state[0].getLine());
	        for(int j=tempA.getLine();j--;)	tempA.set((T)1,j+1,j+1);	        	       	     
	        
	        Mat<T> tempSum((T)0,tempA.getLine(), /*watchout only that : */nbrCommand*R.getLine() );
        	for(int kk=0;kk<=i;kk++)
        	{
        		Mat<T> tempDeriv((T)0, R.getLine(), tempSum.getColumn());
        		for(int j=R.getLine();j--;)	tempDeriv.set( ( (T)(BOZ ? 1 : 1+(floor((kk+1)/nbrStatePerBOZ)+1-(kk+1)) ) ), j+1, 2*floor((kk+1)/nbrStatePerBOZ)+j-1);        		
        		
        		for(int k=i-kk;k--;)	tempA *= A[k];        		        		
        		
        		tempSum += tempA*(B[kk]*tempDeriv);
        		
        	}
        	 
        	costState += tempxi*tempSum;

        }               
        
        
        //computed only for the state,command that are part of the vector being optimized, not for the BOZ-or not states...
        // the same as before, many commandes are related to the optimized commands... :
        for(int i=0;i<nbrState;i++)
        {
        	tempui = transpose(u[i])*( transpose(R)+R );
        	
	        Mat<T> tempDeriv((T)0, R.getLine(), /*watchout : */nbrCommand*R.getLine());
       		for(int j=R.getLine();j--;)	tempDeriv.set( ( (T)(BOZ ? 1 : 1+(floor((i+1)/nbrStatePerBOZ)+1-i) ) ), j+1, 2*floor((i+1)/nbrStatePerBOZ)+j-1);       	       		
        	 
        	costCommand += tempui*tempDeriv;
        }
        
        Jacobian = costState+costCommand;
        //----------------------------------------------
        
        delete[] state;
        delete[] A;
        delete[] B;
        delete[] u;
        
        state = NULL;
        A = NULL;
        B = NULL;
        u = NULL;
        
        return cost;
    }
    
    Mat<T> init()
    {
        //return Mat<T>((T)numeric_limits<T>::epsilon(),nbrCommand*2,1);
        //return Mat<T>((T)3,nbrCommand*2,1);
        //return sqrt((T)numeric_limits<T>::epsilon())*Mat<T>(nbrCommand*2,1,(char)1);
        //return sqrt((T)bound/2)*Mat<T>(nbrCommand*2,1,(char)1);
        
        Mat<T> rim(nbrCommand*2,1);
        NormalRand gen(0, this->bound,(long)10);
        
        for(int i=nbrCommand*2;i--;)	rim.set( gen.dev(), i+1,1);
        
        return rim;
    }
    
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
    
    
    Mat<T>* getCommands(Mat<T> X)
    {
    	u = new Mat<T>[nbrCommand];
    	for(int i=0;i<nbrCommand;i++)
	{            
	    u[i] = extract(X, i*2,1, i*2+1, 1);
	}
	
    	Mat<T>* commands = new Mat<T>[nbrCommand*nbrStatePerBOZ];
    	// have to be freed in the main program....!
    	
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
		        commands[i*nbrStatePerBOZ+j] = u[i];
		else
		{
			commands[i*nbrStatePerBOZ+j] = u[i]+offset;
			offset = offset + stepCommand;
		}
	    }
	}
	
	delete[] u;
	u = NULL;
	
    	return commands;
  
	  
    }
};

#endif
