#ifndef INVERSION_H
#define INVERSION_H

#include "Mat.h"

template<typename T>
Mat<T> gaussjordanInversion( Mat<T>* A_, Mat<T>* b_, int method);
template<typename T>
Mat<T> cramerInversion(const Mat<T>	& A);


template<typename T>
Mat<T> gaussjordanInversion( Mat<T>* A_, Mat<T>* b_, int method)
{
    if( A_ != NULL && b_ != NULL)
    {
        Mat<T> A(*A_), b(*b_);        

        switch(method)
        {
            case 1 :
            {
            /* Full-Pivoting */
            int n, m;
            n = A.getLine();
            m = b.getColumn();
            T big, pivinv;

            /*bookkeeping */
            int irow=1, icol=1;
            Mat<T> indxc(n, 1);
            Mat<T> indxr( n,1);
            Mat<T> ipiv((T)0, n,1);


            Mat<T> I( (T)0, n,n);
            for(int i=n;i--;)
                I.set((T)1, i+1,i+1);

            for(int i=1;i<=n;i++)
            {
                big = (T)0;
                for(int j=1;j<=n;j++)
                {
                    if(ipiv.get(j,1) != (T)1)	/*si cette ligne j ne contient pas deja un pivot... */
                    {
                        for(int k=1;k<=n;k++)
                        {
                            if(ipiv.get(k,1) == 0) /* si cette ligne k ne contient pas deja un pivot..*/
                            {
                                if(abs(A.get(j,k)) > big)
                                {
                                    big = abs(A.get(j,k));
                                    irow = j;
                                    icol = k;
                                }

                            }
                        }

                    }
                }
                ipiv.set(ipiv.get(icol,1)+1, icol,1);
                /*on a bien un pivot à la icol ligne, car c'est la valeur la plus grande de */

                /*Maintenant que l'on a le pivot, il faut le mettre sur la diagonale en interchangeant les colonnes
                 * --> b doit changer de la meme manière que A ; l'inverse obtenue devra etre changee.
                 */

                indxr.set((T)irow, i,1) ;
                indxc.set((T)icol, i,1);
                if(irow != icol)
                {
                    for(int l=1;l<=n;l++)	A.swap( irow, l, icol, l);
                    for(int l=1;l<=n;l++)	I.swap( irow, l, icol, l);
                    for(int l=1;l<=m;l++)	b.swap( irow, l, icol, l);
                }

                /*on peut maintenant faire les divisions et soustractions qui vont bien. */
                if(A.get(icol,icol) == (T)0)
                {
#ifdef verbose
                    cerr << "Singular Matrix : diag element : " << i << endl;
                    A.afficher();
                    //throw("Singular Matrix...");
#endif
                    A.set( numeric_limits<T>::epsilon(), icol, icol);
                }

                pivinv = (T)1.0/A.get(icol,icol);
                for(int l=1;l<=n;l++)	A.set(A.get(icol,l)*pivinv, icol, l);
                for(int l=1;l<=n;l++)	I.set(I.get(icol,l)*pivinv, icol, l);
                for(int l=1;l<=m;l++)	b.set(b.get(icol,l)*pivinv, icol, l);

                for(int ll=1;ll<=n;ll++)
                {
                    if(ll != icol)
                    {
                        T dum = A.get(ll, icol);
                        /*A.set((T)0, ll, icol);*/
                        /*il semble que cette opération soit bien faite
                         * dans la boucle for qui suit.... ?
                         */
                        for(int l=1;l<=n;l++)	A.set( A.get(ll,l)-dum*A.get(icol,l), ll, l);
                        for(int l=1;l<=n;l++)	I.set( I.get(ll,l)-dum*I.get(icol,l), ll, l);
                        for(int l=1;l<=m;l++)	b.set( b.get(ll,l)-dum*b.get(icol,l), ll, l);

                    }
                }

            }

		     Mat<T> Ai(I);
	
            /*------------------------------------*/

            *A_ = I;
            *b_ = b;


            return Ai;

            }
            break;







            case 2 :
            /* backsubstitution : triangular decomposition */

            {
            /* Full-Pivoting */
            int n, m;
            n = A.getLine();
            m = b.getColumn();
            T big, pivinv;

            /*bookkeeping */
            int irow, icol;
            Mat<T> indxc( n, 1);
            Mat<T> indxr( n,1);
            //Mat<T> ipiv((T)0, n,1);

            Mat<T> I( (T)0, n,n);
            for(int i=1;i<=n;i++)
                I.set((T)1, i,i);

            for(int i=1;i<=n;i++)
            {
                big = 0.0;
                for(int j=1;j<=n;j++)
                {
                    //if(ipiv.get(j,1) != (T)1)	/*si cette ligne j ne contient pas deja un pivot... */
                    //{
                        //for(int k=1;k<=n;k++)
                        //{
                            //if(ipiv.get(k,1) == 0) /* si cette ligne k ne contient pas deja un pivot..*/
                            //{
                                if(fabs_(A.get(i,j)) > big)
                                {
                                    big = fabs_(A.get(i,j));
                                    irow = i;
                                    icol = j;
                                }

                            //}
                        //}

                    //}
                }
                //ipiv.set(ipiv.get(icol,1)+1, icol,1);
                /*on a bien un pivot à la irow ligne, car c'est la valeur la plus grande de la colonne i*/

                /*Maintenant que l'on a le pivot, il faut le mettre sur la diagonale en interchangeant les lignes
                 * --> b doit changer de la meme manière que A ; l'inverse obtenue devra etre changee.
                 */

                indxr.set((T)irow, i,1) ;
                indxc.set((T)icol, i,1);
                if(irow != icol)
                {
                    for(int l=1;l<=n;l++)	A.swap( irow, l, icol,l);
                    for(int l=1;l<=n;l++)	I.swap( irow, l, icol,l);
                    for(int l=1;l<=m;l++)	b.swap( irow, l, icol,l);// attention a l'indice et sa limite.
                }
                /*echange des lignes et non des colonnes.*/
                /*du coup il faut faire le changement juste après de icol --> irow car le pivot est en irow,irow !*/

                /*on peut maintenant faire les divisions et soustractions qui vont bien. */
                if(A.get(irow,irow) == 0.0)
                {
                	//throw("Singular Matrix...");
                	A.set( numeric_limits<T>::epsilon(), irow, irow);
                }

                pivinv = 1.0/A.get(irow,irow);
                for(int l=1;l<=n;l++)	A.set(A.get(irow,l)*pivinv, irow, l);
                for(int l=1;l<=n;l++)	I.set(I.get(irow,l)*pivinv, irow, l);
                for(int l=1;l<=m;l++)	b.set(b.get(irow,l)*pivinv, irow, l);

                for(int ll=1;ll<=n;ll++)
                {
                    if(ll > irow)	//petit changement, != --> > de sorte que l'on ne change que les lignes en dessous de la position du pivot du moment.
                    {
                        T dum = A.get(ll, irow);
                        /*A.set((T)0, ll, icol);*/
                        /*il semble que cette opération soit bien faite
                         * dans la boucle for qui suit.... ?
                         */
                        for(int l=1;l<=n;l++)	A.set( A.get(ll,l)-dum*A.get(irow,l), ll, l);
                        for(int l=1;l<=n;l++)	I.set( I.get(ll,l)-dum*I.get(irow,l), ll, l);
                        for(int l=1;l<=m;l++)	b.set( b.get(ll,l)-dum*b.get(irow,l), ll, l);

                    }
                }
             
            }


            
		    Mat<T> Ai(A);
    
            /*------------------------------------*/
            /*------------------------------------*/
            /*------------------------------------*/
            /* BACKSUBSTITUTION */
            /*------------------------------------*/
            /*------------------------------------*/
            /*------------------------------------*/


            Mat<T> X(b);

            for(int j=1;j<=m;j++)
            {
                /*pour chaque colonne de b :*/
                for(int i=n;i>=1;i--)
                {
                    /*backsubstitution*/
                    T temp = 0;
                    T inv = 1.0/A.get(i,i);
                    for(int ii=i+1;ii<=n;ii++)
                        temp += A.get(ii,i)*X.get(ii,j);
                    X.set( inv*(b.get(i,j) - temp), i,j);
                }

            }

            *A_ = A;
            *b_ = b;


            return X;

            }
            break;



            default :
            throw("Mauvaise valeur de method.");
            break;
        }


    }
    else
        throw("Pointeurs arguments non-initializes.");

    return *A_;
}





template<typename T>
Mat<T> cramerInversion(const Mat<T>	& A)
{
	


}



#endif
