/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.1
                          Copyright (2021) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/** \file multiindex.cpp
 * \brief Tools that deal with integer multiindices.
 */

#include <math.h>
#include <assert.h>

#include "gen_defs.h"
#include "tools.h"
#include "arraytools.h"

// Computes the number of PC terms in a total-order expansion
int computeNPCTerms(int ndim,int norder)
{

  if (norder==-1) return 0;

  int enume=1;
  int denom=1;

  int minNO = min(norder,ndim);

  for(int k=0; k < minNO; k++){
      enume = enume*(norder+ndim-k);
      denom = denom*(k+1);
  }

  int nPCTerms = enume/denom;
  
  return nPCTerms;
}

// Computes total-order multiindex set, also return the number of elements
int computeMultiIndex(int ndim, int norder, Array2D<int> &mi)
{
  if ( ndim==0 ) return 1;

  // Compute the number of PC terms
  int npc=computeNPCTerms(ndim,norder);

  // Work arrays
  Array1D<int> ic(ndim,1);

  // Reset multi-index
  int iup  = 0;
  mi.Resize(npc,ndim,0);

  if (norder > 0) {

    //-----------first order terms---------------------------
    for(int idim=0; idim < ndim; idim++){
      iup++;
      mi(iup,idim) = 1; //multiIndex is a kronecker delta
    }

  } // done if order > 0

  if (norder > 1) {

    //-----------higher order terms--------------------------
    for(int iord=2; iord<norder+1; iord++){

      int lessiord = iup;                   //number of terms of order less than iord

      for(int idim=0;idim<ndim;idim++) 
        for(int ii=idim+1;ii<ndim;ii++) 
          ic(idim) += ic(ii);

      for(int idimm=0; idimm<ndim; idimm++) {

        for(int ii=lessiord-ic(idimm)+1; ii<lessiord+1; ii++){
          iup++;
          for(int idim=0; idim < ndim; idim++){
            mi(iup,idim) = mi(ii,idim);
          }
          mi(iup,idimm) += 1;
        }

      }

    } // done loop over orders

  } // done if order > 1

  return npc;
}

// Computes total-order multiindex set in an int* format; also return the number of elements
int computeMultiIndexT(int ndim, int norder, int *mi)
{
  /*
      Note that this function stores the multi-index in column-major
      format, i.e. mi[j*ndim+i] holds the j-th index for dimension i
   */
  
  if ( ndim==0 ) return 1;

  // Compute the number of PC terms
  int npc=computeNPCTerms(ndim,norder);

  // Work array
  int *ic = new int[ndim];
  for (size_t i=0; i < (size_t) ndim; i++ ) ic[i] = 1;

  // Reset multi-index
  int iup  = 0;
  for (size_t i=0; i < (size_t) npc*ndim; i++ ) mi[i]=0;

  if (norder > 0) {

    //-----------first order terms---------------------------
    for(int idim=0; idim < ndim; idim++){
      iup++;
      mi[iup*ndim+idim] = 1; //multiIndex is a kronecker delta
    }

  }

  if (norder > 1) {

    //-----------higher order terms--------------------------
    for(int iord=2; iord<norder+1; iord++){

      int lessiord = iup; //number of terms of order less than iord

      for(int idim=0;idim<ndim;idim++) 
        for(int ii=idim+1;ii<ndim;ii++) 
          ic[idim] += ic[ii];

      for(int idimm=0; idimm<ndim; idimm++) {

        for(int ii=lessiord-ic[idimm]+1; ii<lessiord+1; ii++){
          iup++;
          for(int idim=0; idim < ndim; idim++) {
            mi[iup*ndim+idim] = mi[ii*ndim+idim];
          }
          mi[iup*ndim+idimm] += 1;
        }

      }

    } // done loop over orders

  } // done if order > 1

  delete [] ic;

  return npc;
}

// Computes total-order multiindex with custom ordering
int computeMultiIndex(int ndim, int norder, Array2D<int> &mi, string ordtype)
{

  if ( ndim==0 ) return 1;

  // Test ordtype is known
  bool isLex        = ( ordtype == string("lex")      ) || ( ordtype == ("LEX")      );
  bool isRevLex     = ( ordtype == string("revlex")   ) || ( ordtype == ("REVLEX")   );
  bool isCoLex      = ( ordtype == string("colex")    ) || ( ordtype == ("COLEX")    );
  bool isRevCoLex   = ( ordtype == string("revcolex") ) || ( ordtype == ("REVCOLEX") );
  bool isLexAll     = ( ordtype == string("lexall")      ) || ( ordtype == ("LEXALL")      );
  bool isRevLexAll  = ( ordtype == string("revlexall")   ) || ( ordtype == ("REVLEXALL")   );
  bool isCoLexAll   = ( ordtype == string("colexall")    ) || ( ordtype == ("COLEXALL")    );
  bool isRevCoLexAll= ( ordtype == string("revcolexall") ) || ( ordtype == ("REVCOLEXALL") );

  if ( not isLex && not isRevLex && not isCoLex  && not isRevCoLex 
       && not isLexAll && not isRevLexAll && not isCoLexAll  && not isRevCoLexAll ) {
    throw Tantrum(string("computeMultiIndex: The passed in ordtype \"")
           + ordtype + string("\" does not match any available option"));  
   }

  // Compute the number of PC terms
  int npc=computeNPCTerms(ndim,norder);

  // Work arrays
  Array1D<int> ic(ndim,1);

  // Reset multi-index
  int iup  = 0;
  mi.Resize(ndim,npc,0);

  int npcTmp = computeMultiIndexT(ndim, norder, mi.GetArrayPointer());

  if ( norder == 0 ) {
    Array2D<int> miT;
    transpose(mi, miT);
    mi=miT;
    return npc;
  }  

  if ( isLexAll || isRevLexAll || isCoLexAll || isRevCoLexAll ) {
    int isgn=0, i1=0, j1=0, index=0;
    int iordST=1, iordEN = npc; 
    int nTerms=iordEN-iordST;
    do {

      heap_ext_(&nTerms, &isgn, &i1, &j1, &index);
      if (index < 0) {
        isgn = 0;
        if ( isLexAll ) {
          for ( int j = 0; j < ndim; j++ ) {
            if ( mi(j,iordST+i1-1) < mi(j,iordST+j1-1) ) {
              isgn = -1; break;
            }
            else if ( mi(j,iordST+j1-1) < mi(j,iordST+i1-1) ) {
              isgn = 1;
              break;
            }
          }
        } // end if ordtype = lex
        if ( isRevLexAll ) {
          for ( int j = 0; j < ndim; j++ ) {
            if ( mi(j,iordST+i1-1) > mi(j,iordST+j1-1) ) {
              isgn = -1; break;
            }
            else if ( mi(j,iordST+j1-1) > mi(j,iordST+i1-1) ) {
              isgn = 1;
              break;
            }
          }
        } // end if ordtype = revlex
        else if ( isCoLexAll ) {
          for ( int j = ndim-1; j >= 0; j-- ) {
            if ( mi(j,iordST+i1-1) < mi(j,iordST+j1-1) ) {
              isgn = -1; break;
            }
            else if ( mi(j,iordST+j1-1) < mi(j,iordST+i1-1) ) {
              isgn = 1;
              break;
            }
          }
        } // end if ordtype = colex
        else if ( isRevCoLexAll ) {
          for ( int j = ndim-1; j >= 0; j-- ) {
            if ( mi(j,iordST+i1-1) > mi(j,iordST+j1-1) ) {
              isgn = -1; break;
            }
            else if ( mi(j,iordST+j1-1) > mi(j,iordST+i1-1) ) {
              isgn = 1;
              break;
            }
          }
        } // end if ordtype = revcolex
      }

      if (index > 0) {
        int itmp ;
        for ( int j = 0; j < ndim; j++ ) {
          itmp = mi(j,iordST+i1-1);
          mi(j,iordST+i1-1) = mi(j,iordST+j1-1);
          mi(j,iordST+j1-1) = itmp;
        }
      }

    } while (index != 0);   
 
  } /* done if ( isLexAll || isRevLexAll || isCoLexAll || isRevCoLexAll ) */

  else {

    for (int iord=1; iord<norder+1; iord++){

      int iordST, iordEN=npc;
      for (int i=1; i<npc; i++) {
        if (mi(0,i) == iord) {
          iordST=i;
          break;
	}
      }
      for (int i=1; i<npc; i++) {
        if (mi(0,i) == iord+1) {
          iordEN=i;
          break;
	}
      }

      int nTerms=iordEN-iordST, isgn=0, i1=0, j1=0, index=0;
      do {

        heap_ext_(&nTerms, &isgn, &i1, &j1, &index);
        if (index < 0) {
          isgn = 0;
          if ( ( ordtype == string("lex") ) || ( ordtype == ("LEX") ) ) {
            for ( int j = 0; j < ndim; j++ ) {
              if ( mi(j,iordST+i1-1) < mi(j,iordST+j1-1) ) {
                isgn = -1; break;
              }
              else if ( mi(j,iordST+j1-1) < mi(j,iordST+i1-1) ) {
                isgn = 1;
                break;
              }
            }
          } // end if ordtype = lex
          if ( ( ordtype == string("revlex") ) || ( ordtype == ("REVLEX") ) ) {
            for ( int j = 0; j < ndim; j++ ) {
              if ( mi(j,iordST+i1-1) > mi(j,iordST+j1-1) ) {
                isgn = -1; break;
              }
              else if ( mi(j,iordST+j1-1) > mi(j,iordST+i1-1) ) {
                isgn = 1;
                break;
              }
            }
          } // end if ordtype = revlex
          else if ( ( ordtype == string("colex") ) || ( ordtype == ("COLEX") ) ) {
            for ( int j = ndim-1; j >= 0; j-- ) {
              if ( mi(j,iordST+i1-1) < mi(j,iordST+j1-1) ) {
                isgn = -1; break;
              }
              else if ( mi(j,iordST+j1-1) < mi(j,iordST+i1-1) ) {
                isgn = 1;
                break;
              }
            }
          } // end if ordtype = colex
          else if ( ( ordtype == string("revcolex") ) || ( ordtype == ("REVCOLEX") ) ) {
            for ( int j = ndim-1; j >= 0; j-- ) {
              if ( mi(j,iordST+i1-1) > mi(j,iordST+j1-1) ) {
                isgn = -1; break;
              }
              else if ( mi(j,iordST+j1-1) > mi(j,iordST+i1-1) ) {
                isgn = 1;
                break;
              }
            }
          } // end if ordtype = revcolex
	}

        if (index > 0) {
          int itmp ;
          for ( int j = 0; j < ndim; j++ ) {
            itmp = mi(j,iordST+i1-1);
            mi(j,iordST+i1-1) = mi(j,iordST+j1-1);
            mi(j,iordST+j1-1) = itmp;
          }
        }
      } while (index != 0);

    } /* end loop over orders */

  } /* end else for order by order re-arrangement */

  /* return */
  Array2D<int> miT;
  transpose(mi, miT);
  mi=miT;

  return npc;

}

// Computes tensor-product multiindex set, also return the number of elements
int computeMultiIndexTP(Array1D<int>& maxorders, Array2D<int>& mindex)
{
  int ndim=maxorders.XSize();
  int npc;


  if (ndim==1){
    mindex.Resize(maxorders(0)+1,1,0);
    for(int i=0;i<=maxorders(0);i++)
      mindex(i,0)=i;

    return maxorders(0)+1;
  }

  Array1D<int> maxorders_1(ndim-1,0);
  for (int j=0;j<ndim-1;j++)
    maxorders_1(j)=maxorders(j);
  Array2D<int> mindex_1;
  int npc_1=computeMultiIndexTP(maxorders_1,mindex_1);

  npc=npc_1*(maxorders(ndim-1)+1);
  mindex.Resize(npc,ndim,-1);

  for (int ib=0;ib<=maxorders(ndim-1);ib++){
    for (int k=0; k<npc_1; k++){
      for (int j=0; j< ndim-1; j++){
        mindex(ib*npc_1+k,j)=mindex_1(k,j);
        mindex(ib*npc_1+k,ndim-1)=ib;
      }
    }
  }


  return npc;
}

// Computes the number of terms in an HDMR multiindex set
int computeNPCTermsHDMR(int ndim,  Array1D<int>& maxorders)
{
  int nhdmr=maxorders.XSize()-1;

  int cnt=0;
  for(int i=0;i<=nhdmr;i++){
    cnt+= ( choose(maxorders(i),i)*choose(ndim,i) );
  }

  return cnt;
}

// Computes HDMR multiindex set, also return the number of elements
int computeMultiIndexHDMR(int ndim,Array1D<int>& maxorders, Array2D<int>& mindex)
{
  int nhdmr=maxorders.XSize()-1;
  int npc=computeNPCTermsHDMR(ndim,maxorders);

  mindex.Resize(npc,ndim,0);
  int iup=1;
  for(int i=1; i<=nhdmr;i++){
    Array2D<int> ind;
    chooseComb(ndim,i,ind);
   
    if (maxorders(i)<i)
      return iup;
    else{
      Array2D<int> mi;
      computeMultiIndex(i,maxorders(i)-i,mi);
   
      for(int j=0;j<(int)ind.XSize();j++){
        for(int k=0; k<(int)mi.XSize();k++){
          for(int ii=0;ii<i;ii++){
            mindex(iup,ind(j,ii))=mi(k,ii)+1;
          }
          iup++;
        }
      }
    }
  }

  CHECKEQ(iup,npc);
  return npc; //or iup;
}

// Given multtindex in a sparse format, recover it in the dense format
void decodeMindex(Array1D< Array2D<int> >& sp_mindex, int ndim, Array2D<int>& mindex){
    
    int npc=0;
    for(int i=0;i<(int) sp_mindex.XSize();i++){
        npc+=sp_mindex(i).XSize();
    }
    mindex.Resize(npc,ndim,0);
    int ipc=0;
    for(int i=0;i<(int) sp_mindex.XSize();i++){
      for(int j=0;j<(int) sp_mindex(i).XSize();j++){
        for(int i_effdim=0;i_effdim<i;i_effdim++){
           mindex(ipc,sp_mindex(i)(j,i_effdim))=sp_mindex(i)(j,i_effdim+i);
        }
        ipc++;
      }
    }

    CHECKEQ(ipc,npc);
    
    return;
}


// Checking if a basis is admissible with respect to a given multiindex set
bool is_admis(Array1D<int>& mindex_try,Array2D<int>& mindex){
    
  bool admis=true;
  int npc=mindex.XSize();
  int ndim=mindex.YSize();
  
  Array1D<int> tmp;
  tmp=mindex_try;
  
  for(int j=0;j<ndim;j++){
      
    if(mindex_try(j)>0){
        tmp(j)--;
        for(int ipc=0;ipc<npc;ipc++){
            
          Array1D<int> cur_mindex;
          getRow(mindex,ipc,cur_mindex);
          
          if(is_equal(tmp,cur_mindex)){
              admis=true;
              break;
          }
          admis=false;
        
        }
        
        if(admis==false)
            break;
        
        tmp(j)++;
    }
      
  }
  
  return admis;
}

// Increase a multiindex set by one order with admissible bases
void upOrder(Array2D<int>& mindex,Array2D<int>& new_mindex){
    
  int npc=mindex.XSize();
  int ndim=mindex.YSize();
    
  Array1D<int> orders;
  getOrders(mindex,orders);

  new_mindex=mindex;
  
  int imax;
  int maxOrd=maxVal(orders,&imax);
  
  for(int ipc=0;ipc<npc;ipc++){
    if (orders(ipc)==maxOrd){

      int nzind=0;
      while (mindex(ipc,nzind)==0 && nzind<ndim-1) {
        nzind++;
      }
          
      Array1D<int> new_mindex_try;
      getRow(mindex,ipc,new_mindex_try);
      for(int j=0;j<=nzind;j++){
        new_mindex_try(j)++;
        if(is_admis(new_mindex_try,new_mindex)){
          paddMatRow(new_mindex,new_mindex_try);              
        }
          
        new_mindex_try(j)--;
      } 
    }
  }
  
  return;
}

// Gets the total degree of each multiindex term
void getOrders(Array2D<int>& mindex,Array1D<int>& orders){
    
  int npc=mindex.XSize();
  int ndim=mindex.YSize();

  orders.Resize(npc,0);
    
  for (int ipc=0; ipc<npc ; ipc++) {
      int sum=0;
      for(int id=0;id<ndim;id++)
          sum+=mindex(ipc,id);
      orders(ipc)=sum;
  }

  return;

}

// Given a single multiindex, this returns its relative position in the total-order multiindex set
int get_invmindex(Array1D<int> mi)
{

  int nd=mi.XSize();
  int ss=0;

  for(int id=0;id<nd;id++)
    ss+=mi(id);
      
  int index=computeNPCTerms(nd,ss-1);
  index += get_invmindex_ord(mi);

  return index;

}

// Given a single multiindex, this returns its relative position in the total-order multiindex set among the bases of the same order
int get_invmindex_ord(Array1D<int> mi)
{
  int index=0;
  int nd=mi.XSize();

  int ss=0;
  for(int id=0;id<nd;id++)
    ss+=mi(id);

  if (nd>1){

    for(int ii=ss;ii>mi(0);ii--)
      index+=computeNPCTerms(nd-2,ss-ii);

    Array1D<int> mic(nd-1,0);
    for(int id=0;id<nd-1;id++)
      mic(id)=mi(id+1);
    
    index+=get_invmindex_ord(mic);

  }


  return index;

}
